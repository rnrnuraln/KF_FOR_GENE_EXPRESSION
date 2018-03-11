package common

import java.nio.file.{Files, Paths}
import utils.Utils
import breeze.linalg.DenseVector

/**
  * Optimize parameters so that the model fits data, using
  * EMアルゴリズムを用いて最適化する
  * 交差検証による隠れ変数の次元の決定、など
  * Created by tara on 2017/07/08.
  */
class Optimize(seqs: List[Array[ConObs]], cond: OptimizeCond) {

  val m: Int = {
    val s1 = seqs.head
    def getM(i: Int): Int = {
      s1(i).observe match {
        case Some(b) => b.length
        case None => getM(i + 1)
      }
    }
    getM(0)
  }
  val l: Int = {
    seqs.head(0).control match {
      case Some(c) => c.length
      case None => 0
    }
  }

  private val emCondCommon = EMCond(emTime = cond.emTime(0), delta = cond.delta(0),
    AOpt = cond.AOpt, BOpt = cond.BOpt, HOpt = cond.HOpt, QOpt = cond.QOpt,
    ROpt = cond.ROpt, initStateMeanOpt = cond.initStateMeanOpt,
    initStateCovarianceOpt = cond.initStateCovarianceOpt,
    initkf = null, initStateMeans = List(), initStateCovariances = List())

  /**
    * how to
    * EMCondを元に、どんなlearnが必要かを求める。
    *
    * @return
    */
  def learn(): EMOutputs = {
    cond.emHid._1 match {
      case "synchro" => learnSeveralTimes(cond.emHid._2(0), seqs)
      case "crossValid" => hidDimUnknownLearn()
      case "fixed" => hidDimCandidateLearn(cond.emHid._2)
      case _ => EMOutputs(null, List(), null, null)
    }
  }

  /**
    * 交差検証を用いてlogLikelihoodを擬似的に求める: EMOutputにはここのlikelihoodを入れる。
    * 予測結果をpredictFileに書き込む
    *
    * @param n
    * @return
    */
  def crossValid(n: Int): (EMOutputs, String) = {
    val sb = new StringBuilder()
    val seqNum = seqs.length
    val seqsPerFold = seqNum / cond.foldNum(0)
    val allLikelihood = (0 until cond.foldNum(0)).foldLeft(0.0) {(s, i) =>
      val from = i * seqsPerFold
      val until = if (i == cond.foldNum(0) - 1) seqNum else (i+1) * seqsPerFold
      val validateSeqs = seqs.slice(from, until)
      val seqsWithoutOne = seqs.take(from) ::: seqs.drop(until)
      val learnedEMOutputs = learnSeveralTimes(n, seqsWithoutOne)
      val fws = validateSeqs.map(seq => Forward(learnedEMOutputs.kf, seq, learnedEMOutputs.initStateMeanMean, learnedEMOutputs.initStateCovarianceMean))
      if (cond.crossValidPredictFile != "") {
        val predicts: Seq[Array[DenseVector[Double]]] = fws.map(x => x.predicts)
        predicts.foreach(array => {
          array.foreach(vector => {
            sb.append(utils.Utils.vectorTostring(vector) + "\n")
          })
        })
      }
      s + fws.foldLeft(0.0) {(s2, x) => s2 + x.logLikelihood}
    }
    val allEMOutputs = learnSeveralTimes(n, seqs)
    val allEMOutputsFixed = allEMOutputs.copy(logLikelihoods = List(allLikelihood))
    (allEMOutputsFixed, sb.toString())
  }

  /**
    * 隠れ状態の次元を交差検証を用いて求めるタイプ
    * 1から順番に求めていく奴だけど…… 二部探索とかのが良さげ?
    */
  def hidDimUnknownLearn(): EMOutputs = {
    def hidDimUnknownLearnSub(n: Int, emOutputs: EMOutputs, predict: String): (EMOutputs, String) = {
      if (n > m)
        (emOutputs, predict)
      else {
        val c = crossValid(n)
        if (emOutputs.logLikelihoods.head > c._1.logLikelihoods.head) (emOutputs, predict) else hidDimUnknownLearnSub(n + 1, c._1, c._2)
      }
    }
    val crossValidated = hidDimUnknownLearnSub(1, EMOutputs(logLikelihoods = List(Double.MinValue)), "")
    if (cond.crossValidPredictFile != "") {
      val fp = utils.Utils.makePrintWriter(cond.crossValidPredictFile)
      fp.print(crossValidated._2)
      fp.close()
    }
    crossValidated._1
  }

  /**
    * 隠れ状態の次元の候補から求めるタイプ
    */
  def hidDimCandidateLearn(hids: Array[Int]): EMOutputs = {
    val crossValidated = hids.map(crossValid(_)).max
    if (cond.crossValidPredictFile != "") {
      val fp = utils.Utils.makePrintWriter(cond.crossValidPredictFile)
      fp.print(crossValidated._2)
      fp.close()
    }
    crossValidated._1
  }

  /**
    * In order to handle initial value
    *
    * @return
    */
  def learnSeveralTimes(n: Int, seqs: List[Array[ConObs]]): EMOutputs = {
    val emRandPar = (0 until cond.emRand(0)).par
    emRandPar.tasksupport = new scala.collection.parallel.ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(cond.parallelThreadNum(0)))
    val maxEMOutputs = emRandPar.map(x => {
      val kf = KalmanFilter(cond.Ainit.copy().makeMatrixForA(n), cond.Binit.copy().makeMatrixForB(n, l), cond.Hinit.copy().makeMatrixForH(m, n), cond.Qinit.copy().makeMatrixForCov(n), cond.Rinit.copy().makeMatrixForCov(m))
      val initStateMeans = cond.initStateMeanInit.copy().makeInitStateMeanList(n, seqs.length)
      val initStateCovariances = cond.initStateCovarianceInit.copy().makeInitStateCovarianceList(n, seqs.length)
      val emCond = emCondCommon.copy(emTime = cond.emShallow(0), initkf = kf, initStateMeans = initStateMeans, initStateCovariances = initStateCovariances)
      val emAlgorithm = EMAlgorithm(seqs, emCond)
      val learned = emAlgorithm.learnBasis()
      learned
    }).max
    val emCond = emCondCommon.copy(emTime = cond.emTime(0), initkf = maxEMOutputs.kf, initStateMeans = maxEMOutputs.initStateMeans, initStateCovariances = maxEMOutputs.initStateCovariance)
    val emAlgorithm = EMAlgorithm(seqs, emCond)
    println("log: \n" + emAlgorithm.learnBasis())
    emAlgorithm.learnBasis()
  }

}

object Optimize {
  def learn(inputSeqs: String, output: String, condition: String): Unit = {
    val basedir = Paths.get(output)
    if (Files.notExists(basedir)) Files.createFile(basedir)
    val seqs = Utils.makeSeqs(inputSeqs)
    val optCond = OptimizeCond.apply(condition)
    val optimize = new Optimize(seqs, optCond)
    val outputs = optimize.learn()
    val fp = Utils.makePrintWriter(output)
    fp.print(outputs.toString())
    fp.close()
  }
}