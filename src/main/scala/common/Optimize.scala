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
    initkf = null, initStateMeans = List(), initStateCovariances = List(), seqRegularization = cond.seqRegularization,
    showLikelihood = cond.showLikelihood, RRegularization = cond.seqRegularization
  )

  /**
    * how to
    * EMCondを元に、どんなlearnが必要かを求める。
    *
    * @return
    */
  def learn(): EMOutputs = {
    val learned = cond.AOpt.regParams.foldLeft((EMOutputs(), ""), Array(-1.0)) { (s, x) =>
      cond.AOpt.lambda = x
      val ema = cond.BOpt.regParams.foldLeft((EMOutputs(), ""), Array(-1.0)) { (s, y) =>
        cond.BOpt.lambda = y
        val emb = cond.HOpt.regParams.foldLeft((EMOutputs(), ""), Array(-1.0)) { (s, z) =>
          cond.HOpt.lambda = z
          println("lambda A: " + x + " B: " + y + " H: " + z)
          val emh = learnSub()
          println("CrossValid likelihood: " + emh._1.logLikelihoods(0))
          if (emh._1 > s._1._1) (emh, Array(z)) else s
        }
        if (emb._1._1 > s._1._1) (emb._1, Array(y, emb._2(0))) else s
      }
      if (ema._1._1 > s._1._1) (ema._1, Array(x, ema._2(0), ema._2(1))) else s
    }
    println("Max Lambda A: " + learned._2(0) + "B: " + learned._2(1) + "H: " + learned._2(2))
    if (cond.crossValidPredictFile != "") {
      val fp = utils.Utils.makePrintWriter(cond.crossValidPredictFile)
      fp.print(learned._1._2)
      fp.close()
    }
    learned._1._1
  }

  def learnSub(): (EMOutputs, String) = {
    cond.emHid._1 match {
      case "synchro" => learnSeveralTimes(cond.emHid._2(0), seqs)
      case "crossValid" => hidDimUnknownLearn()
      case "fixed" => hidDimCandidateLearn(cond.emHid._2)
      case _ => (EMOutputs(null, List(), null, null), "")
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
    val allLikelihood = (0 until cond.foldNum(0)).foldLeft(0.0) {(s, i) =>
      val from = seqNum * i / cond.foldNum(0)
      val until = seqNum * (i + 1) / cond.foldNum(0)
      val validateSeqs = seqs.slice(from, until)
      val seqsWithoutOne = seqs.take(from) ::: seqs.drop(until)
      val learnedEMOutputs = learnSeveralTimes(n, seqsWithoutOne)
      val regularizedValidateSeqs = validateSeqs.map(x => {
        x.map(y => {
          y.observe match {
            case Some(o) => ConObs(y.control, Some((o - learnedEMOutputs._1.allMean) /:/ learnedEMOutputs._1.allStd))
            case None => ConObs(y.control, y.observe)
          }
        })
      })
      val fws = regularizedValidateSeqs.map(seq => Forward(learnedEMOutputs._1.kf, seq, learnedEMOutputs._1.initStateMeanMean, learnedEMOutputs._1.initStateCovarianceMean))
      if (cond.crossValidPredictFile != "") {
        val predicts: Seq[Array[DenseVector[Double]]] = fws.map(x => x.predicts)
        predicts.foreach(array => {
          array.foreach(vector => {
            val regularizedVector = vector *:* learnedEMOutputs._1.allStd + learnedEMOutputs._1.allMean
            sb.append(utils.Utils.vectorTostring(regularizedVector) + "\n")
          })
        })
      }
      s + fws.foldLeft(0.0) {(s2, x) => s2 + x.logLikelihood}
    }
    val allEMOutputs = learnSeveralTimes(n, seqs)
    val allEMOutputsFixed = allEMOutputs.copy(_1 = allEMOutputs._1.copy(logLikelihoods = List(allLikelihood)), _2 = sb.toString())
    allEMOutputsFixed
  }

  /**
    * 隠れ状態の次元を交差検証を用いて求めるタイプ
    * 1から順番に求めていく奴だけど…… 二部探索とかのが良さげ?
    */
  def hidDimUnknownLearn(): (EMOutputs, String) = {
    def hidDimUnknownLearnSub(n: Int, emOutputs: EMOutputs, predict: String): (EMOutputs, String) = {
      if (n > m)
        (emOutputs, predict)
      else {
        val c = crossValid(n)
        if (emOutputs.logLikelihoods.head > c._1.logLikelihoods.head) (emOutputs, predict) else hidDimUnknownLearnSub(n + 1, c._1, c._2)
      }
    }
    val crossValidated = hidDimUnknownLearnSub(1, EMOutputs(logLikelihoods = List(Double.MinValue)), "")
    crossValidated
  }

  /**
    * 隠れ状態の次元の候補から求めるタイプ
    */
  def hidDimCandidateLearn(hids: Array[Int]): (EMOutputs, String) = {
    hids.map(crossValid(_)).max
  }

  /**
    * In order to handle initial value
    *
    * @return
    */
  def learnSeveralTimes(n: Int, seqs: List[Array[ConObs]]): (EMOutputs, String) = {
    def learnSeveralTimesSub(iterNum: Int, emOutputList: List[EMOutputs]): (EMOutputs, String) = {
      //処理する
      val randNum: Int = (cond.emRand(0) * Math.pow(2, iterNum - 1)).toInt
      val emRandPar = (0 until randNum).par
      emRandPar.tasksupport = new scala.collection.parallel.ForkJoinTaskSupport(new scala.concurrent.forkjoin.ForkJoinPool(cond.parallelThreadNum(0)))
      val emNum = if (randNum < 5) randNum else 5
      val emOutputs = emRandPar.map(x => {
        val kf = KalmanFilter(cond.Ainit.copy().makeMatrixForA(n), cond.Binit.copy().makeMatrixForB(n, l), cond.Hinit.copy().makeMatrixForH(m, n), cond.Qinit.copy().makeMatrixForCov(n), cond.Rinit.copy().makeMatrixForCov(m))
        val initStateMeans = cond.initStateMeanInit.copy().makeInitStateMeanList(n, seqs.length)
        val initStateCovariances = cond.initStateCovarianceInit.copy().makeInitStateCovarianceList(n, seqs.length)
        val emCond = emCondCommon.copy(emTime = cond.emShallow(0), initkf = kf, initStateMeans = initStateMeans, initStateCovariances = initStateCovariances)
        val emAlgorithm = EMAlgorithm(seqs, emCond)
        val learned = emAlgorithm.learnBasis()
        learned
      }).toList.sorted.reverse.take(emNum)
      val newEmOutputs = (emOutputs ::: emOutputList).sorted.reverse.take(emNum)
      val meanNum = 3 //meanNum個分の平均をとって、topのloglikelihoodとの差分を比較、それが十分小さければ計算をやめる
      val maxLikelihood = newEmOutputs(0).logLikelihoods(0)
      val mean = (newEmOutputs.take(meanNum + 1).map(x => x.logLikelihoods(0)).sum - maxLikelihood) / meanNum
      if (iterNum >= cond.emRand(1) || (mean - maxLikelihood) / maxLikelihood < cond.emRand(2)) {
        println(newEmOutputs.foldLeft("RandLogLikelihoods:") {(s, x) => s + "\t" + x.logLikelihoods(0) })
        println("final RandNum: " + (randNum * 2 - cond.emRand(0)).toInt)
        val maxEMOutputs = newEmOutputs(0)
        val emCond = emCondCommon.copy(emTime = cond.emTime(0), initkf = maxEMOutputs.kf, initStateMeans = maxEMOutputs.initStateMeans, initStateCovariances = maxEMOutputs.initStateCovariance, RRegularization = false)
        val emAlgorithm = EMAlgorithm(seqs, emCond)
        val learned = emAlgorithm.learnBasis()
        val printLoglikelihoodNum = if (learned.logLikelihoods.length > 21) 21 else learned.logLikelihoods.length
        val logLikelihoodLists = learned.logLikelihoods.take(printLoglikelihoodNum).zipWithIndex.foldLeft("") {(s, x) => if (x._2 % 3 == 0) s + "\t" + x._1 else s}
        println("likelihood lists: " + logLikelihoodLists)
        println("log: \n" + learned)
        (emAlgorithm.learnBasis(), "")
      } else {
        learnSeveralTimesSub(iterNum + 1, newEmOutputs)
      }
    }
    learnSeveralTimesSub(1, List())
    //log likelihoodを順番に出力してみる
  }

}

object Optimize {
  def learn(inputSeqs: String, output: String, condition: String): Unit = {
    val start = System.currentTimeMillis()
    val basedir = Paths.get(output)
    if (Files.notExists(basedir)) Files.createFile(basedir)
    val seqs = Utils.readSeqs(inputSeqs)
    val optCond = OptimizeCond.apply(condition)
    val optimize = new Optimize(seqs, optCond)
    val outputs = optimize.learn()
    val fp = Utils.makePrintWriter(output)
    fp.print(outputs.toString())
    fp.close()
    val end = System.currentTimeMillis()
    print((end - start) + ":ms for optimization")
  }
}