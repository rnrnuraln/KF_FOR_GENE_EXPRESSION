package common

import breeze.linalg._
import java.nio.file.{Paths, Files}

/**
  * Created by tara on 2017/06/08.
  */
class Evaluate() {

}

object Evaluate {

  /**
    *
    * @param input  input file
    * @param output output file
    */
  def run(input: String, output: String): Unit = {
    val basedir = Paths.get(output)
    if (Files.notExists(basedir)) Files.createFile(basedir)
    val params: ExperimentCond = ExperimentCond.apply(input)
    val optimizeParamsCommon = OptimizeCond.apply(input)

    def makeControlWithObservePoint(seqNum: Int, seqLen: Int, l: Int): List[Array[(Option[DenseVector[Double]], Boolean)]] = {
      (0 until seqNum).map(x => {
        val control = params.controlGen.makeControlSeq(l, seqLen)
        val observe = (0 until seqLen).map(i => {
          params.observeLocation._1 match {
            case "all" => true
            case "some" =>
              val loc = params.observeLocation._2
              i % (loc(0) + loc(1)) >= loc(0)
            case _ => true
          }
        })
        control.zip(observe)
      }).toList
    }

    val sb = new StringBuilder()
    sb.append(ExperimentOutputCond.apply().classNameToString())
    val maxSeqNum = params.seqNum.max
    val paramStringArray: Array[String] =
      for {n <- params.hDim
           oDims: Array[Int] = if (params.oDim._1 == "synchro") Array(n) else params.oDim._2
           newEMHid = if (params.emHid._1 == "synchro") (params.emHid._1, Array(n)) else params.emHid
           m <- oDims
           l <- params.cDim
           seqLen <- params.seqLen
           testSeqNum <- params.testSeqNum
           a = params.Agen.makeMatrixForA(n)
           b = params.Bgen.makeMatrixForB(n, l)
           h = params.Hgen.makeMatrixForH(m, n)
           q = params.Qgen.makeMatrixForCov(n)
           r = params.Rgen.makeMatrixForCov(m)
           kf = KalmanFilter(a, b, h, q, r)
           aInit = optimizeParamsCommon.Ainit.copy(m = a)
           bInit = optimizeParamsCommon.Binit.copy(m = b.getOrElse(null))
           hInit = optimizeParamsCommon.Hinit.copy(m = h)
           qInit = optimizeParamsCommon.Qinit.copy(m = q)
           rInit = optimizeParamsCommon.Rinit.copy(m = r)
           initStateMeans = params.initStateMeanGen.makeInitStateMeanList(n, maxSeqNum)
           initStateCovariances = params.initStateCovarianceGen.makeInitStateCovarianceList(n, maxSeqNum)
           initInitStateMean = optimizeParamsCommon.initStateMeanInit.copy(v = initStateMeans(0))
           initInitStateCovariance = optimizeParamsCommon.initStateCovarianceInit.copy(m = initStateCovariances(0))
           dummy = println("true: " + EMOutputs(kf, List(0.0), initStateMeans, initStateCovariances))
           conObsBasis = makeControlWithObservePoint(maxSeqNum, seqLen, l)
           kfRuns: List[Array[(DenseVector[Double], Option[DenseVector[Double]])]] = conObsBasis.zip(initStateMeans).zip(initStateCovariances).map(x => kf.run(x._1._1, x._1._2, x._2))
           seqs = conObsBasis.zip(kfRuns).map(x => x._1.zip(x._2).map(y => ConObs(y._1._1, y._2._2)))
           //dummy2 = seqs.foreach(x => x.foreach(y => println({y.control match {case Some(c) => utils.Utils.vectorTostring(c); case None => ""}} + "\t" +
           //  {y.observe match {case Some(c) => utils.Utils.vectorTostring(c); case None => ""}})))
           //テスト用の系列を生成
           conObsBasis2 = makeControlWithObservePoint(seqLen, testSeqNum, l)
           kfTest = conObsBasis2.zip(initStateMeans).zip(initStateCovariances).map(x => kf.run(x._1._1, x._1._2, x._2))
           seqs2 = conObsBasis2.zip(kfTest).map(x => x._1.zip(x._2).map(y => ConObs(y._1._1, y._2._2)))
           seqNum <- params.seqNum
           emTime <- params.emTime
           emShallow <- params.emShallow
           emRand <- params.emRand
           delta <- params.delta
           optimizeCond = optimizeParamsCommon.copy(emHid = newEMHid, emTime = Array(emTime), emShallow = Array(emShallow), emRand = Array(emRand), delta = Array(delta),
             Ainit = aInit, Binit = bInit, Hinit = hInit, Qinit = qInit, Rinit = rInit, initStateMeanInit = initInitStateMean, initStateCovarianceInit = initInitStateCovariance)
           //optimization
           start = System.currentTimeMillis()
           optimized = new Optimize(seqs.take(seqNum), optimizeCond).learn()
           end = System.currentTimeMillis()

           /**
             * Evaluationは、尤度、予測、真の値と推定値とのズレなどによって行う
             * 分散とかも考慮したevaluationもやりたい。
             */
           calcTime = end - start
           optimizeN = optimized.kf.n
           likelipreds = seqs2.foldLeft(0.0, 0.0) {(s, seq) =>
             val fw = Forward(optimized.kf, seq, optimized.initStateMeanMean, optimized.initStateCovarianceMean)
             (s._1 + fw.logLikelihood, s._2 + fw.predCorrEval)
           }
           //真の値でのlikelihood
           realLikepreds = seqs2.zip(initStateMeans.zip(initStateCovariances)).foldLeft(0.0, 0.0) { (s, seq) =>
             val fw = Forward(kf, seq._1, seq._2._1, seq._2._2)
             (s._1 + fw.logLikelihood, s._2 + fw.predCorrEval)
           }

           fixDelta <- params.fixDelta
           kfDif = if (optimizeN != n) Array(-1.0, -1.0, -1.0, -1.0, -1.0) else kf.maxDifs(optimized.kf, fixDelta)
           initStateMeanDif = if (optimizeN != n) -1.0 else initStateMeans.zip(optimized.initStateMeans).map(x => utils.Utils.maxDif(x._1.toArray, x._2.toArray, fixDelta)).max
           initStateCovarianceDif = if (optimizeN != n) -1.0 else initStateCovariances.zip(optimized.initStateCovariance).map(x => x._1.maxDif(x._2, fixDelta)).max
           outputParams = ExperimentOutputCond(
             hDim = n, cDim = l, oDim = m, seqNum = seqNum, seqLen = seqLen, emHidCond = params.emHid._1, emHidDim = optimizeN, foldNum = params.foldNum(0),
             emTime = emTime, emShallow = emShallow, emRand = emRand, delta = delta, fixDelta = fixDelta, Aerror = kfDif(0), Berror = kfDif(1), Herror = kfDif(2), Qerror = kfDif(3), Rerror = kfDif(4),
             initStateMeanError = initStateMeanDif, initStateCovarianceError = initStateCovarianceDif, calcTime = calcTime, logLikelihood = likelipreds._1, predCorr = likelipreds._2,
             realLikelihood = realLikepreds._1, realCorrPred = realLikepreds._2, testSeqNum = testSeqNum, parallelThreadNum = params.parallelThreadNum(0), shFile = params.shFile, jobId = params.jobId, version = params.version)
      } yield outputParams.toString()
    sb.append(paramStringArray.reduce(_ + _))

    val fp = utils.Utils.makePrintWriter(output)
    fp.print(sb.toString())
    fp.close()
  }

}