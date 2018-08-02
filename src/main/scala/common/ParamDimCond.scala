package common

import breeze.linalg._

import scala.collection.immutable.Map
import scala.io.Source

/**
  * Created by tara on 2016/12/21.
  * Obtain parameters that are used in experiments and in a map.
  *
  * 実験などに使うパラメータたちをファイルから取得して、Mapの形で格納。
  * 各実験においては、それぞれにこれを引き継いだクラスを作って、その中で色々こねくり回すつもり
  * 詳しくはexp_formatを参照。
  */
trait ParamDimCond {
  def readParams(paramFile: String): Map[String, String] = {
    Source.fromFile(paramFile).getLines().foldLeft(Map.empty[String, String]) { (m, s) =>
      if (s == "") m
      else {
        if (s.indexOf("=") < 0) {
          System.err.println("invalid input parameters format!")
          m
        } else {
          val s2 = s.split("=", -1)
          m.updated(s2(0), s2(1))
        }
      }
    }
  }

  implicit def StringToArrayInt(s: String): Array[Int] = s.split("\\s+").map(_.toInt)

  implicit def StringToArrayDouble(s: String): Array[Double] = s.split("\\s+").map(_.toDouble)

  implicit def StringToArrayString(s: String): Array[String] = s.split("\\s+")

  //文字,数値のparamsの時に使う。arrayのparamsがない場合
  private def noArrayParams(s: Array[String]) = {
    def stringEmpty(s: => Array[String]): Boolean = s(1) == ""
    s.length < 2 || stringEmpty(s)
  }

  implicit def StringToStringArrayInt(s: String): (String, Array[Int]) = {
    val s2 = s.split(",", -1)
    val a = if (noArrayParams(s2)) Array(): Array[Int] else s2(1).split("\\s+").map(_.toInt)
    (s2(0), a)
  }

  implicit def StringToStringArrayDouble(s: String): (String, Array[Double]) = {
    val s2 = s.split(",", -1)
    val a = if (noArrayParams(s2)) Array(): Array[Double] else s2(1).split("\\s+").map(_.toDouble)
    (s2(0), a)
  }

  implicit def StringToArrayStringArrayDouble(s: String): (Array[String], Array[Double]) = {
    val s2 = s.split(",", -1)
    val a = if (noArrayParams(s2)) Array(): Array[Double] else s2(1).split("\\s+").map(_.toDouble)
    (s2(0).split("\\s+"), a)
  }

  implicit def StringToStringInt(s: String): (String, Int) = {
    val s2 = s.split(",")
    val a = if (noArrayParams(s2)) -1 else s2(1).toInt
    (s2(0), a)
  }

  implicit def StringToGenerateCond(s: String): GenerateCond = GenerateCond(s)

  implicit def StringToBoolean(s: String): Boolean = s.toBoolean

  implicit def StringToOptCond(s: String): OptCond = OptCond(s)

}


/**
  * optimize condition.
  */
trait OptimizeCondTrait {
  val emHid: (String, Array[Int])
  val emTime: Array[Int]
  val emShallow: Array[Int]
  val emRand: Array[Int]
  val Ainit: GenerateCond
  val Binit: GenerateCond
  val Hinit: GenerateCond
  val Qinit: GenerateCond
  val Rinit: GenerateCond
  val initStateMeanInit: GenerateCond
  val initStateCovarianceInit: GenerateCond
  val delta: Array[Double]
  val AOpt: OptCond
  val BOpt: OptCond
  val HOpt: OptCond
  val QOpt: OptCond
  val ROpt: OptCond
  val initStateMeanOpt: OptCond
  val initStateCovarianceOpt: OptCond
  val parallelThreadNum: Array[Int]
  val crossValidPredictFile: String
  val foldNum: Array[Int]
}

//Optimizationの実験条件
//emHid以外の実験条件のarrayは実はお飾りで、一番上に入ってる値しか使われない。
case class OptimizeCond(emHid: (String, Array[Int]),
                        foldNum: Array[Int],
                        emTime: Array[Int],
                        emShallow: Array[Int],
                        emRand: Array[Int],
                        Ainit: GenerateCond,
                        Binit: GenerateCond,
                        Hinit: GenerateCond,
                        Qinit: GenerateCond,
                        Rinit: GenerateCond,
                        initStateMeanInit: GenerateCond,
                        initStateCovarianceInit: GenerateCond,
                        delta: Array[Double],
                        AOpt: OptCond,
                        BOpt: OptCond,
                        HOpt: OptCond,
                        QOpt: OptCond,
                        ROpt: OptCond,
                        initStateMeanOpt: OptCond,
                        initStateCovarianceOpt: OptCond,
                        parallelThreadNum: Array[Int],
                        crossValidPredictFile: String) extends OptimizeCondTrait

/**
  *
  */
object OptimizeCond extends ParamDimCond {

  def apply(paramFile: String): OptimizeCond = {
    val params = readParams(paramFile)
    val emHid = params.getOrElse("emHid", "")
    val foldNum = params.getOrElse("foldNum", "")
    val emTime = params.getOrElse("emTime", "")
    val emShallow = params.getOrElse("emShallow", "")
    val emRand = params.getOrElse("emRand", "")
    val Ainit = params.getOrElse("Ainit", "")
    val Binit = params.getOrElse("Binit", "")
    val Hinit = params.getOrElse("Hinit", "")
    val Qinit = params.getOrElse("Qinit", "")
    val Rinit = params.getOrElse("Rinit", "")
    val initStateMeanInit = params.getOrElse("initStateMeanInit", "")
    val initStateCovarianceInit = params.getOrElse("initStateCovarianceInit", "")
    val delta = params.getOrElse("delta", "")
    val AOpt = params.getOrElse("AOpt", "")
    val BOpt = params.getOrElse("BOpt", "")
    val HOpt = params.getOrElse("HOpt", "")
    val QOpt = params.getOrElse("QOpt", "")
    val ROpt = params.getOrElse("ROpt", "")
    val initStateMeanOpt = params.getOrElse("initStateMeanOpt", "")
    val initStateCovarianceOpt = params.getOrElse("initStateCovarianceOpt", "")
    val parallelThreadNum = params.getOrElse("parallelThreadNum", "")
    val crossValidPredictFile = params.getOrElse("crossValidPredictFile", "")

    OptimizeCond(emHid, foldNum, emTime, emShallow, emRand, Ainit, Binit, Hinit, Qinit, Rinit,
      initStateMeanInit, initStateCovarianceInit, delta, AOpt, BOpt, HOpt, QOpt, ROpt,
      initStateMeanOpt, initStateCovarianceOpt, parallelThreadNum, crossValidPredictFile)
  }
}

//
case class EMCond(emTime: Int = -1,
                  delta: Double = 0.01,
                  AOpt: OptCond = OptCond(""),
                  BOpt: OptCond = OptCond(""),
                  HOpt: OptCond = OptCond(""),
                  QOpt: OptCond = OptCond(""),
                  ROpt: OptCond = OptCond(""),
                  initStateMeanOpt: OptCond = OptCond(""),
                  initStateCovarianceOpt: OptCond = OptCond(""),
                  showLikelihood: Boolean = false,
                  initkf: KalmanFilter, initStateMeans: List[DenseVector[Double]],
                  initStateCovariances: List[Matrices])

//
case class ExperimentCond(hDim: Array[Int],
                           cDim: Array[Int],
                           oDim: (String, Array[Int]),
                           seqNum: Array[Int],
                           seqLen: Array[Int],
                           observeLocation: (String, Array[Int]),
                           Agen: GenerateCond,
                           Bgen: GenerateCond,
                           Hgen: GenerateCond,
                           Qgen: GenerateCond,
                           Rgen: GenerateCond,
                           initStateMeanGen: GenerateCond,
                           initStateCovarianceGen: GenerateCond,
                           controlGen: GenerateCond,
                           emHid: (String, Array[Int]),
                           foldNum: Array[Int],
                           emTime: Array[Int],
                           emShallow: Array[Int],
                           emRand: Array[Int],
                           Ainit: GenerateCond,
                           Binit: GenerateCond,
                           Hinit: GenerateCond,
                           Qinit: GenerateCond,
                           Rinit: GenerateCond,
                           initStateMeanInit: GenerateCond,
                           initStateCovarianceInit: GenerateCond,
                           delta: Array[Double],
                           fixDelta: Array[Double],
                           AOpt: OptCond,
                           BOpt: OptCond,
                           HOpt: OptCond,
                           QOpt: OptCond,
                           ROpt: OptCond,
                           initStateMeanOpt: OptCond,
                           initStateCovarianceOpt: OptCond,
                           testSeqNum: Array[Int],
                           parallelThreadNum: Array[Int],
                           crossValidPredictFile: String,
                           shFile: String,
                           jobId: String,
                           version: String) extends ParamDimCond with OptimizeCondTrait


object ExperimentCond extends ParamDimCond {
  def apply(paramFile: String): ExperimentCond = {
    val params = readParams(paramFile)
    val hDim = params.getOrElse("hDim", "")
    val cDim = params.getOrElse("cDim", "")
    val oDim = params.getOrElse("oDim", "")
    val seqNum = params.getOrElse("seqNum", "")
    val seqLen = params.getOrElse("seqLen", "")
    val observeLocation = params.getOrElse("observeLocation", "")
    val Agen = params.getOrElse("Agen", "")
    val Bgen = params.getOrElse("Bgen", "")
    val Hgen = params.getOrElse("Hgen", "")
    val Qgen = params.getOrElse("Qgen", "")
    val Rgen = params.getOrElse("Rgen", "")
    val initStateMeanGen = params.getOrElse("initStateMeanGen", "")
    val initStateCovarianceGen = params.getOrElse("initStateCovarianceGen", "")
    val controlGen = params.getOrElse("controlGen", "")
    val emHid = params.getOrElse("emHid", "")
    val foldNum = params.getOrElse("foldNum", "")
    val emTime = params.getOrElse("emTime", "")
    val emShallow = params.getOrElse("emShallow", "")
    val emRand = params.getOrElse("emRand", "")
    val Ainit = params.getOrElse("Ainit", "")
    val Binit = params.getOrElse("Binit", "")
    val Hinit = params.getOrElse("Hinit", "")
    val Qinit = params.getOrElse("Qinit", "")
    val Rinit = params.getOrElse("Rinit", "")
    val initStateMeanInit = params.getOrElse("initStateMeanInit", "")
    val initStateCovarianceInit = params.getOrElse("initStateCovarianceInit", "")
    val delta = params.getOrElse("delta", "")
    val fixDelta = params.getOrElse("fixDelta", "")
    val AOpt = params.getOrElse("AOpt", "")
    val BOpt = params.getOrElse("BOpt", "")
    val HOpt = params.getOrElse("HOpt", "")
    val QOpt = params.getOrElse("QOpt", "")
    val ROpt = params.getOrElse("ROpt", "")
    val initStateMeanOpt = params.getOrElse("initStateMeanOpt", "")
    val initStateCovarianceOpt = params.getOrElse("initStateCovarianceOpt", "")
    val testSeqNum = params.getOrElse("testSeqNum", "")
    val parallelThreadNum = params.getOrElse("parallelThreadNum", "")
    val crossValidPredictFile = params.getOrElse("crossValidPredictFile", "")
    val shFile = params.getOrElse("shFile", "")
    val jobId = params.getOrElse("jobId", "")
    val version = params.getOrElse("version", "")
    ExperimentCond(hDim, cDim, oDim, seqNum, seqLen, observeLocation, Agen, Bgen, Hgen, Qgen, Rgen,
      initStateMeanGen, initStateCovarianceGen, controlGen, emHid, foldNum, emTime, emShallow, emRand, Ainit,
      Binit, Hinit, Qinit, Rinit, initStateMeanInit, initStateCovarianceInit, delta, fixDelta,
      AOpt, BOpt, HOpt, QOpt, ROpt, initStateMeanOpt, initStateCovarianceOpt, testSeqNum,
      parallelThreadNum, crossValidPredictFile, shFile, jobId, version)
  }
}