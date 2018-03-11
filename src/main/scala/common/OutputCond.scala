package common

import java.lang.reflect.Field

/**
  * Created by tara on 2017/05/10.
  * outputs of some in tsv.
  */
trait OutputCond {
  val version: String
  val shFile: String

  def classNameToString(): String = toStringCommon(f => f.getName())

  override def toString(): String = toStringCommon(f => {
    f.setAccessible(true)
    f.get(this).toString()
  })

  def toStringCommon(fun: (Field => String)): String = {
    val sbr = new StringBuilder()
    val fields: Array[Field] = this.getClass().getDeclaredFields()
    fields.zipWithIndex.foreach(f => {
      sbr.append(fun(f._1)); if (f._2 < fields.length - 1) sbr.append("\t") else sbr.append("\n")
    })
    sbr.toString()
  }
}

/**
  * the output condition of experiment
  */
case class ExperimentOutputCond(hDim: Int, //dimension of hidden variables
                                cDim: Int, //dimension of control variables
                                oDim: Int, //dimension of observed variables
                                seqNum: Int, //number of sequence
                                seqLen: Int, //
                                emHidCond: String, //hidden condition
                                emHidDim: Int,
                                foldNum: Int,
                                emTime: Int,
                                emShallow: Int,
                                emRand: Int,
                                delta: Double,
                                fixDelta: Double,
                                Aerror: Double,
                                Berror: Double,
                                Qerror: Double,
                                Herror: Double,
                                Rerror: Double,
                                initStateMeanError: Double,
                                initStateCovarianceError: Double,
                                calcTime: Long,
                                logLikelihood: Double,
                                predCorr: Double,
                                realLikelihood: Double,
                                realCorrPred: Double,
                                testSeqNum: Int,
                                parallelThreadNum: Int,
                                shFile: String,
                                jobId: String,
                                version: String
                               ) extends OutputCond

object ExperimentOutputCond {
  def apply(): ExperimentOutputCond = {
    new ExperimentOutputCond(hDim = -1,
      cDim = -1,
      oDim = -1,
      seqNum = -1,
      seqLen = -1,
      emHidCond = "",
      emHidDim = -1,
      foldNum = -1,
      emTime = -1,
      emShallow = -1,
      emRand = -1,
      delta = -1,
      fixDelta = -1,
      Aerror = -1,
      Berror = -1,
      Qerror = -1,
      Herror = -1,
      Rerror = -1,
      initStateMeanError = -1,
      initStateCovarianceError = -1,
      calcTime = -1,
      logLikelihood = -1,
      predCorr = -1,
      realLikelihood = -1,
      realCorrPred = -1,
      testSeqNum = -1,
      parallelThreadNum = -1,
      shFile = "",
      jobId = "",
      version = "")
  }
}