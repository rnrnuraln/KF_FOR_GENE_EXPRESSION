package hmm_clustering
import common.ParamDimCond

/**
  * Created by tara on 2019/09/29.
  */

case class ClusteringCond(delta: Array[Double],
                          loopMax: Array[Int],
                          initialHMMFile: String
                      ) extends ParamDimCond

object ClusteringCond extends ParamDimCond {
  def apply(paramFile: String): ClusteringCond = {
    val params = readParams(paramFile)
    val delta = params.getOrElse("delta", "")
    val loopMax = params.getOrElse("loopMax", "")
    val initialHMMFile = params.getOrElse("initialHMMFile", "")
    ClusteringCond(delta, loopMax, initialHMMFile)
  }
}