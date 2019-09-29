package hmm_clustering

/**
  * Created by tara on 2019/09/29.
  */
class Baum_Welch(initHMM: HMM, geneSeqs: Seq[Array[Double]]) {

  val m = initHMM.m
  val T = geneSeqs(0).length

  def run(): HMM = {
    def learnSub(prevHMM: HMM, i: Int, forwards: Seq[Forward]): HMM = {
      val prevLogLikelihood = forwards.map(_.logLikelihood).sum
      val backwards = geneSeqs.map(x => Backward(prevHMM, x))
      val newHMM = updateParams(forwards.map(_.run), backwards.map(_.run), forwards.map(_.logLikelihood), prevHMM)
      val newForwards = geneSeqs.map(x => Forward(newHMM, x))
      val logLikelihood = newForwards.map(_.logLikelihood).sum
      if (logLikelihood.isNaN || logLikelihood.isInfinity || logLikelihood + 10.0 < prevLogLikelihood) {
        //あり得ない結果
        prevHMM
      } else if (i > 20) {
        //終了条件
        newHMM
      } else {
        learnSub(newHMM, i + 1, newForwards)
      }
    }
    learnSub(initHMM, 0, geneSeqs.map(x => Forward(initHMM, x)))
  }

  def updateParams(forwardSeq: Seq[Array[Array[Double]]], backwardSeq: Seq[Array[Array[Double]]], logLikelihoodSeq: Seq[Double], oldHMM: HMM): HMM = {
    val seqsZip = forwardSeq.zip(backwardSeq).zip(logLikelihoodSeq.zip(geneSeqs))
    val newTransit = (0 until oldHMM.m).map(k => {
      (k until oldHMM.m).map(l => {
        seqsZip.map(x => {
          val forwards = x._1._1
          val backwards = x._1._2
          val ll = x._2._1
          val genes = x._2._2
          val geneAkl = (0 until T).foldLeft(0.0) {(s, i) =>
            val f_k_i = forwards(i)(k)
            val a_k_l = oldHMM.transit(k)(l)
            val e_k_i = oldHMM.calcLikelihood(genes(i))(k)
            val b_l_i_1 = backwards(i+1)(l)
            s + f_k_i * a_k_l * e_k_i * b_l_i_1
          }
          geneAkl / ll
        }).sum
      }).toArray
    }).toArray
    val newMu = (0 until oldHMM.m).map(k => {
      seqsZip.map(x => {
        val forwards = x._1._1
        val backwards = x._1._2
        val genes = x._2._2
        (0 until T).foldLeft(0.0) {(s, i) =>
          val f_k_i = forwards(i)(k)
          val b_l_i = backwards(i)(k)
          s + f_k_i * b_l_i * genes(i)
        }
      }).sum / (T * geneSeqs.length)
    }).toArray

    HMM(newMu, oldHMM.sigma, newTransit)
  }

}

object Baum_Welch {

}