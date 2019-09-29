package hmm_clustering

/**
  * Created by tara on 2019/09/29.
  */
class Baum_Welch(initHMM: HMM, geneSeqs: Seq[Array[Double]]) {

  val m = initHMM.m
  val T = geneSeqs(0).length

  def run(): HMM = {
    def learnSub(prevHMM: HMM, i: Int, forwards: Seq[Forward]): HMM = {
      val prevLogLikelihood = forwards.map(_.likelihood).sum
      // println(i + "th: prev likelihood: " + prevLogLikelihood)
      val backwards = geneSeqs.map(x => Backward(prevHMM, x))
      val newHMM = updateParams(forwards.map(_.run), backwards.map(_.run), forwards.map(_.likelihood), prevHMM)
      val newForwards = geneSeqs.map(x => Forward(newHMM, x))
      val likelihood = newForwards.map(_.likelihood).sum
      if (likelihood.isNaN || likelihood.isInfinity) {
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

  def updateParams(forwardSeq: Seq[Array[Array[Double]]], backwardSeq: Seq[Array[Array[Double]]], likelihoodSeq: Seq[Double], oldHMM: HMM): HMM = {
    val fb_prob: Seq[Array[Double]] = forwardSeq.zip(backwardSeq).map(x => {
      x._1.zip(x._2).map(y => {
        y._1.zip(y._2).map(z => z._1 * z._2).sum
      })
    })

    val fb_prob_2: Seq[Array[Double]] = forwardSeq.zip(backwardSeq).zip(geneSeqs).map(x => {
      val forwards = x._1._1
      val backwards = x._1._2
      val genes = x._2
      (0 until T).map(i => {
        (0 until m).map(k => {
          (0 until m).map(l => {
            val f_k_i = forwards(i)(k)
            val a_k_l = oldHMM.transit(k)(l)
            val e_k_i = oldHMM.calcLikelihood(genes(i))(k)
            val b_l_i_1 = backwards(i + 1)(l)
            f_k_i * a_k_l * e_k_i * b_l_i_1
          }).sum
        }).sum
      }).toArray
    })

    val seqsZip = forwardSeq.zip(backwardSeq).zip(fb_prob_2.zip(geneSeqs))
    val newTransit = (0 until oldHMM.m).map(k => {
      val rawTransit = (0 until oldHMM.m).map(l => {
        if (l < k) 0.0 else {
          seqsZip.map(x => {
            val forwards = x._1._1
            val backwards = x._1._2
            val fbprob = x._2._1
            val genes = x._2._2
            (0 until T).foldLeft(0.0) { (s, i) =>
              val f_k_i = forwards(i)(k)
              val a_k_l = oldHMM.transit(k)(l)
              val e_k_i = oldHMM.calcLikelihood(genes(i))(k)
              val b_l_i_1 = backwards(i + 1)(l)
              val klprob = f_k_i * a_k_l * e_k_i * b_l_i_1 / fbprob(i)
              // println(k + ", " + l + "th prob: " + klprob)
              s + klprob
            }
          }).sum
        }
      }).toArray
      val rtSum = rawTransit.sum
      rawTransit.map(_ / rtSum)
    }).toArray
    val newMu = (0 until oldHMM.m).map(k => {
      val a = seqsZip.map(x => {
        val forwards = x._1._1
        val backwards = x._1._2
        val fbprob = x._2._1
        val genes = x._2._2
        val ak = (0 until T).foldLeft(0.0) {(s, i) =>
          val f_k_i = forwards(i)(k)
          val b_k_i = backwards(i)(k)
          val k_prob = (f_k_i * b_k_i) / fbprob(i)
          // println(k + "th prob: " + k_prob)
          s + k_prob * genes(i)
        }
        ak
      }).sum
      val b = seqsZip.map(x => {
        val forwards = x._1._1
        val backwards = x._1._2
        val fbprob = x._2._1
        val genes = x._2._2
        val ak = (0 until T).foldLeft(0.0) {(s, i) =>
          val f_k_i = forwards(i)(k)
          val b_k_i = backwards(i)(k)
          val k_prob = (f_k_i * b_k_i) / fbprob(i)
          // println(k + "th prob: " + k_prob)
          s + k_prob
        }
        ak
      }).sum
      a / b
    }).toArray

    HMM(newMu, oldHMM.sigma, newTransit)
  }

}

object Baum_Welch {

}