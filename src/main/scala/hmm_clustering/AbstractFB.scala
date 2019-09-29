package hmm_clustering

/**
  * Created by tara on 2019/09/29.
  */
trait AbstractFB {

  val hmm: HMM

  //data set of control and observed variables. It handles only one sequence.
  val seq: Array[Double]

  val T = seq.length

  val n = hmm.transit.length

  val run: Array[Array[Double]]
}

/**
  * This class conducts forward algorithm.
  * @param hmm
  * @param seq
  */
case class Forward(hmm: HMM, seq: Array[Double]) extends AbstractFB {

  lazy val run: Array[Array[Double]] = {
    val forward = Array.fill(T + 1)(Array.fill[Double](n)(0.0))
    //initialization
    forward(0)(0) = 1.0

    (0 until T).foreach(i => {
      val si = seq(i)
      val prev = forward(i)

      val f_e_k = prev.zip(hmm.calcLikelihood(si)).map(x => x._1 * x._2)
      val new_f = (0 until n).map(l => {
        (0 until n).map(k => {
          f_e_k(k) * hmm.transit(k)(l)
        }).sum
      }).toArray
      forward(i + 1) = new_f
    })
    forward
  }

  def likelihood: Double = run(T).sum

}

/**
  * This class conducts forward algorithm.
  * @param hmm
  * @param seq
  */
case class Backward(hmm: HMM, seq: Array[Double]) extends AbstractFB {

  lazy val run: Array[Array[Double]] = {
    val backward = Array.fill(T + 1)(Array.fill[Double](n)(0.0))
    //initialization
    backward(T) = Array.fill(n)(1.0)

    (T - 1 to 0 by -1).foreach(i => {
      val si = seq(i)
      val prev = backward(i + 1)

      val akl_bl = hmm.transit.map(ak => {
        ak.zip(prev).map(x => x._1 * x._2).sum
      })
      val new_b = hmm.calcLikelihood(si).zip(akl_bl).map(x => x._1 * x._2)
      backward(i) = new_b
    })
    backward
  }

  def likelihood: Double = run(0)(0)
}
