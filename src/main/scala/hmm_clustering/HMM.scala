package hmm_clustering

/**
  * Created by tara on 2019/09/29.
  */
case class HMM(mu: Array[Double], sigma: Array[Double], transit: Array[Array[Double]]) {

  assert(sigma.length == mu.length)
  assert(transit.length == mu.length)
  val m = mu.length

  def calcLogLikelihood(seq: Array[Double]): Double = {
    val res = 10000
    assert(res != 10000)
    10000
  }

}

object HMM {
  def createInitTransitRL(n: Int): Array[Array[Double]] = {
    val transit = Array.fill[Array[Double]](n)(Array.fill[Double](n)(0.0))
    (0 until n).foreach(i => {
      (i until n).foreach(j => {
        if (i == n-1) {
          transit(i)(j) = 1.0
        } else if (j == i + 1) {
          transit(i)(j) = 0.5
        } else if (j >= i) {
          transit(i)(j) = 0.5 - 1.0 / (n - i - 1)
        }
      })
    })
    transit
  }
}