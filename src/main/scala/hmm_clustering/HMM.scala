package hmm_clustering

/**
  * Created by tara on 2019/09/29.
  */
case class HMM(mu: Array[Double], sigma: Array[Double], transit: Array[Array[Double]]) {

  assert(sigma.length == mu.length)
  assert(transit.length == mu.length)
  val m = mu.length

  def calcLikelihood(seq: Array[Double]): Double = {
    val forward = Forward(this, seq)
    forward.likelihood
  }

  def calcLikelihood(obs: Double): Array[Double] = {
    mu.zip(sigma).map(x => {
      val mmm = x._1
      val s = x._2
      Math.exp(- (obs - mmm) * (obs - mmm) / (2 * s * s)) / Math.sqrt(2 * Math.PI * s * s)
    })
  }

  def generate(n: Int): (Array[(Double, Int)]) = {
    def decideTransitLoop(transitArray: Array[Double], randomVal: Double, index: Int): Int = {
      if (transitArray(index) >= randomVal || index == m - 1) {
        index
      } else {
        decideTransitLoop(transitArray, randomVal, index + 1)
      }
    }
    val gaussianDistribution = mu.zip(sigma).map(x => breeze.stats.distributions.Gaussian(x._1, x._2))
    val culm = transit.map(x => x.scanLeft(0.0)(_+_).slice(1, m + 1))
    val initState = 0
    val initOutput = gaussianDistribution(0).sample()
    val transitRand = breeze.stats.distributions.Uniform(0.0, 1.0)

    (0 until n).scanLeft(initOutput, initState) {(s, x) =>
      val prevState = s._2
      val tr = transitRand.sample()
      val newState = decideTransitLoop(culm(prevState), tr, 0)
      val outputRand = gaussianDistribution(newState).sample()
      (outputRand, newState)
    }.toArray
  }
}

object HMM {
  def createInitTransitRL(n: Int): Array[Array[Double]] = {
    val transit = Array.fill[Array[Double]](n)(Array.fill[Double](n)(0.0))
    (0 until n).foreach(i => {
      (i until n).foreach(j => {
        if (i == n-1) {
          transit(i)(j) = 1.0
        } else if (i == n-2) {
          transit(i)(i) = 0.8
          transit(i)(i+1) = 0.2
        } else if (j == i + 1) {
          transit(i)(j) = 0.2
        } else if (i == j) {
          transit(i)(j) = 0.6
        } else if (j >= i) {
          transit(i)(j) = 0.2  / (n - i - 2)
        }
      })
    })
    transit
  }
}