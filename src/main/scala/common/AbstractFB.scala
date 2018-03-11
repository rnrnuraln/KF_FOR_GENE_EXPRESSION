package common

import breeze.linalg._

/**
  * Created by tara on 2016/12/23.
  * Abstract class of forward and backward algorithm, using Kalman Filter model.
  * More details will be written in docs/kalman_filter.pdf
  */
trait AbstractFB[V <: AbstractGaussianTransit] {

  //Kalman Filter used for calculation of each algorithm
  val kf: KalmanFilter
  //data set of control and observed variables. It handles only one sequence.
  val seqs: Array[common.ConObs]

  // initial state is considered to be stochastically generated from a Gaussian distribution.
  // parameters below define its mean and covariance
  val initStateMean: DenseVector[Double]
  val initStateCovariance: Matrices

  //these parameters are used common for calculating forward & backward algorithm.
  protected val invR = kf.R.inv
  protected val HR: Matrices = kf.H.t * invR
  protected val HRH = HR * kf.H
  protected val lnRln2pi = Math.log(kf.R.det) + kf.m * Math.log(2 * Math.PI)

  val T = seqs.length
  val n = kf.n
  val m = kf.m

  val run: Array[V]
}

/**
  * This class conducts forward algorithm.
  * In fact, this conducts so-called Kalman Filter.
  * @param kf
  * @param seqs
  * @param initStateMean
  * @param initStateCovariance
  */
case class Forward(kf: KalmanFilter, seqs: Array[common.ConObs],
                   initStateMean: DenseVector[Double],
                   initStateCovariance: Matrices) extends AbstractFB[ForwardGaussian] {

  //this has Forward Gaussian, B and C in docs/kalman_filter.pdf, of each time point
  val run: Array[ForwardGaussian] = {
    val forward = Array.fill(T + 1)(ForwardGaussian())
    //initialization
    forward(0) = ForwardGaussian(initStateCovariance.inv * initStateMean, initStateCovariance)

    for (i <- 0 until T) {
      val si = seqs(i)
      val prev = forward(i)
      val C0_D1 = si.observe match {
        case Some(o) => (prev.C.inv + HRH).inv
        case None => prev.C
      }
      val Bu = si.control match {
        case Some(c) => kf.B.getOrElse(DenseMatrices(DenseMatrix.zeros[Double](n, kf.l))) * c;
        case None => DenseVector.zeros[Double](n)
      }
      val C: Matrices = kf.Q + kf.A * C0_D1 * kf.A.t
      val invC = C.inv
      val B = si.observe match {
        case Some(o) => invC * (Bu + (kf.A * C0_D1 * (prev.B + HR * o)))
        case None => invC * (Bu + kf.A * prev.C * prev.B)
      }
      forward(i + 1) = ForwardGaussian(B, C)
    }
    forward
  }

  def logLikelihood: Double = -0.5 * run.zip(seqs).foldLeft(0.0) { (s, x) =>
    calcLogLikelihood(x._1, x._2.observe, s)
  }

  private def calcLogLikelihood(g: ForwardGaussian, observe: Option[DenseVector[Double]], prev: Double): Double = {
    observe match {
      case Some(o) =>
        val C0_D1 = (g.C.inv + HRH).inv
        val zHCB = o - kf.H * g.C * g.B
        val HRzHCB = HR * zHCB
        prev + Math.log(g.C.det) - Math.log(C0_D1.det) + lnRln2pi + zHCB.t * (invR * zHCB) - HRzHCB.t * (C0_D1 * HRzHCB)
      case None => prev
    }
  }

  //Backwardと合わせたときの検証用に各時刻での対数尤度を求める
  //this is used to
  def logLikelihoodArray: Array[Double] = {
    val array = Array.fill(T + 1)(0.0)
    array(0) = 0.0
    for (i <- 0 until T) {
      array(i + 1) = calcLogLikelihood(run(i), seqs(i).observe, array(i))
    }
    array
  }

  //カルマンフィルタを用いて、逐次観測変数の値を予測する。
  //その時のcovarianceも求める
  def predictsWithCovariance(): Array[(Matrices, DenseVector[Double])] = run.map(x => {
    val invCD = (x.C.inv + HRH).inv
    val K = invR - HR.t * invCD * HR
    (K, K.inv * kf.H * x.C * x.B)
  })

  //covarianceを用いない観測変数の予測。時間に気をつけるのじゃ。t = 0 ~ Tまで
  def predicts: Array[DenseVector[Double]] = run.map(x => kf.H * x.C * x.B)

  //予測値と実測値を比較し、時系列上の相関係数を求める。
  //sum(z(i+1) - z(i))(z2(i+1)-z(i))/sum(z(i+1) - z(i))^2を求める
  //prevsの1番目には相関係数の分子、二番目には分母、z_None: Option[DenseVector[Double]]となっているとこには前の値が入る
  //相関係数は、ベクトルの各値ごとに求めるものとする。
  def predCorrelation: DenseVector[Double] = {
    val prevs = predicts.zip(seqs).foldLeft(DenseVector.zeros[Double](m), DenseVector.zeros[Double](m), None: Option[DenseVector[Double]]) { (s, x) =>
      x._2.observe match {
        case Some(o) =>
          s._3 match {
            case Some(p) =>
              val RealPrevDif: DenseVector[Double] = o - p
              val predictDif: DenseVector[Double] = x._1 - p
              val child = s._1 + (RealPrevDif *:* predictDif)
              val mom = s._2 + RealPrevDif *:* RealPrevDif
              (child, mom, Some(o))
            case None => (s._1, s._2, Some(o))
          }
        case None => s
      }
    }
    prevs._1 /:/ prevs._2
  }

  //prediction correlation
  def predCorrEval: Double = predCorrelation.foldLeft(0.0) { (s, x) => s + (x - 1) * (x - 1) }
}

/**
  * This class conducts backward algorithm.
  * Together with forward algorithm, so-called Kalman Smoothing.
  * @param kf
  * @param seqs
  * @param initStateMean
  * @param initStateCovariance
  */
case class Backward(kf: KalmanFilter, seqs: Array[common.ConObs],
                    initStateMean: DenseVector[Double],
                    initStateCovariance: Matrices) extends AbstractFB[BackwardGaussian] {

  //observeが存在する最後のindex
  val end = {
    def calcEnd(i: Int): Int = {
      seqs(i).observe match {
        case Some(o) => i
        case None => calcEnd(i - 1)
      }
    }
    calcEnd(T - 1)
  }

  val run: Array[BackwardGaussian] = {
    val backward = Array.fill(end + 1)(BackwardGaussian())

    //initialization
    backward(end) = BackwardGaussian(HRH, HR * seqs(end).observe.getOrElse(DenseVector.zeros[Double](m)))

    for (i <- end - 1 to 0 by -1) {
      val si = seqs(i)
      val prev = backward(i + 1)
      val invQD = (kf.Q.inv + prev.D).inv
      //expresses (D.inv + Q).inv using invQD (considering possibility of D being singular)
      val invDQ = prev.D - prev.D * invQD * prev.D
      val D = si.observe match {
        case Some(o) => HRH + kf.A.t * invDQ * kf.A
        case None => kf.A.t * invDQ * kf.A
      }
      val bc = si.control match {
        case Some(c) => kf.B.getOrElse(DenseMatrices(DenseMatrix.zeros[Double](n, kf.l))) * c
        case None => DenseVector.zeros[Double](n)
      }
      val befE = kf.A.t * (kf.Q.inv * (invQD * (prev.E - prev.D * bc)))
      val E = si.observe match {
        case Some(o) => HR * o + befE
        case None => befE
      }
      backward(i) = BackwardGaussian(D, E)
    }
    backward
  }

  def logLikelihoodArray: Array[Double] = {
    val invQ = kf.Q.inv
    val lnQ = Math.log(kf.Q.det)
    val array = Array.fill[Double](T + 1)(0.0)

    val o = seqs(end).observe.getOrElse(DenseVector.zeros[Double](m))
    array(end) = lnRln2pi + o.t * (invR * o)
    for (i <- end - 1 to 0 by -1) {
      val control = seqs(i).control
      val observe = seqs(i).observe
      val back = run(i + 1)
      val invQD = (invQ + back.D).inv
      val bc = control match {
        case Some(c) => kf.B.getOrElse(DenseMatrices(DenseMatrix.zeros[Double](n, kf.l))) * c
        case None => DenseVector.zeros[Double](n)
      }
      val edb = back.E - back.D * bc
      val s2 = array(i + 1) + lnQ - Math.log(invQD.det) - edb.t * (invQD * edb) + bc.t * (back.D * bc) - 2 * (bc.t * back.E)
      array(i) = observe match {
        case Some(o) => s2 + lnRln2pi + o.t * (invR * o)
        case None => s2
      }
    }
    array
  }
}