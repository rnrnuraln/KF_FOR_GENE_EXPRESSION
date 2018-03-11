package common

import breeze.linalg._

/**
  * using EM algorithm.
  *
  * Created by tara on 2017/04/20.
  */
case class EMAlgorithm(seqs: List[Array[common.ConObs]],
                       cond: EMCond) {

  //dimension of control variables
  val l = seqs(0)(0).control match {
    case Some(c) => c.length;
    case None => 0
  }

  //dimension of observed variables
  val m = {
    def defineM(i: Int): Int = seqs(0)(i).observe match {
      case Some(o) => o.length;
      case None => defineM(i + 1)
    }
    defineM(0)
  }

  //dimension of hidden variables
  val n = cond.initkf.n

  implicit def matrixToMatrices(matrix: DenseMatrix[Double]): Matrices = DenseMatrices(matrix)

  // the number in seq which has observed variable
  val dataSize = seqs.foldLeft(0) { (s, x) =>
    s + x.count(y => y.observe.isDefined)
  }

  lazy val allDataMean = seqs.foldLeft(DenseVector.zeros[Double](m)) {(s, x) =>
    s + x.foldLeft(DenseVector.zeros[Double](m)) {(s2, y) =>
      s2 + y.observe.getOrElse(DenseVector.zeros[Double](m))
    }
  } * (1.0 / dataSize)

  lazy val allDataVariance = seqs.foldLeft(DenseVector.zeros[Double](m)) {(s, x) =>
    s + x.foldLeft(DenseVector.zeros[Double](m)) {(s2, y) =>
      y.observe match {
        case Some(o) => s2 + (o - allDataMean) *:* (o - allDataMean)
        case None => s2
      }
    }
  } * (1.0 / dataSize)

  /**
    * Since we cannot update all parameters simultaneously, we have to define which parameters be updated when updating.
    * More precisely, each parameter is divided into three groups; ABQ, HR, initStateMean & Covariance,
    * and parameters in each group cannot be updated simultaneously.
    * This class defines turn for each group.
    * @param abq
    * @param hr
    * @param init
    */
  class TurnForEM(abq: Char = 'A', hr: Char = 'H', init: Char = 'm') {
    val Afixed = cond.AOpt._1(0) == "fixed"
    val Bfixed = cond.BOpt._1(0) == "fixed" || l == 0
    val Qfixed = cond.QOpt._1(0) == "fixed"
    val Hfixed = cond.HOpt._1(0) == "fixed" || (cond.HOpt._1.indexOf("upperIdentity") > 0 && (m == n))
    val Rfixed = cond.ROpt._1(0) == "fixed"
    val initStateMeanfixed = cond.initStateMeanOpt._1(0) == "fixed"
    val initStateCovariancefixed = cond.initStateCovarianceOpt._1(0) == "fixed"
    val abqTurn = if (Afixed && Bfixed && Qfixed) 'F' else abq
    val hrTurn = if (Hfixed && Rfixed) 'F' else hr
    val initTurn = if (initStateCovariancefixed && initStateMeanfixed) 'F' else init

    def rotate(): TurnForEM = {
      def abqRotate(c: Char): Char = {
        c match {
          case 'A' => if (!Bfixed) 'B' else abqRotate('B')
          case 'B' => if (!Qfixed) 'Q' else abqRotate('Q')
          case 'Q' => if (!Afixed) 'A' else abqRotate('A')
          case 'F' => 'F'
          case _ => 'F'
        }
      }
      val abq = abqRotate(abqTurn)
      def hrRotate(c: Char): Char = {
        c match {
          case 'H' => if (!Rfixed) 'R' else hrRotate('R')
          case 'R' => if (!Hfixed) 'H' else hrRotate('H')
          case 'F' => 'F'
          case _ => 'F'
        }
      }
      val hr = hrRotate(hrTurn)
      def initStateRotate(c: Char): Char = {
        c match {
          case 'm' => if (!initStateCovariancefixed) 'S' else initStateRotate('S')
          case 'S' => if (!initStateMeanfixed) 'm' else initStateRotate('m')
          case 'F' => 'F'
          case _ => 'F'
        }
      }
      val init = initStateRotate(initTurn)
      new TurnForEM(abq, hr, init)
    }
  }

  /**
    * learn parameters from data using EM algorithm
    * @return
    */
  def learnBasis(): EMOutputs = {
    def learnSub(prevOutputs: EMOutputs, i: Int, turn: TurnForEM): EMOutputs = {
      val seqsWithInits = seqs.zip(prevOutputs.initStateMeans.zip(prevOutputs.initStateCovariance))
      // Kalman smoothing using forward and backward algorithm
      val forwards = seqsWithInits.map(x => Forward(prevOutputs.kf, x._1, x._2._1, x._2._2))
      val backwards = seqsWithInits.map(x => Backward(prevOutputs.kf, x._1, x._2._1, x._2._2))
      val emDtos = obtainEMDto(forwards.map(_.run).zip(backwards.map(_.run)), prevOutputs.kf)

      val outputs = updateParams(emDtos, prevOutputs, turn)
      if (i > cond.emTime || prevOutputs.calcDif(outputs) < cond.delta) {
        val newOutputs = outputs.copy(logLikelihoods = forwards.map(_.logLikelihood).sum :: prevOutputs.logLikelihoods)
        newOutputs
      } else {
        val newOutputs = if (cond.showLikelihood) outputs.copy(logLikelihoods = forwards.map(_.logLikelihood).sum :: prevOutputs.logLikelihoods) else outputs
        learnSub(newOutputs, i + 1, turn.rotate())
      }
    }
    learnSub(EMOutputs(cond.initkf, List(Double.MinValue), cond.initStateMeans, cond.initStateCovariances), 0, new TurnForEM().rotate())
  }

  /**
    * Using results of forward and backward algorithms,
    * @param fbSeq
    * @param kf
    * @return
    */
  def obtainEMDto(fbSeq: Seq[(Array[ForwardGaussian], Array[BackwardGaussian])], kf: KalmanFilter): Seq[Array[EMDto]] = {
    val invR = kf.R.inv
    val invC = kf.Q.inv
    val AC = kf.A.t * invC
    val ACA = AC * kf.A
    val HR = kf.H.t * invR
    val D = HR * kf.H
    val D_ACA = D + ACA
    fbSeq.zip(seqs).map(x => {
      val fRun = x._1._1
      val bRun = x._1._2
      val seq = x._2
      val end = bRun.length
      val emDtoArray = Array.fill(end + 1)(EMDto())
      for (i <- 0 until end) {
        val f = fRun(i)
        val b = bRun(i)
        val invfC = f.C.inv
        val alpha = (invfC + b.D).inv
        val beta = alpha * (f.B + b.E)
        val gamma = seq(i).observe match {
          case Some(o) => (D_ACA + invfC).inv
          case None => (ACA + invfC).inv
        }
        val delta = gamma * AC
        val e1 = seq(i).observe match {
          case Some(o) => HR * o
          case None => DenseVector.zeros[Double](kf.n)
        }
        val ACb1 = seq(i).control match {
          case Some(c) => AC * kf.B.get * c
          case None => DenseVector.zeros[Double](kf.n)
        }
        val epsilon = gamma * (e1 + f.B - ACb1)
        emDtoArray(i) = EMDto(alpha, beta, gamma, delta, epsilon)
      }
      emDtoArray(end) = EMDto(alpha = fRun(end).C, beta = fRun(end).C * fRun(end).B)
      emDtoArray
    }).toList
  }

  /**
    * updated.
    * @param emDtos
    * @param prevOutputs
    * @param turn
    * @return
    */
  def updateParams(emDtos: Seq[Array[EMDto]], prevOutputs: EMOutputs, turn: TurnForEM): EMOutputs = {
    var newkf = prevOutputs.kf.copy()
    val n = newkf.n
    newkf = turn.abqTurn match {
      case 'A' =>
        val newA = optimizeCoefficientParams(cond.AOpt, getADerivative(emDtos, prevOutputs), prevOutputs.kf.A - DenseMatrix.eye[Double](n), prevOutputs.kf.Q) + DenseMatrix.eye[Double](n)
        newkf.copy(A = newA)
      case 'B' => newkf.copy(B = Some(optimizeCoefficientParams(cond.BOpt, getBDerivative(emDtos, prevOutputs),
        prevOutputs.kf.B.getOrElse(DiagMatrices(DenseVector.zeros[Double](-1))), prevOutputs.kf.Q)))
      case 'Q' => newkf.copy(Q = optimizeCovarianceParams(cond.QOpt, getQDerivative(emDtos, prevOutputs)))
      case 'F' => newkf.copy()
      case _ => newkf.copy()
    }
    newkf = turn.hrTurn match {
      case 'H' =>
        val upperIdentity = cond.HOpt._1.indexOf("upperIdentity") > 0
        val S = if (upperIdentity) {
          prevOutputs.kf.R match {
            case DenseMatrices(matrix) => DenseMatrices(matrix(n until m, n until m))
            case DiagMatrices(vector) => DiagMatrices(vector(n until m))
          }
        } else prevOutputs.kf.R
        val H2 = if (upperIdentity) DenseMatrices(prevOutputs.kf.H.toDenseMatrix(n until m, ::)) else prevOutputs.kf.H
        val newH = optimizeCoefficientParams(cond.HOpt, getHDerivative(emDtos, prevOutputs), H2, S)
        val H3 = if (upperIdentity) DenseMatrices(DenseMatrix.vertcat(DenseMatrix.eye[Double](n), newH.toDenseMatrix)) else newH
        newkf.copy(H = H3)
      case 'R' => newkf.copy(R = optimizeCovarianceParams(cond.ROpt, getRDerivative(emDtos, prevOutputs)))
      case 'F' => newkf.copy()
      case _ => newkf.copy()
    }

    var newEMOutputs = prevOutputs.copy(kf = newkf)
    newEMOutputs = turn.initTurn match {
      case 'm' =>
        val initStateMeansList = if (cond.initStateMeanOpt._1.indexOf("independent") > 0) {
          emDtos.zip(prevOutputs.initStateMeans).zip(prevOutputs.initStateCovariance).map(x => {
            val em = x._1._1
            val mean = x._1._2
            val cov = x._2
            optimizeCoefficientParams(cond.initStateMeanOpt, (em(0).beta.toDenseMatrix.t, DenseMatrix((1.0))), mean.toDenseMatrix.t, cov).toDenseMatrix.toDenseVector
          }).toList
        } else {
          //実装上の難易度の理由から、共分散は全ての系列で同一の物とする。
          require(cond.initStateCovarianceOpt._1.indexOf("independent") < 0)
          val mean = prevOutputs.initStateMeans(0)
          val cov = prevOutputs.initStateCovariance(0)
          val betas = emDtos.foldLeft(DenseVector.zeros[Double](n)) { (s, x) => s + x(0).beta }
          val optInitMean = optimizeCoefficientParams(cond.initStateMeanOpt, getInitStateMeanDerivative(emDtos, prevOutputs), mean.toDenseMatrix.t, cov).toDenseMatrix.toDenseVector
          emDtos.indices.map(x => optInitMean).toList
        }
        newEMOutputs.copy(initStateMeans = initStateMeansList)
      case 'S' =>
        val initStateCovarianceList = if (cond.initStateCovarianceOpt._1.indexOf("independent") > 0) {
          emDtos.zip(prevOutputs.initStateMeans).zip(prevOutputs.initStateCovariance).map(x => {
            val em = x._1._1
            val mean = x._1._2
            val cov = x._2
            val bm = em(0).beta - mean
            val derivs = cov match {
              case DenseMatrices(m) => (1, em(0).alpha + bm * bm.t)
              case DiagMatrices(v) => (1, (em(0).alpha + bm * bm.t).toDiagMatrices)
            }
            optimizeCovarianceParams(cond.initStateMeanOpt, derivs)
          }).toList
        } else {
          val optInitCov = optimizeCovarianceParams(cond.initStateCovarianceOpt, getInitStateCovariancesDerivative(emDtos, prevOutputs))
          emDtos.indices.map(x => optInitCov).toList
        }
        newEMOutputs.copy(initStateCovariance = initStateCovarianceList)
      case 'F' => newEMOutputs.copy()
      case _ => newEMOutputs.copy()
    }
    newEMOutputs
  }


  /**
    * optimizing other parameters than covariance.
    * S ... covariance matrix.
    * Differentiation is S^-1(J - XK).
    * smooth parameters are by ||S^-1||op * ||K||op
    *
    * @param optCond
    * @param derivatives
    * @param X
    * @param S
    * @return
    */
  def optimizeCoefficientParams(optCond: (Array[String], Array[Double]),
                                derivatives: (Matrices, Matrices),
                                X: Matrices, S: Matrices): Matrices = {
    val mode = optCond._1(0)
    val hyperParams = optCond._2
    val K = derivatives._2
    val invS = S.inv
    val J = derivatives._1

    //パラメータを何かしら入れた時の勾配
    def getGradient(M: Matrices): Matrices = {
      (invS * (J - M * K)).rev
    }
    def proximalGradient(softThreshold: (Matrices, Double) => Matrices): Matrices = {
      val initS = 1.0;
      val lambda = hyperParams(0)
      val delta = hyperParams(1)
      val time = hyperParams(2).toInt //打ち切り上限回数
      val kOp = K.operatorNorm
      if (kOp < 0.0) {
        return DenseMatrices(DenseMatrix.zeros[Double](100, 100))
      }
      val ita = kOp * invS.operatorNorm
      val initX = X
      def proximalGradientSub(prevS: Double, i: Int, prevParam: Matrices, prevProximal: Matrices): Matrices = {
        val gradient = getGradient(prevParam)
        val proximal = softThreshold(prevParam - gradient / ita, lambda / ita)
        val newS = (1 + Math.sqrt(1 + 4 * prevS * prevS)) / 2;
        val newParam = proximal + (proximal - prevProximal) * ((prevS - 1) / newS)
        val norm = (prevProximal - proximal).norm
        if (i >= time || norm < delta) {
          proximal
        } else {
          proximalGradientSub(newS, i + 1, newParam, proximal)
        }
      }
      val aho = proximalGradientSub(initS, 0, initX, initX)
      aho
    }
    val kekka: Matrices = mode match {
      case "naive" => J * K.inv
      case "l2_reg" =>
        proximalGradient((m, C) => m / (1 + C))
      case "l1_reg" =>
        proximalGradient((m, C) => {
          val fixedC = if (optCond._1.indexOf("noDataSizeFix") > 0) C else C * dataSize
          m.toDenseMatrix.map(x =>
            if (x < -fixedC) {
              x + fixedC
            } else if (x < fixedC) {
              0.0
            } else {
              x - fixedC
            }
          )
        })
      case _ => DenseMatrix.zeros[Double](-1, -1)
    }
    kekka
  }

  //共分散を最適化。微分はSymmetry(nX - L)の形を取る。
  //たまにアホみたいに大きい共分散が出ることがあり、それが続くとMatrixがSingularになったりするので、それを防ぐために補正したりする
  def optimizeCovarianceParams(optCond: (Array[String], Array[Double]), derivatives: (Int, Matrices)): Matrices = {
    val n = derivatives._1
    val L = derivatives._2
    val mode = optCond._1(0)
    val hyperParams = optCond._2
    val m = L.cols
    mode match {
      case "naive" => L / n
      case "gamma" => {
        //対角行列のみ
        val lowerBound = if (hyperParams.length < 3) 0.01 else hyperParams(2)
        val upperBound: DenseVector[Double] = if (optCond._1.indexOf("automatic") > 0) {
          allDataVariance
        } else if (hyperParams.length < 3) {
          DenseVector.fill[Double](m, 100.0)
        } else {
          DenseVector.fill[Double](m, hyperParams(3))
        }
        val dvec = DenseVector.zeros[Double](m)
        val lDiag = L.toDiagVector
        (0 until m).foreach(i => {
          val x = lDiag(i)
          val a = (x + (1 / hyperParams(1))) / (n + hyperParams(0) - 1)
          dvec(i) = if (a > upperBound(i)) upperBound(i) else if (a < lowerBound) lowerBound else a
        })
        DiagMatrices(dvec)
      }
      case _ => L / (-100000000.0)
    }
  }

  //係数行列の微分の各値を求める
  //S^{-1}(J - XK)のうち、Sは既に求められているので、JとKさえ求めれば良い。
  def getCoefficientDerivative(emDtos: Seq[Array[EMDto]], derivativeFunc: (Array[EMDto], Array[ConObs], Int) => (Matrices, Matrices), N: Int, M: Int): (Matrices, Matrices) = {
    val zeros: (Matrices, Matrices) = (DenseMatrices(DenseMatrix.zeros[Double](N, M)), DenseMatrices(DenseMatrix.zeros[Double](M, M)))
    emDtos.zip(seqs).map(x => {
      val em = x._1
      val seq = x._2
      val len = em.length
      val sumUpd = (0 until len-1).foldLeft(zeros) { (s2, i) =>
        val derivative = derivativeFunc(em, seq, i)
        (s2._1 + derivative._1, s2._2 + derivative._2)
      }
      sumUpd
    }).foldLeft(zeros) {(s, x) => (s._1 + x._1, s._2 + x._2)}
  }

  //共分散行列の微分の各値を求める
  //Sym(nS - L)となるので、nとLさえ求められれば良い。
  def getCovarianceDerivative(emDtos: Seq[Array[EMDto]], derivativeFunc: (Array[EMDto], Array[ConObs], Int) => (Int, Matrices), N: Int): (Int, Matrices) = {
    val zeros: (Int, Matrices) = (0, DiagMatrices(DenseVector.zeros[Double](N)))
    emDtos.zip(seqs).map(x => {
      val em = x._1
      val seq = x._2
      val len = em.length
      val sumUpd = (0 until len-1).foldLeft(zeros) { (s2, i) =>
        val derivative = derivativeFunc(em, seq, i)
        (s2._1 + derivative._1, s2._2 + derivative._2)
      }
      sumUpd
    }).foldLeft(zeros) { (s, x) => (s._1 + x._1, s._2 + x._2) }
  }

  def getADerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Matrices, Matrices) = {
    val n = emOutputs.kf.n
    val derivativeFunc = (em: Array[EMDto], seq: Array[ConObs], i: Int) => {
      val dbe: DenseVector[Double] = em(i).delta * em(i + 1).beta + em(i).epsilon
      val K = em(i).gamma + em(i).delta * em(i + 1).alpha * em(i).delta.t + (dbe * dbe.t)
      val bc: DenseVector[Double] = seq(i).control match {
        case Some(c) => emOutputs.kf.B.get * c;
        case None => DenseVector.zeros[Double](n)
      }
      val J = em(i + 1).alpha * em(i).delta.t + ((em(i + 1).beta - bc) * dbe.t)
      (J - K, K)
    }
    getCoefficientDerivative(emDtos, derivativeFunc, n, n)
  }

  //Noneの場合はそもそも呼び出さない
  def getBDerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Matrices, Matrices) = {
    require(l > 0)
    val n = emOutputs.kf.n
    val derivativeFunc = (em: Array[EMDto], seq: Array[ConObs], i: Int) => {
      val dbe: DenseVector[Double] = em(i).delta * em(i + 1).beta + em(i).epsilon
      val J: Matrices = seq(i).control match {
        case Some(c) => (em(i + 1).beta - emOutputs.kf.A * dbe) * c.t;
        case None => DenseMatrix.zeros[Double](n, l)
      }
      val K: Matrices = seq(i).control match {
        case Some(c) => c * c.t;
        case None => DenseMatrix.zeros[Double](l, l)
      }
      (J, K)
    }
    getCoefficientDerivative(emDtos, derivativeFunc, n, l)
  }

  def getHDerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Matrices, Matrices) = {
    val upperIdentity = cond.HOpt._1.indexOf("upperIdentity") > 0
    val n = emOutputs.kf.n
    val derivativeFunc = (em: Array[EMDto], seq: Array[ConObs], i: Int) => {
      val dbe: DenseVector[Double] = em(i).delta * em(i + 1).beta + em(i).epsilon
      val J: Matrices = seq(i).observe match {
        case Some(o) => o * em(i).beta.t
        case None => DenseMatrix.zeros[Double](m, n)
      }
      val K: Matrices = seq(i).observe match {
        case Some(o) => em(i).alpha + em(i).beta * em(i).beta.t
        case None => DenseMatrix.zeros[Double](n, n)
      }
      val J2 = if (upperIdentity) DenseMatrices(J.toDenseMatrix(n until m, ::)) else J
      (J2, K)
    }
    if (upperIdentity) getCoefficientDerivative(emDtos, derivativeFunc, m-n, n) else getCoefficientDerivative(emDtos, derivativeFunc, m, n)
  }

  //R^{-1}の微分
  def getRDerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Int, Matrices) = {
    val derivativeFunc = (em: Array[EMDto], seq: Array[ConObs], i: Int) => {
      val d: (Int, Matrices) = emOutputs.kf.R match {
        case DenseMatrices(matrix) =>
          seq(i).observe match {
            case Some(o) =>
              val hbo = emOutputs.kf.H * em(i).beta - o
              (1, emOutputs.kf.H * em(i).alpha * emOutputs.kf.H.t + hbo * hbo.t)
            case None =>
              (0, DenseMatrix.zeros[Double](m, m))
          }
        case DiagMatrices(vector) =>
          seq(i).observe match {
            case Some(o) =>
              val hbo = emOutputs.kf.H * em(i).beta - o
              (1, DiagMatrices(hbo *:* hbo + diagABAt(emOutputs.kf.H, em(i).alpha)))
            case None =>
              (0, DiagMatrices(DenseVector.zeros[Double](m)))
          }
      }
      d
    }
    getCovarianceDerivative(emDtos, derivativeFunc, m)
  }


  def getQDerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Int, Matrices) = {
    val n = emOutputs.kf.n
    val derivativeFunc = (em: Array[EMDto], seq: Array[ConObs], i: Int) => {
      val bc = seq(i).control match {
        case Some(c) => emOutputs.kf.B.get * c
        case None => DenseVector.zeros[Double](n)
      }
      val badbaebc = em(i + 1).beta - emOutputs.kf.A * em(i).delta * em(i + 1).beta - emOutputs.kf.A * em(i).epsilon - bc
      val iad = DenseMatrices(DenseMatrix.eye[Double](emOutputs.kf.n)) - emOutputs.kf.A * em(i).delta
      val d: (Int, Matrices) = emOutputs.kf.Q match {
        case DenseMatrices(matrix) => (1, emOutputs.kf.A * em(i).gamma * emOutputs.kf.A.t + badbaebc * badbaebc.t + (iad * em(i + 1).alpha * iad.t))
        case DiagMatrices(vector) => (1, DiagMatrices(diagABAt(emOutputs.kf.A, em(i).gamma) + badbaebc *:* badbaebc + diagABAt(iad, em(i + 1).alpha)))
      }
      d
    }
    getCovarianceDerivative(emDtos, derivativeFunc, n)
  }

  def getInitStateMeanDerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Matrices, Matrices) = {
    val mean = emOutputs.initStateMeans(0)
    val cov = emOutputs.initStateCovariance(0)
    val betas = emDtos.foldLeft(DenseVector.zeros[Double](emOutputs.kf.n)) { (s, x) => s + x(0).beta }
    (betas.toDenseMatrix.t, DenseMatrix((emDtos.length.toDouble)))
  }

  def getInitStateCovariancesDerivative(emDtos: Seq[Array[EMDto]], emOutputs: EMOutputs): (Int, Matrices) = {
    val cov = emOutputs.initStateCovariance(0)
    val abms = emDtos.zip(emOutputs.initStateMeans).foldLeft(DiagMatrices(DenseVector.zeros[Double](emOutputs.kf.n)): Matrices) { (s, x) =>
      val bm = x._1(0).beta - x._2
      cond.initStateCovariances(0) match {
        case DenseMatrices(m) => s + x._1(0).alpha + bm * bm.t
        case DiagMatrices(v) => s + x._1(0).alpha.toDiagMatrices + DiagMatrices(bm *:* bm)
      }
    }
    (emDtos.length, abms)
  }

  //used in diagMatrices. diagonal elements of A*B*A^t.
  private def diagABAt(A: Matrices, B: Matrices): DenseVector[Double] = sum(((A * B) *:* A).toDenseMatrix(*, ::))
}

/**
  * EMDtos.
  * @param alpha
  * @param beta
  * @param gamma
  * @param delta
  * @param epsilon
  */
case class EMDto(alpha: Matrices = null, beta: DenseVector[Double] = null,
                 gamma: Matrices = null, delta: Matrices = null,
                 epsilon: DenseVector[Double] = null)

/**
  * The parameters to be optimized in EM algorithm and the likelihood.
  * log likelihood
  *
  * @param kf
  * @param logLikelihoods
  * @param initStateMeans
  * @param initStateCovariance
  */
case class EMOutputs(kf: KalmanFilter = null, logLikelihoods: List[Double] = List(Double.MinValue),
                     initStateMeans: List[DenseVector[Double]] = List(), initStateCovariance: List[Matrices] = List())
  extends Ordered[EMOutputs] {
  //logLikelihoodを除いた部分の差のnorm
  def calcDif(emOutputs: EMOutputs): Double = {
    val kfDifDouble = Math.pow(kf.calcDif(emOutputs.kf), 2)
    val initStateMeanDifDouble = initStateMeans.zip(emOutputs.initStateMeans).foldLeft(0.0) { (s, x) => s + Math.pow(norm(x._1 - x._2), 2) }
    val initStateCovDifDouble = initStateCovariance.zip(emOutputs.initStateCovariance).foldLeft(0.0) { (s, x) => s + Math.pow((x._1 - x._2).norm, 2) }
    Math.sqrt(kfDifDouble + initStateMeanDifDouble + initStateCovDifDouble)
  }

  override def compare(that: EMOutputs): Int = {
    this.logLikelihoods.head.compare(that.logLikelihoods.head)
  }

  def initStateMeanMean: DenseVector[Double] = initStateMeans.foldLeft(DenseVector.zeros[Double](kf.n)) { (s, x) => s + x } * (1.0 / initStateMeans.length)

  def initStateCovarianceMean: Matrices = initStateCovariance(0) match {
    case DiagMatrices(vector) => DiagMatrices(initStateCovariance.foldLeft(DenseVector.zeros[Double](kf.n)) { (s, x) => s + x.toDiagVector }) / initStateMeans.length
    case DenseMatrices(matrix) => DenseMatrices(initStateCovariance.foldLeft(DenseMatrix.zeros[Double](kf.n, kf.n)) { (s, x) => s + x.toDenseMatrix }) / initStateMeans.length
  }

  override def toString() = {
    val kfSt = kf.toString()
    val covSt = initStateCovarianceMean.toString()
    val meanSt = utils.Utils.vectorTostring(initStateMeanMean)
    kfSt + "\n\n" + logLikelihoods(0) + "\n\n" + meanSt + "\n\n" + covSt
  }
}

