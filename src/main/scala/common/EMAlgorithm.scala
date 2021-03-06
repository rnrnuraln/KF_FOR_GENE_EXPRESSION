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

  val allDataMean = seqs.foldLeft(DenseVector.zeros[Double](m)) {(s, x) =>
    s + x.foldLeft(DenseVector.zeros[Double](m)) {(s2, y) =>
      s2 + y.observe.getOrElse(DenseVector.zeros[Double](m))
    }
  } * (1.0 / dataSize)


  val allDataVariance = seqs.foldLeft(DenseVector.zeros[Double](m)) {(s, x) =>
    s + x.foldLeft(DenseVector.zeros[Double](m)) {(s2, y) =>
      y.observe match {
        case Some(o) => s2 + (o - allDataMean) *:* (o - allDataMean)
        case None => s2
      }
    }
  } * (1.0 / dataSize)

  val allDataStd = allDataVariance.map(x => Math.sqrt(x))

  val regularizedSeqs: List[Array[common.ConObs]] = seqs.map(x => {
    x.map(y => {
      y.observe match {
        case Some(o) => ConObs(y.control, Some((o - allDataMean) /:/ allDataStd))
        case None => ConObs(y.control, y.observe)
      }
    })
  })

  val trainSeqs: List[Array[common.ConObs]] = if (cond.seqRegularization) regularizedSeqs else seqs


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
    val Afixed = cond.AOpt.mode == "fixed"
    val Bfixed = cond.BOpt.mode == "fixed" || l == 0
    val Qfixed = cond.QOpt.mode == "fixed"
    val Hfixed = cond.HOpt.mode == "fixed" || (cond.HOpt.upperIdentity && (m == n))
    val Rfixed = cond.ROpt.mode == "fixed"
    val initStateMeanfixed = cond.initStateMeanOpt.mode == "fixed"
    val initStateCovariancefixed = cond.initStateCovarianceOpt.mode == "fixed"
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
    val stopCalcVal = 1000.0// 正則化項こそあれ、あまりに大きい項が出てきた場合はそこで打ち止め。ここではとりあえず1000以上が出てきたらやめる方針で
    def learnSub(prevOutputs: EMOutputs, i: Int, turn: TurnForEM, forwards: Seq[Forward]): EMOutputs = {
      val trainSeqsWithInits = trainSeqs.zip(prevOutputs.initStateMeans.zip(prevOutputs.initStateCovariance))
      // Kalman smoothing using forward and backward algorithm
      val backwards = trainSeqsWithInits.map(x => Backward(prevOutputs.kf, x._1, x._2._1, x._2._2))
      val emDtos = obtainEMDto(forwards.map(_.run).zip(backwards.map(_.run)), prevOutputs.kf)
      val outputs = updateParams(emDtos, prevOutputs, turn)
      val newForwards = trainSeqsWithInits.map(x => Forward(outputs.kf, x._1, x._2._1, x._2._2))
      val logLikelihood = newForwards.map(_.logLikelihood).sum
      val newOutputs = outputs.copy(logLikelihoods = logLikelihood :: prevOutputs.logLikelihoods)
      if (logLikelihood.isNaN || logLikelihood.isInfinity || logLikelihood + 10.0 < prevOutputs.logLikelihoods(0) || outputs.kf.max > stopCalcVal) {
        //あり得ない結果
        prevOutputs
      } else if (i > cond.emTime || prevOutputs.calcDif(outputs) < cond.delta) {
        //終了条件
        newOutputs
      } else {
        learnSub(newOutputs, i + 1, turn.rotate(), newForwards)
      }
    }
    //どのturnから始めるかをrandomとして設定するためのメソッド
    def initTurn(i: Int, N: Int, turn: TurnForEM): TurnForEM = {
      if (i < N) initTurn(i+1, N, turn.rotate()) else turn
    }
    val turn = initTurn(0, (Math.random() * 6).toInt, new TurnForEM().rotate())
    val mean = if (cond.seqRegularization) allDataMean else DenseVector.zeros[Double](m)
    val std = if (cond.seqRegularization) allDataStd else DenseVector.ones[Double](m)
    //rのregularization
    val r = if (cond.seqRegularization && cond.RRegularization) cond.initkf.R * DiagMatrices(allDataVariance).inv else cond.initkf.R
    val kf = cond.initkf.copy(R = r)
    val trainSeqsWithInits = trainSeqs.zip(cond.initStateMeans.zip(cond.initStateCovariances))
    val forwards = trainSeqsWithInits.map(x => Forward(kf, x._1, x._2._1, x._2._2))
    learnSub(EMOutputs(kf, List(Double.MinValue), cond.initStateMeans, cond.initStateCovariances, mean, std), 0, turn, forwards)
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
    fbSeq.zip(trainSeqs).map(x => {
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
        val upperIdentity = cond.HOpt.upperIdentity
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
        val initStateMeansList = if (cond.initStateMeanOpt.isIndependent) {
          emDtos.zip(prevOutputs.initStateMeans).zip(prevOutputs.initStateCovariance).map(x => {
            val em = x._1._1
            val mean = x._1._2
            val cov = x._2
            optimizeCoefficientParams(cond.initStateMeanOpt, (em(0).beta.toDenseMatrix.t, DenseMatrix((1.0))), mean.toDenseMatrix.t, cov).toDenseMatrix.toDenseVector
          }).toList
        } else {
          //実装上の難易度の理由から、共分散は全ての系列で同一の物とする。
          require(!cond.initStateCovarianceOpt.isIndependent)
          val mean = prevOutputs.initStateMeans(0)
          val cov = prevOutputs.initStateCovariance(0)
          val betas = emDtos.foldLeft(DenseVector.zeros[Double](n)) { (s, x) => s + x(0).beta }
          val optInitMean = optimizeCoefficientParams(cond.initStateMeanOpt, getInitStateMeanDerivative(emDtos, prevOutputs), mean.toDenseMatrix.t, cov).toDenseMatrix.toDenseVector
          emDtos.indices.map(x => optInitMean).toList
        }
        newEMOutputs.copy(initStateMeans = initStateMeansList)
      case 'S' =>
        val initStateCovarianceList = if (cond.initStateCovarianceOpt.isIndependent) {
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
    * smooth parameters are calculated by ||S^-1||op * ||K||op
    *
    * @param optCond
    * @param derivatives
    * @param X
    * @param S
    * @return
    */
  def optimizeCoefficientParams(optCond: OptCond,
                                derivatives: (Matrices, Matrices),
                                X: Matrices, S: Matrices): Matrices = {
    val mode = optCond.mode
    val K = derivatives._2
    val invS = S.inv
    val J = derivatives._1

    //パラメータを何かしら入れた時の勾配
    def getGradient(M: Matrices): Matrices = {
      (invS * (J - M * K)).rev
    }
    def proximalGradient(softThreshold: (Matrices, Double) => Matrices): Matrices = {
      val initS = 1.0;
      val lambda = optCond.lambda
      val delta = optCond.subParams(0)
      val time = optCond.subParams(1).toInt //打ち切り上限回数
      val kOp = K.operatorNorm
      if (kOp < 0.0) {
        print("K operator norm is minus" + K)
        return DenseMatrices(DenseMatrix.zeros[Double](J.rows, K.rows))
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
          val fixedC = if (optCond.noDataSizeFix) C else C * dataSize
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
  def optimizeCovarianceParams(optCond: OptCond, derivatives: (Int, Matrices)): Matrices = {
    val n = derivatives._1
    val L = derivatives._2
    val mode = optCond.mode
    val m = L.cols
    mode match {
      case "naive" => L / n
      case "gamma" => {
        //対角行列のみ
        val lowerBound = if (optCond.subParams.length < 1) 0.01 else optCond.subParams(0)
        val upperBound: DenseVector[Double] = if (optCond.automatic) {
          allDataVariance
        } else if (optCond.subParams.length < 2) {
          DenseVector.fill[Double](m, 100.0)
        } else {
          DenseVector.fill[Double](m, optCond.subParams(1))
        }
        val dvec = DenseVector.zeros[Double](m)
        val lDiag = L.toDiagVector
        (0 until m).foreach(i => {
          val x = lDiag(i)
          val a = (x + (1 / optCond.regParams(1))) / (n + optCond.regParams(0) - 1)
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
    emDtos.zip(trainSeqs).map(x => {
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
    emDtos.zip(trainSeqs).map(x => {
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
    val upperIdentity = cond.HOpt.upperIdentity
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
  *
  */
case class EMOutputs(kf: KalmanFilter = null, logLikelihoods: List[Double] = List(Double.MinValue),
                     initStateMeans: List[DenseVector[Double]] = List(), initStateCovariance: List[Matrices] = List(),
                     allMean: DenseVector[Double] = DenseVector(), allStd: DenseVector[Double]  = DenseVector[Double]())
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
    val allMeanSt = utils.Utils.vectorTostring(allMean)
    val allVarianceSt = utils.Utils.vectorTostring(allStd)
    kfSt + "\n\n" + logLikelihoods(0) + "\n\n" + meanSt + "\n\n" + covSt + "\n\n" + allMeanSt + "\n\n" + allVarianceSt + "\n"
  }
}

object EMOutputs {
  //外部からファイルによるinputでEMOutputを作成 skip行飛ばす
  def apply(input: String, skip: Int): EMOutputs = {
    val source = scala.io.Source.fromFile(input).getLines()
    (0 until skip).foreach(x => source.next()) //n回skip
    val n = source.next().toInt
    val m = source.next().toInt
    val l = source.next().toInt
    source.next()
    //対角行列の時の対処がめんどくさいですね……
    def readMatrices(matrixList: List[Array[Double]], rowNum: Int): Matrices = {
      val s = source.next()
      if (s != "") {
        val ss = s.split("\t").map(_.toDouble)
        readMatrices(ss :: matrixList, rowNum)
      } else {
        if (matrixList.length == rowNum) {
          DenseMatrices(DenseMatrix(matrixList.reverse: _*))
        } else {
          //対角行列の場合
          DiagMatrices(DenseVector(matrixList(0): _*))
        }
      }
    }
    def readVector(): DenseVector[Double] = {
      val vector = source.next().split("\t").map(_.toDouble)
      if (source.hasNext) {
        source.next()
      }
      DenseVector(vector: _*)
    }
    val A = readMatrices(List(), n)
    val B = if (l == 0) {source.next(); None} else Some(readMatrices(List(), n))
    val H = readMatrices(List(), m)
    val Q = readMatrices(List(), n)
    val R = readMatrices(List(), m)
    val kf = KalmanFilter(A, B, H, Q, R)
    source.next()
    val logLikelihood = source.next().toDouble
    source.next()
    val mu = readVector()
    val S = readMatrices(List(), n)
    source.next()
    val mean = readVector()
    val std = readVector()
    EMOutputs(kf, List(logLikelihood), List(mu), List(S), mean, std)
  }

  def readLogLikelihoods(input: String): Unit = {
    val em = EMOutputs(input, 0)
    println(em.logLikelihoods(0))
  }
}