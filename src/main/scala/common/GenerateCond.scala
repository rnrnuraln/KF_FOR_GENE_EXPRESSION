package common

import breeze.linalg._
import breeze.stats.distributions._
import org.apache.commons.math3.random.MersenneTwister
/**
  * Created by tara on 2017/05/11.
  * 実験パラメータや初期パラメータをGenerateするためのCondition
  * sにはどうgenerateするかの情報が含まれている。
  * s = mode(必須) 何かしらの情報,fileNameまたはrandomなmodeの時のパラメータ
  * s0は,で区切られる前、s1は,で区切られた後の情報をそれぞれスペースでsplitしたもの
  * initState固定の場合は、全ての奴を同じように固定する。
  */
case class GenerateCond(s: String, m: Matrices = null, v: DenseVector[Double] = null, seed: Long = Long.MaxValue) {
  protected val s0 = s.split(",")(0).split("\\s+")
  protected val s1: Array[String] = if (s.split(",").length > 1) s.split(",")(1).split("\\s+") else Array()
  protected val s2: Array[String] = if (s.split(",").length > 2) s.split(",")(2).split("\\s+") else Array()

  val seed2 = if (seed == Long.MaxValue) System.currentTimeMillis + System.identityHashCode(this) else seed
  implicit val randBasis: RandBasis = new RandBasis(new ThreadLocalRandomGenerator(new MersenneTwister(seed2)))
  val mode = s0(0)
  val fileName: Option[String] = if (mode == "file") Some(s1(0)) else None
  val generateRand: Option[Rand[Double]] = {
    mode match {
      case "gaussian" => Some(Gaussian(s1(0).toDouble, s1(1).toDouble))
      case "gamma" => Some(Gamma(s1(0).toDouble, s1(1).toDouble))
      case "uniform" => Some(Uniform(s1(0).toDouble, s1(1).toDouble))
      case "synchro" => None
      case "identity" => None
      case "zero" => None
      case _ => None
    }
  }
  val isDiag = s0.indexOf("diag") > 0

  //initStateのconditionが各系列によってindependentかどうか
  val initStateIndependent = s0.indexOf("initStateIndependent") > 0

  val isSparse = s0.indexOf("sparse") > 0

  def toSparse(mat: DenseMatrix[Double]): Matrices = {
    val n = mat.rows
    val m = mat.cols
    val r = DenseMatrix.rand[Double](n, m)
    for (i <- 0 until n) {
      for (j <- 0 until m) {
        mat(i, j) = if (r(i, j) > s2(0).toDouble) {
          mat(i, j)
        } else {
          0.0
        }
      }
    }
    DenseMatrices(mat)
  }

  //モデルから曖昧性をなくすため、行列の上の部分をidentity matrixとする(Hに用いる)
  val upperIdentity = s0.indexOf("upperIdentity") > 0

  val matrix = readMatrix(fileName).getOrElse(m)
  val vector = readVector(fileName).getOrElse(v)

  def abstractMakeMatrix(n: Int, m: Int, f: Rand[Double] => (Matrices, String)): Matrices = {
    generateRand match {
      case Some(x) => {
        val matrices = f(x)._1
        val s = f(x)._2
        //sparse性を考慮する
        val realMatrices = if (!isSparse || s == "A") matrices else toSparse(matrices.toDenseMatrix)
        if (realMatrices.max > 10.0) abstractMakeMatrix(n, m, f) else realMatrices
      }
      case None =>
        mode match {
          case "synchro" => matrix
          case "identity" => DiagMatrices(DenseVector.ones[Double](n))
          case "zero" => DenseMatrices(DenseMatrix.zeros[Double](n, m))
          case _ => DenseMatrices(DenseMatrix.zeros[Double](n, m))
        }
    }
  }

  //diagMatrixを作る場合は、ほとんどの場合逆ガンマ分布を仮定していると思われるので逆数にする(逆ガンマ分布の代わりみたいな感じ)
  def makeDiagMatrix(n: Int) = {
    require(!isSparse, "diag matrix must not be sparse")
    abstractMakeMatrix(n, n, x => {
      def makeA: DiagMatrices = {
        val a = DenseVector.rand[Double](n, x).map(y => 1.0 / y)
        val aMin = a.foldLeft(Double.MaxValue) { (s, y) => if (s > y) y else s }
        val aMax = a.foldLeft(-1.0) { (s, y) => if (s > y) y else s }
        //極端な値の排除
        if (aMin < 0.01 || aMax > 10.0) makeA else DiagMatrices(a)
      }
      (makeA, "diag")
    })
  }

  def makeDenseMatrix(n: Int, m: Int, s: String = "dense") = abstractMakeMatrix(n, m, x => (DenseMatrices(DenseMatrix.rand[Double](n, m, x)), s))

  def makeMatrixForA(n: Int) = {
    def makeA(rand: Rand[Double]): (Matrices, String) = {
      def makeP: DenseMatrix[Double] = {
        val P = DenseMatrix.rand(n, n, rand)
        if (det(P) > 0.01) P else makeP
      }
      def makepI: DenseMatrix[Double] = {
        val prevP = DenseMatrix.rand(n, n, rand)
        val P = if (isSparse) toSparse(prevP).toDenseMatrix else prevP
        val pI = P + DenseMatrix.eye[Double](n)
        if (det(pI) > 0.01) pI else makeP
      }
      //lambdaFirstの場合は、lambdaを先に生成してから後で調節する。こっちではsparseはdo not consider
      val A = if (s0.indexOf("lambdaFirst") > 0) {
        val P = makeP
        //固有値が-1 ~ 1の値を取るように調節
        val lambda = diag(DenseVector.rand[Double](n, Uniform(0, 1)).map(x => 2 * Math.abs(x) - 1))
        inv(P) * lambda * P
      } else {
        //Aの生成方法について、固有値を先に計算するのではなく、適当に生成してから固有値の絶対値の最大値を1に揃える方法を用いる。
        //固有値が虚数の場合とか色々考えることができるはず。
        //仮生成した行列PにIを加えた行列の固有値の絶対値の最大値を求める
        val pI = makepI
        val eigpI = eig(pI)
        //各固有値の絶対値を求める。虚部まで含めた時の値。
        val pIabs: DenseVector[Double] = ((eigpI.eigenvalues *:* eigpI.eigenvalues) +
          (eigpI.eigenvaluesComplex *:* eigpI.eigenvaluesComplex)).map(x => Math.sqrt(x))
        val pMax = pIabs.foldLeft(Double.MinValue) { (s, x) => if (s > x) s else x }
        pI / pMax
      }
      //値が変に大きくなりすぎないよう制限
      if (breeze.numerics.abs(max(A)) > 5.0) makeA(rand) else (DenseMatrices(A), "A")
    }
    abstractMakeMatrix(n, n, makeA)
  }

  //本当はこの実装だと良くないんだけど、これで行ってしまう
  def makeMatrixForB(n: Int, m: Int): Option[Matrices] = {
    if (m == 0) None else Some(makeDenseMatrix(n, m, "B"))
  }

  def makeMatrixForH(m: Int, n: Int): Matrices = {
    if (upperIdentity) {
      if (m > n) {
        val lower = makeDenseMatrix(m - n, n)
        DenseMatrices(DenseMatrix.vertcat(DenseMatrix.eye[Double](n), lower.toDenseMatrix))
      } else {
        DenseMatrices(DenseMatrix.eye[Double](m))
      }
    } else {
      makeDenseMatrix(m, n, "H")
    }
  }

  //covarianceのために正定値対称なmatrixを生成 -> のはずだったが、diagMatrixを作るだけにした。それ以外は考えない
  def makeMatrixForCov(n: Int): Matrices = {
    makeDiagMatrix(n)
    /**
    else if (mode == "synchronized") matrix
    else {
      def makeCov(rand: Rand[Double]): Matrices = {
        def makePrev: DenseMatrix[Double] = {
          val P = DenseMatrix.rand(n, n, rand)
          if (det(P) > 0.01) P else makePrev
        }
        //まずランダムに行列を生成
        val prevPrevCov = makePrev
        //対称化
        val prevCov = prevPrevCov + prevPrevCov.t
        //固有値のベクトルの最小値を取って来る
        val l: Double = eig(prevCov).eigenvalues.map(_ + 1).foldLeft(0.0) {
          (s, x) => if (x > s) s else x
        }
        //固有値が0~の値を取るように調節
        DenseMatrices(if (l < 0) prevCov + DenseMatrix.eye[Double](n) * (l + Math.abs(rand.sample())) else prevCov)
      }
      abstractMakeMatrix(n, n, makeCov)
    }
      */
  }

  protected def readMatrix(fileName: Option[String]): Option[Matrices] = {
    fileName match {
      case Some(f) =>
        val s = scala.io.Source.fromFile(f).getLines()
        val list = s.map(s => s.split("\t").map(_.toDouble)).toList
        Some(DenseMatrices(DenseMatrix(list: _*)))
      case None =>
        None
    }
  }

  protected def readVector(fileName: Option[String]): Option[DenseVector[Double]] = {
    fileName match {
      case Some(f) =>
        val s = scala.io.Source.fromFile(f).getLines()
        if (s.hasNext) {
          val d = DenseVector(s.next().split("\t").map(_.toDouble): _*)
          Some(d)
        } else {
          println("nothing written on vector file!")
          None
        }
      case None => None
    }
  }

  /**
    * control variableはこれで大体生成する
    *
    * @param n
    * @return
    */
  def makeVector(n: Int): DenseVector[Double] = {
    generateRand match {
      case Some(x) => DenseVector.rand(n, x)
      case None =>
        if (mode == "synchro") {
          vector
        } else {
          null
        }
    }
  }

  /**
    * initStateの平均のリストを作成する
    *
    * @param n
    * @param listLen
    * @return
    */
  def makeInitStateMeanList(n: Int, listLen: Int): List[DenseVector[Double]] =
    makeInitStates[DenseVector[Double]](n, listLen, makeVector)

  /**
    * initStateのcovarianceのリストを作成する
    */
  def makeInitStateCovarianceList(n: Int, listLen: Int): List[Matrices] =
    makeInitStates[Matrices](n, listLen, makeMatrixForCov)

  /**
    * initStateを作る抽象的な関数
    *
    * @param n
    * @param listLen
    * @param function
    * @tparam T
    * @return
    */
  def makeInitStates[T](n: Int, listLen: Int, function: (Int => T)): List[T] = {
    if (initStateIndependent) {
      (0 until listLen).foldLeft(List(): List[T]) { (s, x) =>
        function(n) :: s
      }.reverse
    } else {
      val madeT = function(n)
      (0 until listLen).foldLeft(List(): List[T]) { (s, x) =>
        madeT :: s
      }.reverse
    }
  }

  //control項の値が前後のそれに依存するという仮定の元、その系列を生成するのでございます。
  //具体的には、ある地点での値に次の値をrandとして加える。
  //初期値は他のrandと同じ生成の仕方だが、その他の項に関してはs1(2)の値に応じて大きくしたり小さくしたりする
  def makeControlSeq(n: Int, seqLen: Int): Array[Option[DenseVector[Double]]] = {
    if (n == 0) Array.fill(seqLen)(None) else {
      val init = makeVector(n)
      val array: Array[Option[DenseVector[Double]]] = Array.fill(seqLen)(None)
      array(0) = Some(init)
      for (i <- 1 until seqLen) {
        val fluctuation = if (s1.length < 3) makeVector(n) else makeVector(n) * s1(2).toDouble
        array(i) = Some(array(i-1).get + fluctuation)
      }
      array
    }
  }
}