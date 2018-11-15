package common

import breeze.linalg.{MatrixSingularException, _}

/**
  * This class is the wrapper of Matrix class used in this program.
  * There are two classes now; DenseMatrices and DiagMatrices.
  * Maybe we will use breeze.linalg.CSCMatrix, in order to handle with sparse matrix efficiently.
  * Created by tara on 2017/06/28.
  */
abstract class Matrices {
  val cols: Int
  val rows: Int

  def toArray: Array[Double]

  def *(matrices: Matrices) = matrixCaseFunc(matrices, timeDense, timeDiag)

  protected def timeDense(matrix: DenseMatrix[Double]): Matrices

  protected def timeDiag(vector: DenseVector[Double]): Matrices

  def +(matrices: Matrices) = matrixCaseFunc(matrices, plusDense, plusDiag)

  protected def plusDense(matrix: DenseMatrix[Double]): Matrices

  protected def plusDiag(vector: DenseVector[Double]): Matrices

  def *(vector: DenseVector[Double]): DenseVector[Double]

  def *(d: Double): Matrices

  def /(d: Double): Matrices

  def t: Matrices

  def inv: Matrices

  //reverse of the sign of matrix
  def rev: Matrices

  def -(matrices: Matrices): Matrices = this + matrices.rev

  private def abstractCaseFunc[U](u: U, denseFunc: (DenseMatrix[Double] => U), diagFunc: (DenseVector[Double] => U)): U = {
    u match {
      case DenseMatrices(x) => denseFunc(x)
      case DiagMatrices(x) => diagFunc(x)
    }
  }

  def matrixCaseFunc(matrices: Matrices, denseFunc: (DenseMatrix[Double] => Matrices), diagFunc: (DenseVector[Double] => Matrices)): Matrices =
    abstractCaseFunc[Matrices](matrices, denseFunc, diagFunc)

  //see maximum relative error of another matrices A against this matrix
  def maxDif(A: Matrices, fixDelta: Double): Double = {
    utils.Utils.maxDif(this.toArray, A.toArray, fixDelta)
  }

  def norm: Double

  def toDenseMatrix: DenseMatrix[Double]

  //take diagonal elements of the matrices
  def toDiagVector: DenseVector[Double]

  def operatorNorm: Double

  def \(m: Matrices): Matrices

  def *:*(m: Matrices): Matrices

  def det: Double

  //leave diagonal elements of the matrices and take off another
  def toDiagMatrices: DiagMatrices = DiagMatrices(toDiagVector)

  def trace: Double

  def max: Double
}

/**
  * DenseMatrices includes breeze.linalg.DenseMatrix as constructor variable
  * Usually uses this class
  * @param denseMatrix
  */
case class DenseMatrices(denseMatrix: DenseMatrix[Double] = null) extends Matrices {
  val cols = denseMatrix.cols
  val rows = denseMatrix.rows

  def toArray = denseMatrix.toArray

  def timeDense(matrix: DenseMatrix[Double]) = DenseMatrices(denseMatrix * matrix)

  def timeDiag(vector: DenseVector[Double]) = {
    val matrix = denseMatrix.copy
    val len = vector.length
    require(denseMatrix.cols == len)
    for (j <- 0 until denseMatrix.cols) {
      matrix(::, j) := matrix(::, j) * vector(j)
    }
    DenseMatrices(matrix)
  }

  def plusDense(matrix: DenseMatrix[Double]) = DenseMatrices(denseMatrix + matrix)

  def plusDiag(vector: DenseVector[Double]) = DenseMatrices(denseMatrix + breeze.linalg.diag(vector))

  def *(vector: DenseVector[Double]) = denseMatrix * vector

  def inv = {
    try {
      DenseMatrices(breeze.linalg.inv(denseMatrix))
    } catch {
      case m: MatrixSingularException => {
        DenseMatrices(breeze.linalg.pinv(denseMatrix));
      }
    }
  }

  def t = DenseMatrices(denseMatrix.t)

  def rev = DenseMatrices(-denseMatrix)

  override def toString: String = {
    val s = new StringBuilder()
    for (i <- 0 until denseMatrix.rows) {
      for (j <- 0 until denseMatrix.cols - 1) {
        s.append(denseMatrix(i, j) + "\t")
      }
      s.append(denseMatrix(i, denseMatrix.cols - 1) + "\n")
    }
    s.toString
  }

  def norm = Math.sqrt(denseMatrix.data.foldLeft(0.0) { (s, x) => s + x * x })

  override def toDenseMatrix: DenseMatrix[Double] = denseMatrix

  override def toDiagVector: DenseVector[Double] = {
    require(denseMatrix.cols == denseMatrix.rows)
    DenseVector((for {i <- 0 until denseMatrix.cols} yield denseMatrix(i, i)).toArray)
  }

  def *(d: Double) = DenseMatrices(denseMatrix * d)

  def /(d: Double) = DenseMatrices(denseMatrix / d)

  def operatorNorm: Double = {
    if (denseMatrix.valuesIterator.exists(_.isNaN)) {
      println("m: " + denseMatrix + "end\n")
      -1.0
    } else {
      eig(denseMatrix).eigenvalues.foldLeft(Double.MinValue) { (s, x) => if (s > x) s else x }
    }
  }

  def \(m: Matrices): Matrices = DenseMatrices(denseMatrix \ m.toDenseMatrix)

  def *:*(m: Matrices): Matrices = {
    m match {
      case DenseMatrices(matrix) => DenseMatrices(denseMatrix *:* matrix)
      case DiagMatrices(vector) => DiagMatrices(this.toDiagVector *:* vector)
    }
  }

  def det = breeze.linalg.det(denseMatrix)

  def trace = breeze.linalg.trace(denseMatrix)

  def max = breeze.linalg.max(denseMatrix)

}

/**
  * DiagMatrices includes Densevector.
  * This handles only diagonal matrix.
  * Mainly used for Covariance Matrix.
  * @param diagVector
  */
case class DiagMatrices(diagVector: DenseVector[Double] = null) extends Matrices {
  val rows = diagVector.length
  val cols = diagVector.length

  def toArray = diagVector.toArray

  def timeDense(matrix: DenseMatrix[Double]) = {
    val matrix2 = matrix.copy
    val len = diagVector.length
    require(matrix.rows == len)
    for (i <- 0 until matrix.rows) {
      matrix2(i, ::) := matrix(i, ::) * diagVector(i)
    }
    DenseMatrices(matrix2)
  }

  def timeDiag(vector: DenseVector[Double]) = DiagMatrices(diagVector *:* vector)

  def plusDense(matrix: DenseMatrix[Double]) = DenseMatrices(matrix) + this

  def plusDiag(vector: DenseVector[Double]) = DiagMatrices(diagVector + vector)

  def *(vector: DenseVector[Double]) = diagVector *:* vector

  def t = this

  def inv = DiagMatrices(diagVector.map(x => 1.0 / x))

  def rev = DiagMatrices(-diagVector)

  override def toString: String = {
    utils.Utils.vectorTostring(diagVector) + "\n"
  }

  def norm = breeze.linalg.norm(diagVector)

  override def toDenseMatrix: DenseMatrix[Double] = diag(diagVector)

  override def toDiagVector: DenseVector[Double] = diagVector

  def *(d: Double) = DiagMatrices(diagVector * d)

  def /(d: Double) = DiagMatrices(diagVector / d)

  def operatorNorm: Double = diagVector.foldLeft(Double.MinValue) { (s, x) => if (s > x) s else x }

  def \(m: Matrices): Matrices = DenseMatrices(diag(diagVector) \ m.toDenseMatrix)

  def *:*(m: Matrices): Matrices = {
    DiagMatrices(diagVector *:* m.toDiagVector)
  }

  def det = diagVector.foldLeft(1.0) { (s, x) => s * x }

  def trace = sum(diagVector)

  def max = breeze.linalg.max(diagVector)
}

