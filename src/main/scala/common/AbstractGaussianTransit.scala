package common

import breeze.linalg.{DenseVector, Matrix}
import utils.Utils

trait AbstractGaussianTransit

case class ForwardGaussian(B: DenseVector[Double] = null, C: Matrices = null) extends AbstractGaussianTransit

case class BackwardGaussian(D: Matrices = null, E: DenseVector[Double] = null) extends AbstractGaussianTransit
