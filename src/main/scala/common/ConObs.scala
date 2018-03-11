package common
import breeze.linalg.DenseVector

/**
  * Created by tara on 2017/06/05.
  * this class has control and observed variables.
  * array of this class expresses sequences of all available variables
  */
case class ConObs(control: Option[DenseVector[Double]], observe: Option[DenseVector[Double]])