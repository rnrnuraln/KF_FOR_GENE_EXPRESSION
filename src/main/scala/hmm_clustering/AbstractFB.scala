package hmm_clustering

/**
  * Created by tara on 2019/09/29.
  */
trait AbstractFB {

  //Kalman Filter used for calculation of each algorithm
  val hmm: HMM
}
