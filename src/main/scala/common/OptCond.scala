package common

/**
  * Created by tara on 2018/07/30.
  */
case class OptCond(s: String) {
  protected val s0 = s.split(",")(0).split("\\s+")
  //正則化に用いるパラメータたち
  val regParams: Array[Double] = if (s.split(",").length > 1) s.split(",")(1).split("\\s+").map(_.toDouble) else Array(1.0)
  //その他の最適化等に用いるパラメータたち
  val subParams: Array[Double] = if (s.split(",").length > 2) s.split(",")(2).split("\\s+").map(_.toDouble) else Array()

  val mode = s0(0)

  var lambda = regParams(0)

  //モデルから曖昧性をなくすため、行列の上の部分をidentity matrixとする(Hに用いる)
  val upperIdentity = s0.indexOf("upperIdentity") > 0

  val isIndependent = s0.indexOf("independent") > 0

  val noDataSizeFix = s0.indexOf("noDataSizeFix") > 0

  val automatic = s0.indexOf("automatic") > 0

}

