package utils
import java.io.PrintWriter

import breeze.linalg._
import common.{DenseMatrices, DiagMatrices}
import common.ConObs

import scala.collection.immutable.Map
import scala.io.Source

object Utils {
  def readInitState(fileName: String) = {
    val source = scala.io.Source.fromFile(fileName).getLines()
    val x = source.next().split("\t").map(x => x.toDouble)
    val xVec = DenseVector(x:_*)
    xVec
  }
  def readObserve(fileName: String) = {
    val source = scala.io.Source.fromFile(fileName)
    source.getLines.map(x => DenseVector(x.split("\t").map(_.toDouble):_*)).toArray
  }
  def readControl(fileName: String): Option[Array[DenseVector[Double]]] = {
    val source = scala.io.Source.fromFile(fileName)
    if (source.hasNext) {
      Some(source.getLines.map(x => DenseVector(x.split("\t").map(_.toDouble): _*)).toArray)
    } else {
      None
    }
  }
  def vectorTostring(vector: Vector[Double]): String = {
    val n = vector.length
    val s = new StringBuilder
    for (i <- 0 until n-1) {
      s.append(vector(i).toString + "\t")
    }
    s.append(vector(n-1))
    s.toString()
  }
  def makePrintWriter(output: String): PrintWriter = {
    if (output == "stdout" || output == "") {
      new PrintWriter(System.out)
    } else {
      new PrintWriter(output)
    }
  }

  def makeSeqs(input: String): List[Array[ConObs]] = {
    val source = scala.io.Source.fromFile(input).getLines()
    val header = source.next().split("\t")
    val m = header.count(x => x.startsWith("obs"))
    val l = header.count(x => x.startsWith("con"))
    //必ずconの方から始まる奴を想定している。改良したいですね。
    if (l > 0) {
      require(header(0) == "con1")
    }

    //madeListにはすでに終わった系列、makingListには作り途中の系列があるぞい
    //入力はcon -> obsの順 系列は1-origin
    def makeSeqsSub(madeList: List[Array[ConObs]], makingList: List[ConObs]): List[Array[ConObs]] = {
      if (source.hasNext) {
        val s = source.next()
        //何もない場合はそのまま次に移行
        if (s == "") {
          makeSeqsSub(madeList, makingList)
        } else {
          val ss = s.split("\t")
          val con = if (l > 0) Some(DenseVector[Double](ss.slice(0, l).map(x => x.toDouble))) else None
          val obs = if (ss(l) == "") None else Some(new DenseVector[Double](ss.slice(l, l + m).map(x => x.toDouble)))
          //new sequence
          if (ss(m + l).toInt > madeList.length + 1) {
            makeSeqsSub(makingList.reverse.toArray :: madeList, List(ConObs(con, obs)))
          } else {
            makeSeqsSub(madeList, ConObs(con, obs) :: makingList)
          }
        }
      } else {
        //終了処理
        (makingList.reverse.toArray :: madeList).reverse
      }
    }
    makeSeqsSub(List(), List())
  }

  /**
    * calculate max(||v - w||)
    * @param v
    * @param w
    * @param fixDelta
    * @return
    */
  def maxDif(v: Array[Double], w: Array[Double], fixDelta: Double): Double = {
    v.zip(w).foldLeft(0.0) {(max, x) =>
      val dif = Math.abs(x._1 - x._2) / (Math.abs(x._1) + fixDelta)
      if (dif > max) dif else max
    }
  }

}

