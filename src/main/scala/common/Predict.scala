package common
import breeze.linalg._
import utils.Utils

/**
  * Created by tara on 2017/07/21.
  * カルマンフィルタモデルを用いた隠れ状態や観測状態の推定
  */
class Predict(emOutputs: EMOutputs, seqs: Seq[Array[ConObs]], predictCond: PredictCond) {

  /**
    * 多分ここに隠れ状態とかを予測っつーか推定っつーかするメソッドを作る
    */
  def run: Seq[Array[(DenseVector[Double], DenseVector[Double])]] = {
    predictCond.mode match {
      case "obs" => predictObs(predictCond.forwardEstimateNum(0))
      case "hid" => predictHid()
      case _ => Seq()
    }
  }

  //n個先を共分散と共に予測。ただし、共分散求めるのは大変なので、分散only
  def predictObs(n: Int): Seq[Array[(DenseVector[Double], DenseVector[Double])]] = {
    seqs.map(x => {
      val regularizedSeq = x.map(y => {
          y.observe match {
            case Some(o) => ConObs(y.control, Some((o - emOutputs.allMean) /:/ emOutputs.allStd))
            case None => ConObs(y.control, y.observe)
          }
      })
      val fw = Forward(emOutputs.kf, regularizedSeq, emOutputs.initStateMeanMean, emOutputs.initStateCovarianceMean)
      val predicted = fw.predictObs(n)
      predicted.map(y => {
        val trueMean = y._1 *:* emOutputs.allStd + emOutputs.allMean
        val trueVar = y._2 *:* (emOutputs.allStd *:* emOutputs.allStd)
        (trueMean, trueVar)
      })
    })
  }

  //カルマンフィルタリングによる内部状態の推定
  def predictHid(): Seq[Array[(DenseVector[Double], DenseVector[Double])]] = {
    seqs.map(x => {
      val regularizedSeq = x.map(y => {
        y.observe match {
          case Some(o) => ConObs(y.control, Some((o - emOutputs.allMean) /:/ emOutputs.allStd))
          case None => ConObs(y.control, y.observe)
        }
      })
      val fw = Forward(emOutputs.kf, regularizedSeq, emOutputs.initStateMeanMean, emOutputs.initStateCovarianceMean)
      fw.estimateHid()
    })
  }

  //カルマンスムーシングによる内部状態の推定


}

object Predict {
  def predict(emOutputFile: String, seqFile: String, predictCondFile: String, outputFile: String): Unit = {
    val emOutputs = EMOutputs(emOutputFile, 0)
    val seqs = utils.Utils.readSeqs(seqFile)
    val predictCond = PredictCond(predictCondFile)
    val predict = new Predict(emOutputs, seqs, predictCond)
    val predictSeqs = predict.run
    writeSeqs(predictSeqs, outputFile)
  }

  def writeSeqs(seqs: Seq[Array[(DenseVector[Double], DenseVector[Double])]], outputFile: String): Unit = {
    val sb = new StringBuilder()
    seqs.zipWithIndex.foreach(x => {
      val seqNumber = x._2
      x._1.foreach(y => {
        sb.append(utils.Utils.vectorTostring(y._1) + "\t" + utils.Utils.vectorTostring(y._2) + "\t" + (seqNumber + 1) + "\n")
      })
    })
    val fp = Utils.makePrintWriter(outputFile)
    fp.print(sb.toString())
    fp.close()
  }

}

