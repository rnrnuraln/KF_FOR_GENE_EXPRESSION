package analysis

import common._
import utils._
import breeze.linalg._

/**
  * Created by tara on 2018/08/23.
  * 適当に置いておく的な
  */
class Analyze {



}

object Analyze {

  //predictionにまつわるanalysis
  def analyzePredict(seqFile: String, logFile: String, skip: Int, dif: Int, foldNum: Int): ((PredictCond, String) => Unit) = {
    val seqs = utils.Utils.readSeqs(seqFile)
    val seqNum = seqs.length
    val emOutputSeq = (0 until foldNum).map(i => EMOutputs(logFile, skip + i * dif))

    (predCond: PredictCond, output: String) => {
      val predicted = (0 until foldNum).foldLeft(List(): List[Array[(DenseVector[Double], DenseVector[Double])]]) { (s, i) =>
        val from = seqNum * i / foldNum
        val until = seqNum * (i + 1) / foldNum
        val validateSeqs = seqs.slice(from, until)
        val emOutputs = emOutputSeq(i)
        val predict = new Predict(emOutputs, validateSeqs, predCond)
        s ::: predict.run.toList
      }
      Predict.writeSeqs(predicted, output)
    }
  }
}