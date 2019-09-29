package hmm_clustering

import java.nio.file.{Files, Paths}

import common.ConObs
import utils.Utils

/**
  * Created by tara on 2019/09/29.
  * Cluster genes with K-means
  */
class Clustering(seqs: Array[Array[Double]], clusteringCond: ClusteringCond, initialHMMs: Array[HMM]) {

  val m = seqs.length // 遺伝子の数

  // HMMとそれに対応する遺伝子のペアを返す
  def run(): Array[(HMM, Set[Int])] = {
    def runloop(hmmWithGenes: Array[(HMM, Set[Int])], loopNum: Int, allLogLikelihood: Double): Array[(HMM, Set[Int])] = {
      println("log likelihood in clustering: " + allLogLikelihood)
      // 各モデルの計算
      val newModels = hmmWithGenes.map(x => {
        val model = x._1
        val genes = x._2
        val geneSeqs = x._2.foldLeft(List(): List[Array[Double]]) {(s, x) => seqs(x) :: s}
        val bw = new Baum_Welch(model, geneSeqs)
        bw.run()
      })
      val newGeneSetsWithLL = assignMaxLikelihoodModel(newModels)
      val newHMMWithGenes = newModels.zip(newGeneSetsWithLL._1)
      val newLoglikelihood = newGeneSetsWithLL._2
      if (newLoglikelihood- allLogLikelihood < clusteringCond.delta(0) || loopNum > clusteringCond.loopMax(0)) {
        return newHMMWithGenes
      }
      runloop(newHMMWithGenes, loopNum + 1, newLoglikelihood)
    }
    val maxArraySetWithLogLikelihood = assignMaxLikelihoodModel(initialHMMs)
    runloop(initialHMMs.zip(maxArraySetWithLogLikelihood._1), 0, maxArraySetWithLogLikelihood._2)
  }

  // HMMから, それぞれに属するクラスタに属する遺伝子のindexのsetを返す 同時にその時の全てのlikelihoodのlogの積も返す
  // それぞれの対数尤度を計算して最も大きいものを割り当てる
  private def assignMaxLikelihoodModel(HMMs: Array[HMM]): (Array[Set[Int]], Double) = {
    val eachGeneIndex = seqs.map(geneSeq => {
      HMMs.zipWithIndex.foldLeft(-1.0E10, -1) {(s, x) =>
        val likelihood = x._1.calcLikelihood(geneSeq)
        if (likelihood > s._1) (likelihood, x._2) else s
      }
    })
    val allLogLikelihood = eachGeneIndex.foldLeft(0.0) {(s, x) => s + Math.log(x._1)}
    val maxArraySet = Array.fill[Set[Int]](m)(Set())
    eachGeneIndex.zipWithIndex.foreach(x => maxArraySet(x._1._2) = maxArraySet(x._1._2) + x._2)
    (maxArraySet, allLogLikelihood)
  }
}

object Clustering {
  def run(inputSeqs: String, output: String, condition: String): Unit = {
    val start = System.currentTimeMillis()
    val seqs = Clustering.readSeqsForClustering(inputSeqs)
    val clustCond = ClusteringCond.apply(condition)
    val initialHMMs = Clustering.readInitialHMMFile(clustCond.initialHMMFile)
    val clustering = new Clustering(seqs, clustCond, initialHMMs)
    val outputs = clustering.run()
    outputs.zipWithIndex.foreach(x => {
      println("gene in " + x._2 + "th cluster: " + x._1._2.toString())
    })
    val end = System.currentTimeMillis()
    print((end - start) + ":ms for clustering")
  }

  def readSeqsForClustering(input: String): Array[Array[Double]] = {
    val source = scala.io.Source.fromFile(input).getLines()
    source.map(s => {
      s.split("\t").map(_.toDouble)
    }).toArray
  }

  // 形式は一行ごとに各モデルの各stateの mu,sigma\t mu,sigma \t mu, sigma......
  def readInitialHMMFile(initialHMMFile: String): Array[HMM] = {
    val source = scala.io.Source.fromFile(initialHMMFile).getLines()
    source.map(s => {
      val mu_sig = s.split("\t").foldLeft(List(): List[Double], List(): List[Double]) {(s, x) => (x.split(",")(0).toDouble :: s._1, x.split(",")(1).toDouble :: s._2)}
      HMM(mu_sig._1.reverse.toArray, mu_sig._2.reverse.toArray, HMM.createInitTransitRL(mu_sig._1.length))
    }).toArray
  }
}