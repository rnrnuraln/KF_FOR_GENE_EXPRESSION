

/**
  * 引数として何を与えるかによってモードが変化する
  *
  */
object Main {

  case class Config(input: String = "", output: String = "", mode: String = "", condition: String = "", inputSeq: String = "")

  def main(args: Array[String]): Unit = {
    val parser = new scopt.OptionParser[Config]("kalman") {
      head("kalman", "1.0")
      opt[String]('i', "input").action((x, c) => c.copy(input = x)).text("specify the input directory")
      opt[String]('o', "output").required().action((x, c) => c.copy(output = x)).text("the output file")
      opt[String]('m', "mode").required().action((x, c) => c.copy(mode = x)).text("the mode")
      opt[String]('c', "condition").required().action((x, c) => c.copy(condition = x)).text("condition for algorithm, experiment, etc...")
      opt[String]('s', "inputSeq").action((x, c) => c.copy(inputSeq = x)).text("input of seqs for prediction")
    }

    parser.parse(args, Config()) match {
      case Some(Config(_, o, "evaluate", c, s)) => common.Evaluate.run(c, o)
      case Some(Config(i, o, "learn", c, s)) => common.Optimize.learn(i, o, c)
      case Some(Config(i, o, "predict", c, s)) => common.Predict.predict(i, s, c, o)
      case Some(Config(i, o, "read", _, _)) => common.EMOutputs.readLogLikelihoods(i)
      case Some(c) =>
      case None =>
      // arguments are bad, error message will have been displayed
    }
  }
}

