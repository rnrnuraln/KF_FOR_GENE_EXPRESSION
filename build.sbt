lazy val root = (project in file(".")).
  settings(
    name := "KF_FOR_GENE_EXPRESSION",
    version := "1.0",
    scalaVersion := "2.11.8"
  )

libraryDependencies  ++= Seq(
  "org.scalanlp" %% "breeze" % "0.13.1",
  "org.scalanlp" %% "breeze-natives" % "0.13.1",
  "org.scalanlp" %% "breeze-viz" % "0.13.1"
)
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1"
libraryDependencies += "com.github.scopt" %% "scopt" % "3.5.0"

libraryDependencies += "org.scalactic" %% "scalactic" % "3.0.1"
libraryDependencies += "org.scalautils" % "scalautils_2.11" % "2.1.5"
libraryDependencies += "org.scalatest" %% "scalatest" % "3.0.1" % "test"
resolvers ++= Seq(
  "Sonatype Snapshots" at "https://oss.sonatype.org/content/repositories/releases/"
)


libraryDependencies +=
  "com.typesafe.akka" %% "akka-actor" % "2.5.0"

