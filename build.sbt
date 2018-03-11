import AssemblyKeys._ // put this at the top of the file
lazy val root = (project in file(".")).
  settings(
    name := "Kalman",
    version := "1.0",
    scalaVersion := "2.11.8"
  )

libraryDependencies  ++= Seq(
  // Last stable release
  "org.scalanlp" %% "breeze" % "0.13.1",

  // Native libraries are not included by default. add this if you want them (as of 0.7)
  // Native libraries greatly improve performance, but increase jar sizes.
  // It also packages various blas implementations, which have licenses that may or may not
  // be compatible with the Apache License. No GPL code, as best I know.
  "org.scalanlp" %% "breeze-natives" % "0.13.1",

  // The visualization library is distributed separately as well.
  // It depends on LGPL code
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

assemblySettings