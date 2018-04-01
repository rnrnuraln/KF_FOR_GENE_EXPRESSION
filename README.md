KF_FOR_GENE_EXPRESSION
===


## Description

Kalman Filter is a widely used filtering method in the field of control engineering. It estimates hidden state of the system at a time point from time series of observed variable until the time point.

This tool aims to predict gene expression using Kalman filter, still being under construction.

There are some varieties of Kalman Filter model, and here we consider the model that consists of three main variables (hidden (x), observed (z), and environmental (u)), and two noise variables (v, w) normally distributed with mean 0. They are all real vector. Temporal behavior of each variable is described by linear stochastic difference equation shown below.

![equation for Kalman Filter](https://github.com/rnrnuraln/KF_FOR_GENE_EXPRESSION/blob/master/images/Screen%20Shot%202018-03-17%20at%2010.11.31.png)

![Explanation of Kalman Filter](https://github.com/rnrnuraln/KF_FOR_GENE_EXPRESSION/blob/master/images/kalman_picture.png)

Since we did not assume what exactly the hidden states are, we made a tool to estimate some parameters of Kalman Filter from data, so that observed variable fits the data. Thanks to the simplicity of the model, we can estimate the parameters using EM algortihm, just like as Baum-Welch algorithm used in parameter estimation of HMM. EM algorithm can find parameters that make log-likelihood high in a short time.

## Feature

Some original features are added, .

* To avoid over fitting, L1 or L2 regularization can be selected.

* Missing values of observed data are permitted.

* To make the model simple, covariane of Gaussian noise can be diagonal.

* Dimension of the hidden state can be determined by cross validation.

And more...


## Usage

Simple usage will be written here. See more details on docs/format.md.


* [-m mode -i input -o output -c condition]
 There are two mode in this code; mode and evaluation mode. We can change mode using -m option.
 Input sequence file (tsv file) is specified using -i option.
 Condition for EM algorithm is written in -c option.
 You can specify output directory using -o option.

** "learn" mode
 Estimates the parameters using EM algorithm.

** "evaluate" mode
 This mode is used to evaluate the performance of the tool.

** "predict" mode
 Using Kalman Filter, predict gene expression. This have not yet implemented.

## Install

To install the program, sbt is required.

https://www.scala-sbt.org/

Compile using sbt by assembly command.

Then, jar file for this program will be generated under target/scala-2.11.

## Tutorial

 Install sbt and enter the code below at the command line.

```
 sbt run -m learn -c tutorial/tutorial.cnd -i tutorial/observe_control_seq.tsv -o tutorial/result.kalf
```

 Then, results would be appeared under tutorial.

## Performance

Based on this algorithm, We developed a tool to estimate them and watched how well it estimates parameters from simulation data. The simulation data were some sequences of vector generated according to Kalman Filter model. We changed the number of the sequences in several ways and watched how relative error of estimated value changed with that. 

![log likelihood ratio for simulation data](https://github.com/rnrnuraln/KF_FOR_GENE_EXPRESSION/blob/master/images/log_likelihood.png)

Log-likelihood ratio (estimated/true) for the test data decreased with data increasing, which means that estimated parameters becomes near to true one. Hence, we can say that our program worked well. 


## Licence


Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, 
  this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, 
  this list of conditions and the following disclaimer in the documentation 
  and/or other materials provided with the distribution.
* The names of its contributors 
  may be used to endorse or promote products derived from this software 
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Author

These codes are written by Kentaro Tara, a graduated student in University of Tokyo.
