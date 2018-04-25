KF_FOR_GENE_EXPRESSION
===


## Description

Kalman Filter is a widely used filtering method in the field of control engineering. It estimates hidden state of the system at a time point from time series of observed variable until the time point.

This tool aims to predict gene expression using Kalman filter, still being under construction.

There are some varieties of Kalman Filter model, and here we consider the model that consists of three main variables (hidden (x), observed (z), and environmental (u)), and two noise variables (v, w) normally distributed with mean 0. 
They are all real vector. Temporal behavior of each variable is described by linear stochastic difference equation shown below.

![equation for Kalman Filter](https://github.com/rnrnuraln/KF_FOR_GENE_EXPRESSION/blob/master/images/Screen%20Shot%202018-03-17%20at%2010.11.31.png)

![Explanation of Kalman Filter](https://github.com/rnrnuraln/KF_FOR_GENE_EXPRESSION/blob/master/images/kalman_picture.png)

Since we did not assume what exactly the hidden states are, we made a tool to estimate some parameters of Kalman Filter from data, so that observed variable fits the data. 
Thanks to the simplicity of the model, we can estimate the parameters using [EM algortihm](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm), just like as Baum-Welch algorithm used in parameter estimation of HMM. 
Using EM algorithm, parameters that increase log-likelihood of the model can be found efficiently.

## Estimation

 The parameters to estimate are described below.
 
* Matrix A, B, H, in the formula above.
* Covariance of noise variables of hidden and observed variables.
* Mean and covariances of the initial state of hidden variables.
 Since the model is described by difference equation, we have to specify initial state of the model in some way.

 Estimating the parameters of this model are written [here](http://mlg.eng.cam.ac.uk/zoubin/course04/tr-96-2.pdf).
 Based on the algorithm, I added some original features.

* To avoid overfitting, L1 or L2 regularization can be selected.
* Missing values of observed data are permitted.
* Since EM algorithm tends to convergence to a local optimum, initial values for EM algorithm are chosen randomly many times and the best parameters are taken.
* To make the model simple, covariance of Gaussian noise can be diagonal.
* Dimension of the hidden state can be determined by cross validation.
* ...

 We do not write detail of the algorithm here, because it is so complicated. 
 
## Prediction

 After estimating, we can predict gene expression level using Kalman Filtering method. 
 This mode has not yet implemented.

## Usage

Simple usage is written here. 

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
 Using Kalman Filter, predict gene expression. This has not yet implemented.
 
 In "learn" and "evaluate" mode, this tool has many parameters.
 We can arrange them using -c option.
 However, it is so complicated that we cannot describe in this space.
 The explanation of the condition would be written in docs/format.txt later.
 
## Install

To install the program, sbt is required.

https://www.scala-sbt.org/

After install, start sbt at this directory. It takes some time for initializing. Then, the interactive mode of sbt will be appeared. 

You can run the program in the interactive mode by typing "run".

You can also generate packaged jar file by typing "assembly" in the interactive mode.

## Tutorial

 Install and run sbt, and copy and paste the words below.

```
> run -m learn -c tutorial/tutorial.cnd -i tutorial/con_obs_seq.tsv -o tutorial/result.kalf
```

 "tutorial/con_obs_seq.tsv" shows the data to learn parameters.
 Dimension of observed variable(z in the formula above) is 5 and that of environmental variable(u in the formula above) is 3. 
 The number of independent sequence is 5 and the length of each sequence is 100.

 Condition for this program is written in "tutorial/tutorial.cnd". 
 The format of condition file is complicated and is hard to explain all of that.
 Let me pick up some feature of the condition.
 
* Matrix A, B, H are estimated using L1-regularization.
* Covariance matrix of the noise variable of the hidden layer is fixed to be identical.
* Covariance matrix of the noise variable of the observed layer is designed to be diagonal.

 It takes some time to complete compiling and running the program.
 After running, results would be appeared under tutorial directory.
 
 Results of the formats are generated at "tutorial/result.kalf".
 It shows results of estimation of each parameter and likelihood of the model.
 Explanation of the formats are written in "docs/format".

## Performance

Based on the estimation algorithm, We developed a tool to estimate them and watched how well it estimates parameters from simulation data. 
The simulation data were some sequences of vector generated according to Kalman Filter model. 
We changed the number of the sequences in several ways and watched how relative error of estimated value changed with that. 

![log likelihood ratio for simulation data](https://github.com/rnrnuraln/KF_FOR_GENE_EXPRESSION/blob/master/images/log_likelihood.png)

Log-likelihood ratio (estimated/true) for the test data decreased with data increasing, which means that estimated parameters becomes near to true one. 
Hence, we can say that our program worked well in the simulation data. 


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
