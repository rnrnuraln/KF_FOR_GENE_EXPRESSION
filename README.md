KF_FOR_GENE_EXPRESSION
===

Overview

Gene prediction
KF_FOR_BIOLOGY provides a good prediction of temporal phenomena on gene expression, using Kalman Filter.

## Description

Kalman Filter is a method that predicts the system using HMM.
See here for more details here(the paper or URL or somewhat for Kalman Filter).

These codes are written by Kentaro Tara, a university student in University of Tokyo.

## Usage

Simple usage will be written here. See more details on docs/format.md

 Create params and control-seq randomly.
 n: the dimension of x, hidden variable
 m: the dimension of u, control variable
 l: the dimension of z, measured variable
 T: the length of control sequence.
 All args must be integer.
 Outputs must be specified by

* learn [-u dir]
 Run the Kalman filter, according to parameters and control-seq.
 You must specify the directory that includes their data.
 In the directory, you must prepare directories shown below.
 A
 B
 H
 R
 Q
 control_seq
 All of them must be CSV.

* learn [-u dir]
 Estimates the parameters for .
 We use Baum-Welch-like algorithm to
 We described that

## Install

Command line

## Result

I am great! This program worked well.

## Author

Kentaro Tara
