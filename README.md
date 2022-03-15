# SSTMetric-Manuscript-Code
This is the repository archiving code for the manuscript -- Self-Supervised Metric Learning in Multi-View Data: A Downstream Task Perspective

## File Introduction

The code provides all functions and scripts necessary to reproduce the results of the numerical experiments in the manuscript. In particular, the R scripts are divided into two collections: Simulation and RealData. In Simulation, the code includes four R scripts to reproduce the results in Section 6.1: 

* sampleidentification.R for Table 3 
* twosample.R for Figure 3
* kmeans.R for Table 4
* knn.R for Figure 4

In RealData, the code includes two R scripts to reproduce the results in Section 6.2: 

* MNIST.R for the left part of Table 5
* Fashion_MNIST.R for the right part of Table 5

In these R scripts, data are generated or accessed by R functions. All results can be reproduced if the random seeds are chosen as in the R scripts.
