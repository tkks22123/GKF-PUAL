# Implementation of group kernel-free PU classifier with asymmetric loss (GKF-PUAL)
This is the code for GKF-PUAL in the paper "GKF-PUAL: A group kernel-free approach to positive-unlabelled learning with variable selection" [1].

* ```EX_fun_PU.R``` is a library of the functions for model optimisation via Alternating Direction Method of Multipliers (ADMM) [2].
* ```XXX_RRGLLC_0.5_1.R``` is the original script to train GKF-PUAL on the dataset XXX from UCI Machine Learning Repository with irrelevant features embedded. How the datasets were preprocessed can be find in Section 4.1 of [1].
* ```GKF_GLPUAL_pen_0.25.R``` is a pure demo on dataset Pen-Based Recognition of Handwritten Digits (Pen) [3] with structure optimised and detailed notes added.

## Requirements
* R >= 4.1.0
* mvtnorm >= 1.1.3
* expm >= 0.999.6


## Preparation
Before you run the scripts to train GKF-PUAL on PEN to assess the performance of GKF-PUAL, there are some preparation to do:
* Download the datasets form UCI Machine Learning Repository: https://archive.ics.uci.edu/.
* Substitute the path of ```EX_fun_PU.R``` in the script with the path of ```EX_fun_PU.R``` in your device. 
* Substitute the path of the dataset in the script with the path of the corresponding dataset in your device. 


## Quick start
Run ```XXX_RRGLLC_0.5_1.R``` or ```GKF_GLPUAL_pen_0.25.R``` and wait for a while. Then the F1-score and error rate for classification on the test set will be output on R console. The hyper-parameters in ```GKF_GLPUAL_pen_0.25.R``` has been determined previously via the method in Section 4.3 of paper [1]. You can also set the value of them manually in the block "The setting of hyper-parameters" from row 27 to row 40.


## Contact
For more details to discuss, please email me at xiaoke.wang.18@alumni.ucl.ac.uk.

## Reference

[1] Wang X, Zhu R, Xue J H. GKF-PUAL: A group kernel-free approach to positive-unlabeled learning with variable selection[J]. Information Sciences, 2025, 690: 121574.

[2] Boyd S, Parikh N, Chu E, et al. Distributed optimization and statistical learning via the alternating direction method of multipliers[J]. Foundations and Trends® in Machine learning, 2011, 3(1): 1-122.

[3] Alpaydin, E. & Alimoglu, F. (1996). Pen-Based Recognition of Handwritten Digits [Dataset]. UCI Machine Learning Repository. https://doi.org/10.24432/C5MG6K.
