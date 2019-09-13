# FEOpt

FEOpt is an optimization framework for the tuning of patient/subject-specific msterial parameters with finite element models. It is implemented in Python, Bash, and Matlab and uses [Continuity](https://continuity.ucsd.edu) for model evaluations.

Meshes of patient/subject heart scans are passively inflated to the measured pressure. Since the scans are of a loaded (inflated) ventricle, an iterative process is used to fetermine the unloaded geometry. The quality of differrent material parameter values is assessed via fitment with the Klotz curve.

## Features
This project contains a number of features that support the automation of model tuning: 
* The iterative process - consiting of a variable number of inflation and deflation steps - is handled via a Bash script. Stopping criteria are also implemented.
* Model parameters to be tuned are listed in a struct. The initialization of reference values and scaling is automatically handled.
* Models are evaluated in parallel during the optimization. Additionally, each model is evaluated on up to 8 cores; more cores don't provide significant runtime improvement.
* Results are logged and queried to prevent duplicate model evaluations.

## Dependencies
The core of the cardiovascular finite element model is [Continuity](https://continuity.ucsd.edu), which can be downloaded directly from their website.

This program also uses solvers from the [Matlab Optimization Toolbox](https://www.mathworks.com/products/optimization.html).

## Installation
1. Clone the FEOpt repository
2. Download and install Continuity

## Using CircOpt
1. Follow the example syntax in Run.m to tune a model
2. Change the '' and 'data...m' files (and make sure to reference that file in Run.m) to provide a patient-specific mesh and data for tuning

## File Organization
* TBD

## Publications
1. Joshua Mineroff, Andrew D. McCulloch, David Krummen, Baskar Ganapathysubramanian, and Adarsh Krishnamurthy (Accepted 2019). Optimization Framework for Patient-Specific Cardiac Modeling. _Cardiovascular Engineering and Technology_
2. Joshua Mineroff (2018). An Optimization and Uncertainty Quantification Framework for Patient-Specific Cardiac Modeling. _Graduate Theses and Dissertations_

## Licensing
This code is licensed under MIT license.
