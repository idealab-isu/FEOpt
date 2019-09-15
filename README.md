# FEOpt

FEOpt is an optimization framework for the tuning of patient/subject-specific msterial parameters with finite element models. It is implemented in Python, Bash, and Matlab and uses [Continuity](https://continuity.ucsd.edu) for model evaluations.

Meshes of patient/subject heart scans are passively inflated to the measured pressure. Since the scans are of a loaded (inflated) ventricle, an iterative process is used to fetermine the unloaded geometry. The quality of differrent material parameter values is assessed via fitment with the Klotz curve.

## Features
This project contains a number of features that support the automation of model tuning: 
* Seperate submission scripts manage optimization and evaluation. Status and result updates are passed via worker-specific folders.
* The iterative process - consisting of a variable number of inflation and deflation steps - is handled via a Bash script. Stopping criteria are also implemented.
* Model parameters to be tuned are listed in a struct. The initialization of reference values and scaling is automatically handled.
* Models are evaluated in parallel during the optimization. Additionally, each model is evaluated on up to 8 cores; more cores don't provide significant runtime improvement.
* Results are logged and queried to prevent duplicate model evaluations.

## Dependencies
The core of the cardiovascular finite element model is [Continuity](https://continuity.ucsd.edu), which can be downloaded directly from their website.

This program also uses solvers from the [Matlab Optimization Toolbox](https://www.mathworks.com/products/optimization.html).

## Installation
1. Clone the FEOpt repository
2. Download and install Continuity
3. Set correct home and Continuity directories in Python, Bash, and Matlab files

## Using FEOpt
1. Submission of 'evaluate.script' and 'optimize.script' to a SLURM scheduler will tune the material parameters for the subject in that folder
2. Running 'clear' manually resets the folder for a new optimization run
3. 'PaperInflatePlot.m' can be run to produce the same plot style in the paper

## File Organization
* Each folder contains the files for one of four canine subjects
* '...cont6': Contain the geometric ventricle models and Continuity settings for each canine subject
* '...xls': Contains the node connectivity for a given subject
* 'optimize.script': SLURM job submission script for the optimizer
* 'InflationOpt...m': Matlab optimization manager and solver
* 'Inflate...vec.m': Vectorized model evaluation - sets flags and logs results of Continuity model
* 'evaluate.script': SLURM job submission script for the Continuity workers
* 'runCase.script': Bash script managing a specific Continuity evaluation
* 'Full_...py': Python script for inflation and deflation of the ventricle

## Publications
1. Joshua Mineroff, Andrew D. McCulloch, David Krummen, Baskar Ganapathysubramanian, and Adarsh Krishnamurthy (Accepted 2019). Optimization Framework for Patient-Specific Cardiac Modeling. _Cardiovascular Engineering and Technology_
2. Joshua Mineroff (2018). An Optimization and Uncertainty Quantification Framework for Patient-Specific Cardiac Modeling. _Graduate Theses and Dissertations_

## Licensing
This code is licensed under MIT license.
