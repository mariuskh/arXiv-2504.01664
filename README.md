# arXiv-2504.01664

This repository contains code to reproduce the results given in https://arxiv.org/abs/2504.01664. 

The code is written in Julia and is primarily using the QuantumOptics.jl package. 

The code is organized as follows:

- The src folder contains all the code needed to run the simulations.

- The simulation folder contains scripts for running example simulations. The results are written to files in the jld2 format.

- The data folder contains data corresponding to each simulation, i.e. the output of the simulation scripts have been moved to the data folder.

- The analysis folder contains notebooks where the simulation data is analyzed and plots similar to the figures in the preprint can be produced. 


The code has been tested with the following packages installed: 

```
Pkg.installed()
> 
Dict{String, VersionNumber} with 13 entries:
  "CSV"                     => v"0.10.15"
  "HypergeometricFunctions" => v"0.3.28"
  "LaTeXStrings"            => v"1.4.0"
  "IJulia"                  => v"1.26.0"
  "Plots"                   => v"1.40.9"
  "PyPlot"                  => v"2.11.5"
  "FileIO"                  => v"1.16.6"
  "Dates"                   => v"1.11.0"
  "JLD2"                    => v"0.5.11"
  "DifferentialEquations"   => v"7.16.1"
  "SpecialFunctions"        => v"2.5.0"
  "QuantumOptics"           => v"1.2.1"
  "PlotlyJS"                => v"0.18.15"
```
