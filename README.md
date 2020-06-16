[![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-Ready--to--Code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/gpwhs/drifty-gr) 

# drifty-gr

University of Manchester: 1 semester long masters project

Objective: To parameterise the elliptical carousel model for pulsar emission through use of a bi-drifting pulse profile. 

Background: PSRSALA is an open source pulsar analysis package. If provided with geometric parameters, it is able to simulate pulsar pulse profiles using the elliptical carousel model. These profiles can be expressed as p3folds which give an average drift band over many pulsar periods. Simulated p3folds can be compared to experimental p3folds to gain insight into the properties of pulsars. 

Methodology: The elliptical carousel model takes 20+ individual parameters, so Monte Carlo was chosen to search the parameter space. Each pulsar generated was compared against the experimental data for goodness of fit. P3folds are images, and so the fit was defined via a difference image analysis method. The Monte Carlo method then optimised for the minimum pixel difference between the simulated and experimental images.  

Monte Carlo is compuationally inefficient, so in this project multithreading was utilised to parallelise the search. 
