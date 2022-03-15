This repository host the source codes for the project "Gene expression noise accelerates the evolution of a biological oscillator", authored by Yen Ting Lin (CCS-3, LANL) and Nicolas E. Buchler (North Carolina State University). 

The source codes were reviewed and approved with a LANL C-number C21109 by the Richard P. Feynman Center for Innovation (FCI) at the Los Alamos National Laboratory. 

# General structure of the codes

We developed two types of simulators for generating ***deterministic*** and ***stochastic*** trajectories of the gene dynamics. For the ***deterministic*** dynamics, we used ```scipy.integrate.solve_ivp``` to solve the differential equations. For ***stochastic dynamics***, we used an in-house Gillespie simulation implemented in c++ to enhance the computational efficiency. The model structures are hard-coded in the provided c++ source code (```/Ccode/repressilator.cpp``` and ```/Ccode/titration-oscilattor.cpp```). The users will have to compile the c++ code, and the executable files will be called by Python scripts. 

The implementation of the evolutionary process is implemented in Python. The process entails repetative (1) initiation of a population of cells with different biophysical parameters, (2) calling and executing a swarm of simulators, either ***deterministic*** or ***stochastic***, (3) computation of the fitness of each cell, (4) selecting the fitter ones for mutation, and (5) initiation of the next-generation of the cells.

<p align="center">
<img src='https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/figures/schematics.png' width='60%'>
</p>
