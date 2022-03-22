This repository host the source codes for the project "Gene expression noise accelerates the evolution of a biological oscillator", authored by Yen Ting Lin (CCS-3, LANL) and Nicolas E. Buchler (North Carolina State University). 

The source codes were reviewed and approved with a LANL C-number C21109 by the Richard P. Feynman Center for Innovation (FCI) at the Los Alamos National Laboratory. 

# General structure of the codes

We developed two types of simulators for generating ***deterministic*** and ***stochastic*** trajectories of the gene dynamics. For the ***deterministic*** dynamics, we used ```scipy.integrate.solve_ivp``` to solve the differential equations. For ***stochastic dynamics***, we used an in-house Gillespie simulation implemented in c++ to enhance the computational efficiency. The model structures are hard-coded in the provided c++ source code ([```/Ccode/repressilator.cpp```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/CTMC_simulations/repressilator.cpp) and [```/Ccode/titration-oscilattor.cpp```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/CTMC_simulations/titrationOscillator.cpp)). The users will have to compile the c++ code, and the executable files will be called by Python scripts. 

The implementation of the evolutionary process is implemented in Python. The process entails repetative (1) initiation of a population of cells with different biophysical parameters, (2) calling and executing a swarm of simulators, either ***deterministic*** or ***stochastic***, (3) computation of the fitness of each cell, (4) selecting the fitter ones for mutation, and (5) initiation of the next-generation of the cells.

<p align="center">
<img src='https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/figures/schematics.png' width='60%'>
</p>

All the evolutionary processes were coded in the Python language and we provided the Jupyter notebooks, located in [```/Evolutionary_Processes```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/tree/main/Evolutionary_Processes). For those evolutionary processes whose underlying simulator is stochastic, the users need to compile the c++ source code ([```/Ccode/repressilator.cpp```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/CTMC_simulations/repressilator.cpp) and [```/Ccode/titration-oscilattor.cpp```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/CTMC_simulations/titrationOscillator.cpp)) to executables ([```/Evolutionary_processes/r_evolution.out```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/r_evolution.out) and [```/Evolutionary_processes/t_evolution.out```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/t_evolution.out)) in the same folder where the evolutionary codes are. A folder named [```evoBuffer```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/evoBuffer/) for storing temporary trajectories needs to be created. 

A set of evolutionary processes include:

#### Repressilator:

[```/Evolutionary_processes/repressilator-deterministicEvolution-r0.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/repressilator-deterministicEvolution-r0.ipynb): 1D evolutionary process of the biophysical parameter $r_0$ in the repressilator model with the deterministic gene dynamics; Fig. 3 in the preprint 
[```/Evolutionary_processes/repressilator-deterministicEvolution-r1.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/repressilator-deterministicEvolution-r1.ipynb): 1D evolutionary process of the biophysical parameter $r_1$ in the repressilator with the deterministic gene dynamics; Fig. 3 in the preprint 
[```/Evolutionary_processes/repressilator-stochasticEvolution-r0.ipynb.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/repressilator-stochasticEvolution-r0.ipynb): 1D evolutionary process of the biophysical parameter $r_0$ in the repressilator with the stochastic gene dynamics; Fig. 3 in the preprint 
[```/Evolutionary_processes/repressilator-stochasticEvolution-r1.ipynb.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/repressilator-stochasticEvolution-r1.ipynb): 1D evolutionary process of the biophysical parameter $r_1$ in the repressilator with the stochastic gene dynamics; Fig. 3 in the preprint 

#### Titration oscillator:

[```/Evolutionary_processes/titrationOscillator-deterministicEvolution-betaBX.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-deterministicEvolution-betaBX.ipynb): 1D evolutionary process of the biophysical parameter $\beta^B_X$ in the titration model with the deterministic gene dynamics; Fig. 4 in the preprint 

[```/Evolutionary_processes/titrationOscillator-deterministicEvolution-betaFX.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-deterministicEvolution-betaFX.ipynb): 1D evolutionary process of the biophysical parameter $\beta^F_X$ in the titration model with the deterministic gene dynamics; Fig. 4 in the preprint 

[```/Evolutionary_processes/ttitrationOscillator-stochasticEvolution-betaBX.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-stochasticEvolution-betaBX.ipynb): 1D evolutionary process of the biophysical parameter $\beta^B_X$ in the titration model with the stochastic gene dynamics; Fig. 4 in the preprint 

[```/Evolutionary_processes/titrationOscillator-stochasticEvolution-betaFX.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-stochasticEvolution-betaFX.ipynb): 1D evolutionary process of the biophysical parameter $\beta^F_X$ in the titration model with the stochastic gene dynamics; Fig. 4 in the preprint 

[```/Evolutionary_processes/titrationOscillator-2D-deterministic-landscape.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-2D-deterministic-landscape.ipynb): Grid-based computation of the 2D fitness landscape of the titration model with the deterministic gene dynamics. The output is saved to [```detLandscape.npz```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/detLandscape.npz) for the visualization in the following 2D/3D evolutionary codes. 

[```/Evolutionary_processes/titrationOscillator-2D-stochastic-landscape.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-2D-stochastic-landscape.ipynb): Grid-based computation of the 2D fitness landscape of the titration model with the stochastic gene dynamics. The output is saved to [```stoLandscape.npz```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/stoLandscape.npz) for the visualization in the following 2D/3D evolutionary codes. 

[```/Evolutionary_processes/titrationOscillator-deterministicEvolution-betaBX-betaFX.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-deterministicEvolution-betaBX-betaFX.ipynb): 2D evolutionary process of the biophysical parameters $\beta^B_X$ and $\beta^F_X$ in the titration model with the deterministic gene dynamics; Fig. 5 in the preprint 

[```/Evolutionary_processes/titrationOscillator-stochasticEvolution-betaBX-betaFX.ipynb.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-stochasticEvolution-betaBX-betaFX.ipynb): 2D evolutionary process of the biophysical parameters $\beta^B_X$ and $\beta^F_X$ in the titration model with the stochastic gene dynamics; Fig. 5 in the preprint 

[```/Evolutionary_processes/titrationOscillator-stochasticEvolution-betaBX-betaFX-Omega.ipynb.ipynb```](https://github.com/lanl/In-Silico-Evolution-of-oscillatory-gene-dynamics/blob/main/Evolutionary_Processes/titrationOscillator-stochasticEvolution-betaBX-betaFX-Omega.ipynb): 3D evolutionary process of the biophysical parameters $\beta^B_X$, $\beta^F_X$, and $\Omega$ in the titration model; Fig. 6 in the preprint 


