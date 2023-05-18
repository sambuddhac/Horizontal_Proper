# GameCapEx (Game-Theoretic Capacity Expansion Planning Model)

## Introduction

GameCapEx is a mathematical model and a software tool-suite (based on the afore-mentioned mathematical models) for strategic investment coordination among multiple regions for expanding both transmission (which is referred to as `Horizontal Coordination`) and generation (which is referred to as `Vertical Coordination`: currently, under development) capacities. It models the problem as a Mixed Integer Linear Programming (MILP) problem and solves both the optimal investment at building new transmission lines and generators, as well as, the operational problem of optimally dispatching the generators to meet demand, while obeying the power-flow network constraints. It does so, using a linearized or "DC-" OPF (Optimal Power Flow) approximation. The entire model for this capacity-expansion planning has been implemented in three modes, which can be selected by the user for the purposes of comparison or suitability of use. The different modes are: (1) Centralized, (2)Cooperative, and (3) Distributed Stochastic Optimization based algorithmic market mechanism design which includes "Stage1" and "Stage2".  All of them are written in Julia using JuMP, with ability to choose optimization solvers from a range of available options (Gurobi, MOSEK, GLPK etc.).

## Running code

In order to run the code, you need to have Julia and Conda installed in your computer. Use [this link](https://julialang.org/downloads/) to download and install Julia. Conda full package can be also downloaded and installed by following the instructions in [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html) or [this website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html). Then, you need to open the Example folder. Select the file you want to run. Open the file and run the codes in the Jupyter Notebook.

### Settings



### Example notebooks



### Command line interface


## Contributing

