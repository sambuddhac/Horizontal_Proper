# GameCapEx (Game-Theoretic Capacity Expansion Planning Model)

## Introduction

GameCapEx is a mathematical model and a software tool-suite (based on the afore-mentioned mathematical model) for strategic investment coordination among multiple regions for expanding both transmission (which is referred to as `Horizontal Coordination`) and generation (which is referred to as `Vertical Coordination`: currently, under development) capacities. It models the problem as a Mixed Integer Linear Programming (MILP) problem and solves both the optimal investment at building new transmission lines and generators, as well as, the operational problem of optimally dispatching the generators to meet demand, while obeying the power-flow network constraints. It does so, using a linearized or "DC-" OPF (Optimal Power Flow) approximation. The entire model for this capacity-expansion planning has been implemented in three modes, which can be selected by the user for the purposes of comparison or suitability of use. The different modes are: (1) Centralized, (2) Non-Cooperative Nash-Equilibrium (NE) based game-theoretic (to be implemented), and (3) Distributed Stochastic Optimization based algorithmic market mechanism design. GameCapEx is written in Python and Julia/JuMP, with ability to choose optimization solvers from a range of available options (Gurobi, MOSEK, GLPK etc.).

## Installation
In order to install and run the GameCapEx software, go through the following steps: 
1. Clone this repository to your local machine and navigate to the top level Horizontal_Proper folder.

```sh
git clone https://github.com/sambuddhac/Horizontal_Proper.git
```

2. If you are interested to run the Centralized model, go inside the folder, `CentralControl_PyJuMP` and create a conda environment named `centralHorCoord` by the `environment.yml` file, by typing

```sh
conda env create -f environment.yml
```

3. To activate `centralHorCoord` environment, please type

```sh
conda activate centralHorCoord
```

4. If you are interested to run either the NE based game-theoretic model or the market-mechanism desing model, go inside the `DistMechDes_PyJuMP` and repeat step 2. For activating the conda environment `distHorCoord`, this time, type

```sh
conda activate distHorCoord
```

5. For the first-time use, you will also need to create the Julia virtual enviroment (for which the first prerequisite is to have Julia installed on your machine), which can be done by navigating inside the `Julia_src` folder (within each of `CentralControl_PyJuMP` and `DistMechDes_PyJuMP` depending on which model you want to run) and by typing

```sh
julia julenv.jl
```

## Running code

### Settings



### Example notebooks



### Command line interface


## Contributing

