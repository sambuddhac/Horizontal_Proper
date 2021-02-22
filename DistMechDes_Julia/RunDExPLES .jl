# current working directory
working_path = pwd()

# DExPLES path
dexples_path = pwd()

#Settings file path
settings_path = joinpath(pwd(), "DExPLES_settings.yml")

# Load GenX modules
push!(LOAD_PATH, dexples_path)
println(settings_path)
println("Loading packages")

using HorMILPDist
using YAML
using Dates
using DataFrames
using Gurobi
using CPLEX

println(now())

# Load inputs
push!(LOAD_PATH, working_path)

inpath="$working_path/Example_1Zone_REGas_test_short"
setup = YAML.load(open(settings_path))

inputs=Dict()

# KickStart the model
println("Loading inputs and starting the model")
inputs = HorMILPCentral(setup,inpath)