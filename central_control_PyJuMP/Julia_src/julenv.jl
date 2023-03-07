# import all of the packages (and specific versions) needed
import Pkg
using Pkg
#Pkg.activate("CentJulEnv")
Pkg.add("IJulia")

# core packages
Pkg.add(Pkg.PackageSpec(name="BenchmarkTools"))
Pkg.add(Pkg.PackageSpec(name="CSV")) 
Pkg.add(Pkg.PackageSpec(name="DataFrames"))
Pkg.add(Pkg.PackageSpec(name="DataStructures"))
Pkg.add(Pkg.PackageSpec(name="Dates"))
Pkg.add(Pkg.PackageSpec(name="LinearAlgebra"))
Pkg.add(Pkg.PackageSpec(name="MathOptInterface"))
Pkg.add(Pkg.PackageSpec(name="MathProgBase"))
Pkg.add(Pkg.PackageSpec(name="RecursiveArrayTools"))
Pkg.add(Pkg.PackageSpec(name="Statistics"))
Pkg.add(Pkg.PackageSpec(name="StatsBase"))
Pkg.add(Pkg.PackageSpec(name="YAML"))
Pkg.add(Pkg.PackageSpec(name="Plots"))
Pkg.add(Pkg.PackageSpec(name="XLSX"))
Pkg.add(Pkg.PackageSpec(name="GeometryBasics"))
# optimization packages
# Pkg.add(Pkg.PackageSpec(name="Cbc"))
# Pkg.add(Pkg.PackageSpec(name="Clp"))
Pkg.add(Pkg.PackageSpec(name="GLPK"))
# Pkg.add(Pkg.PackageSpec(name="Gurobi")) # uncomment this if you want to use Gurobi
Pkg.add(Pkg.PackageSpec(name="HiGHS"))
Pkg.add(Pkg.PackageSpec(name="Ipopt"))
Pkg.add(Pkg.PackageSpec(name="JuMP"))
Pkg.add(Pkg.PackageSpec(name="SCIP"))