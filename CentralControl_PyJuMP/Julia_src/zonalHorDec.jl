import Pkg
Pkg.add("Gurobi")
Pkg.add("GLPK")
Pkg.add("MathOptInterfaceMosek")
Pkg.add("MathOptInterface")
Pkg.add("Cbc")
Pkg.add("Ipopt")
using JuMP
using Gurobi
using GLPK
using MathOptInterfaceMosek
using Cbc
using Ipopt
using MathOptInterface

function zonalHorDec(
   )
   model = Model(with_optimizer(Gurobi.Optimizer))
   model = Model(with_optimizer(GLPK.Optimizer))
   model = Model(with_optimizer(MathOptInterfaceMosek.Optimizer))
   model = Model(with_optimizer(Cbc.Optimizer))
   model = Model(with_optimizer(Cbc.Optimizer))

   @variables model begin
   end 

   @constraints model begin
   end

   @objective()
   end

   results = Dict()

   return results
end
