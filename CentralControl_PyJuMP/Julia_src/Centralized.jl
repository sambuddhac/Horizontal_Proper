using Pkg
Pkg.activate("GTheoryJulEnv")
#Defining sets

Z=network["zoneNumber"] #get the number of Zones from network Json files ###Should we aggregate these files for our centralized model??
###I couldn't find zone numbers in this file 
K=network["sharedCLines"]+network["internalCLines"] #The number of candidate Lines 
S=network["CountofScenarios"] #Scenarios
H=network["tranNumber"]+network["sharedELines"]    #Existing Lines
G=network["genNumber"] #The number of generators

### prob[s] should be defined, I can't find the related data
#Defining model
model = Model(GLPK.Optimizer)

#Defining variables
@variable(model, candLineDecision[1:K,1:Z], Bin) #Decision variable 
@variable(model, candFlowMW[1:S, 1:K])  #Power flowing on candidate candLineDecision 
@variable(model, 0 <= candPhaseAngleTo[1:S, 1:K]) #Phase angle decision for candidate line "To node"
@variable(model, 0 <= candPhaseAngleFrom[1:S, 1:K]) #Phase angle decision for candidate line "From node"
@variable(model, EFlowMW[1:S, 1:H])  #Power flowing on existing lines
@variable(model, 0 <= EPhaseAngleTo[1:S, 1:H]) #Phase angle decision for existing lines "To node"
@variable(model, 0 <= EPhaseAngleFrom[1:S, 1:H]) #Phase angle decision for existing  lines "From node"
@variable(model, 0 <= P_gen[1:S, 1:G]) #Power of generator
@variable(model, <= P_demand[1:S, 1:N])) #Power of demand
@variable(model, F) #Objective function

@constraint(model,sum(sum(prob[s].*sum(P_gen[s,:].*Gen["linCostCoeff"][:]).+sum(candLineDecision[:,z].*CandLine["interestRate"][:].*((CandLine["interestRate"][:]).^CandLine["lifeTime"][:])./((1+CandLine["interestRate"][:]).^CandLine["lifeTime"][:])-1))for s in S)for z in Z)) #Defining the OF
for z in 1:Z
    for s in 1:S
        for h in 1:H
            if Tran["tranZoneID"][h]==z
                @constraint(model, EFlowMW[s,h] .== EPhaseAngleFrom[s,h]./Tran["Reactance"][h] .- EPhaseAngleTo[s,h]./Tran["Reactance"][h]) #Constraint regarding the power flowing on existing lines
                @constraint(model, EFlowMW[s,h]<= Tran["lineLimit"][h])  # Line capacity constraints
                @constraint(model, -Tran["lineLimit"][h]<= EFlowMW[s,h])  # Line capacity constraints
            end
        end
        for k in 1:K
            if CandLine["tranZoneID"][h]==z
                @constraint(model, -10000 * (1-candLineDecision[k,z]) .<= candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./CandLine["Reactance"][k] .- candPhaseAngleTo[s,k]./CandLine["Reactance"][k]]) #Constraint regarding the power flowing on shared existing lines
                @constraint(model, candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./CandLine["Reactance"][k] .- candPhaseAngleTo[s,k]./CandLine["Reactance"][k]] .<= 10000 * (1-candLineDecision[k,z])) #Constraint regarding the power flowing on shared existing lines
                @constraint(model, candFlowMW[s,k]<= candLineDecision[k,z].*CandLine["lineLimit"][k])
                @constraint(model, -candLineDecision[k,z].*CandLine["lineLimit"][k]<=candFlowMW[s,k])
            end
        end
        for n in N
            if Node["nodeZoneID"][n]==z ####I think we should define a zoneID for each node
                @constraint(model, sum(P_gen[s,g] for g in G if Gen["genNodeID"][g]==n).- P_demand[s,n].== sum(candFlowMW[s,k] for k in K if CandLine["fromnode"][k]==n)-sum(candFlowMW[s,k] for k in K if CandLine["tonode"][k]==n)+sum(EFlowMW[s,h] for h in H if Tran["fromnode"][h]==n)-sum(EFlowMW[s,h] for h in H if Tran["tonode"][h]==n))
            end
        end
        for g in G
            if Node["genZoneID"][g]==z ###I think we should define a zoneID for each generator 
                @constraint(P_gen[s,g]<=Gen["PgMax"][g])
                @constraint(Gen["PgMin"][g]<=P_gen[s,g])
            end
        end
    end
end
@objective(model, Min, F)
status = optimize!(model)
@show(status, value.(candLineDecision))
