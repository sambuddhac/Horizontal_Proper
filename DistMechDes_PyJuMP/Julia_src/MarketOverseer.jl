using Pkg
Pkg.activate("GTheoryJulEnv")
function milp_marketOverseer(network::Dict, sharedCLines::Dict,internalCLines::Dict, , sharedELines::Dict, Tran::Dict, lagrangeMultPi::Dict, lagrangeMultXi::Dict, setup::Dict)
    K=network["sharedCLines"] #The number of shared candidate lines
    S=network["CountofScenarios"] #Scenarios
    H=network["sharedELines"]#Shared existing lines
    T=network["Hours"] #Total number of hours considered (whole year)
    N=network["nodeNumber"]
    zoneList = Array{Int8}(undef, Z) #Array of zones ####I think we should define that and add zone ID for each network
    Z=size(zoneList) #Total number of zones 
    flag=1 #If flag is 1 variable candLineDecision is binary, otherwise it is not (Determining relaxed LP or MIP)
"""These parameters should be obtained or defined
lagrangeMultPi=[z,k] #Binary decision variable lagrange multiplier
K1lagrangeMultXi=[z,k,s] #Lagrange multiplier associated with "From node" candidate lines
K2lagrangeMultXi=[z,k,s] #Lagrange multiplier associated with "To node" candidate lines
H1lagrangeMultXi=[z,h,s] #Lagrange multiplier associated with "From node" shared existing lines
H2lagrangeMultXi=[z,h,s] #Lagrange multiplier associated with "To node" shared existing lines
SEReactance=[h] #Reactance of shared existing lines
candReactance=[k] #Reactance of candidate lines
SECapacity=[h]
candCapacity=[k]
"""
    MOMod = Model(GLPK.Optimizer)
    @variable(MOMod, 0 <= candLineDecision[1:K] <= 1) #Decision variable 
    @variable(MOMod, 0 <= candFlowMW[1:S, 1:K])  #Power flowing on candidate candLineDecision 
    @variable(MOMod, 0 <= candPhaseAngle[1:S, 1:K]) #Phase angle decision for candidate line 
    @variable(MOMod, 0 <= SEFlowMW[1:S, 1:H])  #Power flowing on existing shared lines
    @variable(MOMod, 0 <= SEPhaseAngle[1:S, 1:H]) #Phase angle decision for existing shared lines 
    @variable(MOMod, F[1:S,1:Z])
    
    for z in 1:S
        for s in 1:Z
            @expression(MOMod, F[1:S,1:Z], -[sum(lagrangeMultPi[z,:].*candLineDecision[:]) 
                    .+ sum(lagrangeMultXi[z,:,s].*candPhaseAngle[s,:])
                    .+ sum(lagrangeMultXi[z,:,s].*SEPhaseAngle[s,:])])
    
            
            for h in 1:H
                @constraint(model, SEFlowMW[s,h] .== SEPhaseAngleFrom[s,h]./SEReactance[h] .- SEPhaseAngleTo[s,h]./SEReactance[h]) #Constraint regarding the power flowing on shared existing lines
                @constraint(model, SEFlowMW[s,h]<= SECapacity[h])  # Capacity constraints
            end
            for k in 1:K
                if flag==1
                    @constraint(model, candLineDecision[k] in MOI.Integer())
                end
                @constraint(model, -10000 * (1-candLineDecision[k]) .<= candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./candReactance[k] .- candPhaseAngleTo[s,k]./candReactance[k]]) #Constraint regarding the power flowing on shared candidate lines
                @constraint(model, candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./candReactance[k] .- candPhaseAngleTo[s,k]./candReactance[k]] .<= 10000 * (1-candLineDecision[k])) #Constraint regarding the power flowing on shared candidate lines
                @constraint(model, candFlowMW[s,k]<= candCapacity[k])
            end
        end
    end
    @objective(model, Min, sum(F[:,:]))
    status = optimize!(model)
    @show(status, value.(candLineDecision))
    end
