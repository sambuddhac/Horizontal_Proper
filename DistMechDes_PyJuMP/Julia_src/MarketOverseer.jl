using Pkg
Pkg.activate("GTheoryJulEnv")
function marketOverseerOpt(lagrangeMultPi, K1lagrangeMultXi, K2lagrangeMultXi, H1lagrangeMultXi, H2lagrangeMultXi)
    size(zoneList)[1] #Zones
    K=size(candLine)[1] #Candidate shared Lines 
    S=size(scenVector)[1] #Scenarios
    H=size(SELine)[1]    #Shared existing Lines
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
    model = Model(GLPK.Optimizer)
    @variable(model, 0 <= candLineDecision[1:K] <= 1) #Decision variable 
    @variable(model, 0 <= candFlowMW[1:S, 1:K])  #Power flowing on candidate candLineDecision 
    @variable(model, 0 <= candPhaseAngleTo[1:S, 1:K]) #Phase angle decision for candidate line "To node"
    @variable(model, 0 <= candPhaseAngleFrom[1:S, 1:K]) #Phase angle decision for candidate line "From node"
    @variable(model, 0 <= SEFlowMW[1:S, 1:H])  #Power flowing on existing shared lines
    @variable(model, 0 <= SEPhaseAngleTo[1:S, 1:H]) #Phase angle decision for existing shared lines "To node"
    @variable(model, 0 <= SEPhaseAngleFrom[1:S, 1:H]) #Phase angle decision for existing shared lines "From node"
    @variable(model, F[1:S,1:Z])
    for z in 1:S
        for s in 1:Z
            @constraint(model, F[s,z] .==-[sum(lagrangeMultPi[z,:].*candLineDecision[:]) .+ sum(K2lagrangeMultXi[z,:,s].*candPhaseAngleTo[s,:]) .+ sum(K1lagrangeMultXi[z,:,s].*candPhaseAngleFrom[s,:]) .+ sum(H2lagrangeMultXi[z,:,s].*SEPhaseAngleTo[s,:]) .+ sum(H1lagrangeMultXi[z,:,s].*SEPhaseAngleFrom[s,:])])  #Objective Function
            for h in 1:H
                @constraint(model, SEFlowMW[s,h] .== SEPhaseAngleFrom[s,h]./SEReactance[h] .- SEPhaseAngleTo[s,h]./SEReactance[h]) #Constraint regarding the power flowing on shared existing lines
                @constraint(model, SEFlowMW[s,h]<= SECapacity[h])  # Capacity constraints
            end
            for k in 1:K
                if flag==1
                    @constraint(model, candLineDecision[k] in MOI.Integer())
                end
                @constraint(model, -10000 * (1-candLineDecision[k]) .<= candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./candReactance[k] .- candPhaseAngleTo[s,k]./candReactance[k]]) #Constraint regarding the power flowing on shared existing lines
                @constraint(model, candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./candReactance[k] .- candPhaseAngleTo[s,k]./candReactance[k]] .<= 10000 * (1-candLineDecision[k])) #Constraint regarding the power flowing on shared existing lines
                @constraint(model, candFlowMW[s,k]<= candCapacity[k])
            end
        end
    end
    @objective(model, Min, sum(F[:,:]))
    status = optimize!(model)
    @show(status, value.(candLineDecision))
    end
