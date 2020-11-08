#= **These packages won't be needed here as I have already created a julia virtual environment in master, which we will be activating here. See below
using JuMP
using Xpress
using Dates
=#
using Pkg
Pkg.activate("GTheoryJulEnv")
function marketOverseerOpt(LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum)
    Z=size(zoneList)[1] #Zones
    K=size(candLine)[1] #Candidate shared Lines 
    S=size(scenVector)[1] #Scenarios
    H=size(SELine)[1]    #Shared existing Lines
    KN1=size(candFromRank)[1] #Candidate line "from nodes"
    KN2=size(candToRank)[1]  #Candidate line "to nodes"
    HN1=size(SEFromRank)[1] #Shared existing line "from nodes"
    HN2=size(SEToRank)[1]  #Shared existing line "to nodes"
    N=KN1+KN2+HN1+HN2   #nodes

    flag[k]=1         #Please define a flag (0 for finding lower bound and 1 for finding optimal value)
    #These parameters should be obtained 
    lagrangeMultPi[1:Z, 1:K] #Binary decision variable lagrange multiplier
    K1lagrangeMultXi[1:Z, 1:KN1, 1:S] #Lagrange multiplier associated with "From node" candidate lines
    K2lagrangeMultXi[1:Z, 1:KN2, 1:S] #Lagrange multiplier associated with "To node" candidate lines
    H1lagrangeMultXi[1:Z, 1:HN1, 1:S] #Lagrange multiplier associated with "From node" shared existing lines
    H2lagrangeMultXi[1:Z, 1:HN1, 1:S] #Lagrange multiplier associated with "To node" shared existing lines
    SEReactance[1:H] #Reactance of shared existing lines
    candReactance[1:K] #Reactance of candidate lines

    model = Model(GLPK.Optimizer)


    @variables model begin
        [flag[k]==1], candLineDecision[1:K], Bin #Binary decision variable indicating capacity expansion
        [flag[k]==0], 0 <= candLineDecision[1:K] <= 1 #Decision variable in relaxed optimization
        0 <= candFlowMW[1:S, 1:K]  #Power flowing on candidate candLineDecision 
        0 <= canPhaseAngleTo[1:S, 1:KN2] #Phase angle decision for candidate line "To node"
        0 <= canPhaseAngleFrom[1:S, 1:KN1] #Phase angle decision for candidate line "From node"
        0 <= SEFlowMW[1:S, 1:H]  #Power flowing on existing shared lines
        0 <= SEPhaseAngleTo[1:S, 1:HN2] #Phase angle decision for existing shared lines "To node"
        0 <= SEPhaseAngleFrom[1:S, 1:HN1] #Phase angle decision for existing shared lines "From node"

    end

    @constraints model begin

        [z=1:Z, s=1:S], F0[s,z] ==-[sum(lagrangeMultPi[z,:]*candLineDecision[:]) +sum(K2lagrangeMultXi[z,:,s]*canPhaseAngleTo[s,:])+
        sum(K1lagrangeMultXi[z,:,s]*canPhaseAngleFrom[s,:]) + sum(H2lagrangeMultXi[z,:,s]*SEPhaseAngleTo[s,:])+
        sum(H1lagrangeMultXi[z,:,s]*SEPhaseAngleFrom[s,:])]  #Objective Function
        [s=1:S, h=1:H], SEFlowMW[s,h] == SEPhaseAngleFrom[s,h]/SEReactance[h] - SEPhaseAngleTo[s,h]/SEReactance[h] #Constraint regarding the power flowing on shared existing lines
        [s=1:S, k=1:K], -10000 * (1-candLineDecision[k]) <= candlowMW[s,k] - [candPhaseAngleFrom[s,k]/candReactance[k] - candPhaseAngleTo[s,k]/candReactance[k]]<=10000 * (1-candLineDecision[k]) #Constraint regarding the power flowing on shared existing lines
    
        [h = 1:H], SEFlowMW[:, h]<= SECapacity[h]  # Capacity constraints
        [k= 1:K], candFlowMW[:,k]<= candCapacity[k]
    end

    @objective(model, Min, sum(F0[:,:]))

    optimize!(model)
end
