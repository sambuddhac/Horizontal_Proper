using Pkg
Pkg.activate("GTheoryJulEnv")
Z=size(zoneList)[1] #Zones
S=size(scenVector)[1] #Scenarios
H=size(SELine)[1]    #Shared Lines
N=size(NodeList)[1] #It should be equal to H*2
Iteration=100 #Maximum number of iteration
FF=0


"""These parameters should be obtained or defined
epsilon
Eta
Gamma
MarginalCost[z1,n] #Marginal cost of genrator g located in zone z
Prob[s] #Pobability of scenario s
lagrangeMulti[it=0,z1,z2,h] #initial value of Second-stage Lagrange Multiplier
lagrangeMultj[it=0,z1,z2,h] #initial value of Second-stage Lagrange Multiplier
SReactance[z1,h] #Reactance of shared line h between z1 and z2
SCapacity[z1,h] #Capacity of shared line h between z1 and z2
genCapacity[z1,n] #Generating limits of generator at node n and zone z1
P_demand[z1,n] #Demand of node n located at zone z1
SL[z1,z2,h] #A matrix with 1 and 0 elements, 1 if h is a shared line between z1 and z2, otherwise 0
"""

model = Model(GLPK.Optimizer)
@variable(model,0 <= ThetaSLj[1:S,1:Iteration,1:Z,1:Z,1:H]) #Belief of each zone about the voltage phase angle of the node belonging to other Zones
@variable(model,0 <= ThetaSLi[1:S,1:Iteration,1:Z,1:H]) #Belief of each zone about the voltage phase angle of the node belonging to itself
@variable(model,0 <= ThetaToNode[1:S,1:N,1:Z,1:H]) #Phase angle decision for existing shared lines "To node"
@variable(model,0 <= ThetaFromNode[1:S,1:N,1:Z,1:H]) #Phase angle decision for existing shared lines "From node"
@variables(model,begin
0 <= P_gen[1:S,1:Z,1:N] #Generator capacity at each node
0 <= flow[1:S,1:Z,1:H] #Flow of each line 
0 <= OF[1:Z,1:Iteration] #Objective Function
end)


for it in 1:Iteration
    for z1 in 1:Z
        for s in 1:S
            @constraint(model, sum(P_gen[s,z1,:]).- sum(P_demand[s,z1,:]).== sum(flow[s,z1,:])) #(3b) The balance-related constraint
            for h in 1:H
                if SL[z1,:,h]==0
                    @constraint(model, flow[s,z1,h].==(ThetaFromNode[s,z1,h].-ThetaToNode[s,z1,h])./SReactance[z1,h]) #(3c) Flow of non-shared lines
                else
                    for z2 in 1:Z
                        @constraint(model, flow[s,z1,h].==(ThetaSLi[s,it,z1,h].-ThetaSLj[s,it,z1,z2,h])./SReactance[z1,h]) #(3d) Flow of shared lines
                    end
                    @constraint(model,OF[z1,it] .== sum(Prob[s].*sum(sum(sum((MarginalCost[z1,n].*P_gen[s,z1,n].+Eta*(ThetaSLi[s,it,z1,h].*ThetaSLi[s,it-1,z1,h].-ThetaSLi[s,it,z1,h].*ThetaSLj[s,it-1,z2,z1,h].+ThetaSLj[s,it,z1,z2,h].*ThetaSLj[s,it-1,z1,z2,h].+ThetaSLj[s,it,z1,z2,h].+ThetaSLi[s,it-1,z2,h]).+(Gamma/2)*((ThetaSLi[s,it,z1,h].-ThetaSLi[s,it-1,z1,h])^2+(ThetaSLj[s,it,z1,z2,h].-ThetaSLj[s,it-1,z1,z2,h])^2).+lagrangeMulti[it,z1,z2,h].*ThetaSLi[s,it,z1,h].+lagrangeMultj[it,z1,z2,h].*ThetaSLj[s,it,z1,z2,h]) for n in N) for h in H) for z2 in 1:Z) for s in 1:S)) #Objective Function
                end
                @constraint(model,-SCapacity[z1,h] .<= flow[s,z1,h] .<= SCapacity[z1,h]) #Capacity constraints of lines
            end
            for n in 1:N
                @constraint(model,-genCapacity[z1,n] .<= P_gen[s,z1,n] .<= genCapacity[z1,n]) #Capacity constraints of generators
            end
        end
        @objective(model, Min, OF[z1,it])
        status = optimize!(model)
    end
    for z1 in 1:Z
        for  z2 in 1:Z  #Updating lagrangemultiplier
            lagrangeMulti[it,z1,z2,h]=lagrangeMulti[it-1,z1,h]+ThetaSLi[s,it,z1,h].-ThetaSLj[s,it,z2,z1,h]
            lagrangeMultj[it,z1,z2,h]=lagrangeMultj[it-1,z1,z2,h]+ThetaSLj[s,it,z1,z2,h].-ThetaSLi[s,it,z2,h]
            if (ThetaSLi[s,it,z1,h].-ThetaSLj[s,it,z2,z1,h])^2+(ThetaSLj[s,it,z1,z2,h].-ThetaSLi[s,it,z2,h])^2 .<= epsilon
                FF=FF+1 #Breaking iteration
                if FF==size(SL)[1]
                    @goto escape_label
                end
            end
        end
    end
end
@label escape_label









