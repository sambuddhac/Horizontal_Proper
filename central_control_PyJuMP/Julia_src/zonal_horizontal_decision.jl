function model_var_creation()
   Mod3 = Model(GLPK.Optimizer) #use different solver e.g. IPOPT, GLPK
   @variable(Mod3,0 <= gen_var_1[1:14,1:nrow(scen_prob)]) # generation at node 1
   @variable(Mod3,0 <= gen_var_2[1:30,1:nrow(scen_prob)]) # generation at node 2
   @variable(Mod3,0 <= gen_var_3[1:5,1:nrow(scen_prob)]) # generation at node 3

   @variable(Mod3,shared_line_decision_var[1:nrow(shared_cand)], Bin) #Decision variable for shared candidate lines

   @variable(Mod3,int_line_decision_var_1[1:nrow(int_c(1))], Bin) #Decision variable for internal candidate lines of zone 1 (Binary)
   @variable(Mod3,int_line_decision_var_2[1:nrow(int_c(2))], Bin) #Decision variable for internal candidate lines of zone 2 (Binary)
   @variable(Mod3,int_line_decision_var_3[1:nrow(int_c(3))], Bin) #Decision variable for internal candidate lines of zone 3 (Binary)

   @variable(Mod3,shared_cand_flow[1:nrow(shared_cand),1:nrow(scen_prob)])  #Power flowing on shared candidate lines 

   @variable(Mod3, int_cand_flow_1[1:nrow(int_c(1)),1:nrow(scen_prob)]) #Power flowing on internal candidate lines of zone 1 
   @variable(Mod3, int_cand_flow_2[1:nrow(int_c(2)),1:nrow(scen_prob)]) #Power flowing on internal candidate lines of zone 2
   @variable(Mod3, int_cand_flow_3[1:nrow(int_c(3)),1:nrow(scen_prob)]) #Power flowing on internal candidate lines of zone 3

   @variable(Mod3, 0 <= node_voltage_phase_angle[i=1:nrow(zone_summary),j=1:zone_summary.Nodes_Total[i], 1:nrow(scen_prob)]<= 2*pi)

   @variable(Mod3, shared_ex_flow[1:nrow(shared_ex),1:nrow(scen_prob)])  #Power flowing on shared existing shared lines

   @variable(Mod3, int_ex_flow_1[1:nrow(int_e(1)),1:nrow(scen_prob)])  #Power flowing on internal existing lines of zone 1 
   @variable(Mod3, int_ex_flow_2[1:nrow(int_e(2)),1:nrow(scen_prob)])  #Power flowing on internal existing lines of zone 2
   @variable(Mod3, int_ex_flow_3[1:nrow(int_e(3)),1:nrow(scen_prob)])  #Power flowing on internal existing lines of zone 3
end

function objective_function_creation()
   @expression(Mod3,total_cost,0)

   @expression(Mod3, gen_costNS1[n=1:14,s=1:nrow(scen_prob)], (gen_var_1[n,s] .* sum((g(1).gNodeID .== n) .* MC(1))))
   @expression(Mod3, gen_cost_prob_weightedS1[s=1:nrow(scen_prob)], (scen_weight(s)[1]).*sum(gen_costNS1[n,s] for n in 1:14))
   @expression(Mod3, gen_cost_expected1, sum(gen_cost_prob_weightedS1[s] for s in 1:nrow(scen_prob)))
   Mod3[:total_cost] += gen_cost_expected1
   @expression(Mod3, gen_costNS2[n=1:30,s=1:nrow(scen_prob)], (gen_var_2[n,s] .* sum((g(2).gNodeID .== n) .* MC(2))))
   @expression(Mod3, gen_cost_prob_weightedS2[s=1:nrow(scen_prob)], (scen_weight(s)[1]).*sum(gen_costNS2[n,s] for n in 1:30))
   @expression(Mod3, gen_cost_expected2, sum(gen_cost_prob_weightedS2[s] for s in 1:nrow(scen_prob)))
   Mod3[:total_cost]+=gen_cost_expected2
   @expression(Mod3, gen_costNS3[n=1:5,s=1:nrow(scen_prob)], (gen_var_3[n,s] .* sum((g(3).gNodeID .== n) .* MC(3))))
   @expression(Mod3, gen_cost_prob_weightedS3[s=1:nrow(scen_prob)], (scen_weight(s)[1]).*sum(gen_costNS3[n,s] for n in 1:5))
   @expression(Mod3, gen_cost_expected3, sum(gen_cost_prob_weightedS3[s] for s in 1:nrow(scen_prob)))
   Mod3[:total_cost]+=gen_cost_expected3
   @expression(Mod3, investment_cost , (sum(shared_line_decision_var[c] .* shared_cand.costPerCap[c] .* shared_cand.interestRate[c] 
            .*((1 + shared_cand.interestRate[c]) .^ shared_cand.lifeTime[c]) ./ (((1 + shared_cand.interestRate[c]) .^ shared_cand.lifeTime[c])-1) for c in 1:nrow(shared_cand))
            .+ sum(int_line_decision_var_1[c] .* int_c(1).costPerCap[c] .* int_c(1).interestRate[c] 
            .*((1 + int_c(1).interestRate[c]) .^ int_c(1).lifeTime[c]) ./ (((1 + int_c(1).interestRate[c]) .^ int_c(1).lifeTime[c])-1) for c in 1:nrow(int_c(1)))
            .+ sum(int_line_decision_var_2[c] .* int_c(2).costPerCap[c] .* int_c(2).interestRate[c] 
            .*((1 + int_c(2).interestRate[c]) .^ int_c(2).lifeTime[c]) ./ (((1 + int_c(2).interestRate[c]) .^ int_c(2).lifeTime[c])-1) for c in 1:nrow(int_c(2)))
            .+ sum(int_line_decision_var_3[c] .* int_c(3).costPerCap[c] .* int_c(3).interestRate[c] 
            .*((1 + int_c(3).interestRate[c]) .^ int_c(3).lifeTime[c]) ./ (((1 + int_c(3).interestRate[c]) .^ int_c(3).lifeTime[c])-1) for c in 1:nrow(int_c(3)))))
   Mod3[:total_cost]+=investment_cost
end

function zonal_power_balance()
   #The first zone
   for n in 1:14
      for s in 1:nrow(scen_prob)
         # Power balance constraint for each node
         @constraint(Mod3, sum(g(1).gNodeID .== n) .* gen_var_1[n,s] .+ sum(l(1,s)[:,3].* (l(1,s).lNodeID .== n))./100 .==
         sum((shared_cand.tNodeID1 .== n) .*bin_c(1) .* shared_cand_flow[:,s]) .- sum((shared_cand.tNodeID2 .== n) .* bin_c(1) .* shared_cand_flow[:,s]) .+
         sum((shared_ex.tNodeID1 .== n) .* bin_e(1) .* shared_ex_flow[:,s]) .- sum((shared_ex.tNodeID2 .== n) .* bin_e(1) .* shared_ex_flow[:,s]) .+
         sum((int_c(1).tNodeID1 .== n) .* int_cand_flow_1[:,s]) .- sum((int_c(1).tNodeID2 .== n) .* int_cand_flow_1[:,s]) .+
         sum((int_e(1).tNodeID1 .== n) .* int_ex_flow_1[:,s]) .- sum((int_e(1).tNodeID2 .== n) .* int_ex_flow_1[:,s]))
         #Lower limit for generation of each node
         @constraint(Mod3, sum(g(1).gNodeID .== n) .* gen_var_1[n,s] .<= sum((g(1).gNodeID .== n) .* g(1).PgMax)./100)
         #Upper limit for generation of each node
         @constraint(Mod3, sum((g(1).gNodeID .== n) .* g(1).PgMin)./100 .<= sum(g(1).gNodeID .== n) .* gen_var_1[n,s])
      end
   end

   #The second zone
   for n in 1:30
      for s in 1:nrow(scen_prob)
         # Power balance constraint for each node
         @constraint(Mod3, sum(g(2).gNodeID .== n) .* gen_var_2[n,s] .+ sum(l(2,s)[:,3] .* (l(2,s).lNodeID .== n))./100 .==
         sum((shared_cand.tNodeID1 .== n) .*bin_c(2) .* shared_cand_flow[:,s]) .- sum((shared_cand.tNodeID2 .== n) .* bin_c(2) .* shared_cand_flow[:,s]) .+
         sum((shared_ex.tNodeID1 .== n) .* bin_e(2) .* shared_ex_flow[:,s]) .- sum((shared_ex.tNodeID2 .== n) .* bin_e(2) .* shared_ex_flow[:,s]) .+
         sum((int_c(2).tNodeID1 .== n) .* int_cand_flow_2[:,s]) .- sum((int_c(2).tNodeID2 .== n) .* int_cand_flow_2[:,s]) .+
         sum((int_e(2).tNodeID1 .== n) .* int_ex_flow_2[:,s]) .- sum((int_e(2).tNodeID2 .== n) .* int_ex_flow_2[:,s]))
         #Lower limit for generation of each node
         @constraint(Mod3, sum(g(2).gNodeID .== n) .* gen_var_2[n,s] .<= sum((g(2).gNodeID .== n) .* g(2).PgMax)./100)
         #Upper limit for generation of each node
         @constraint(Mod3, sum((g(2).gNodeID .== n) .* g(2).PgMin)./100 .<= sum(g(2).gNodeID .== n) .* gen_var_2[n,s])
      end
   end

   #The third zone
   for n in 1:5
      for s in 1:nrow(scen_prob)
         # Power balance constraint for each node
         @constraint(Mod3, sum(g(3).gNodeID .== n) .* gen_var_3[n,s] .+ sum(l(3,s)[:,3] .* (l(3,s).lNodeID .== n))./100 .==
         sum((shared_cand.tNodeID1 .== n) .*bin_c(3) .* shared_cand_flow[:,s]) .- sum((shared_cand.tNodeID2 .== n) .* bin_c(3) .* shared_cand_flow[:,s]) .+
         sum((shared_ex.tNodeID1 .== n) .* bin_e(3) .* shared_ex_flow[:,s]) .- sum((shared_ex.tNodeID2 .== n) .* bin_e(3) .* shared_ex_flow[:,s]) .+
         sum((int_c(3).tNodeID1 .== n) .* int_cand_flow_3[:,s]) .- sum((int_c(3).tNodeID2 .== n) .* int_cand_flow_3[:,s]) .+
         sum((int_e(3).tNodeID1 .== n) .* int_ex_flow_3[:,s]) .- sum((int_e(3).tNodeID2 .== n) .* int_ex_flow_3[:,s]))
         #Lower limit for generation of each node
         @constraint(Mod3, sum(g(3).gNodeID .== n) .* gen_var_3[n,s] .<= sum((g(3).gNodeID .== n) .* g(3).PgMax)./100)
         #Upper limit for generation of each node
         @constraint(Mod3, sum((g(3).gNodeID .== n) .* g(3).PgMin)./100 .<= sum(g(3).gNodeID .== n) .* gen_var_3[n,s])
      end
   end

end

function shared_candidate_line_constraint()
   #Shared candidate lines
   M = 5 # this should be changed/translated to PU value...
   for c in 1:nrow(shared_cand)
      for s in 1:nrow(scen_prob)
         @constraint(Mod3,-M .* (1 .- shared_line_decision_var[c]) .<= shared_cand_flow[c,s] .- ((1 ./ shared_cand.reacT[c]) .* (node_voltage_phase_angle[shared_cand.nodeZone1[c],shared_cand.tNodeID1[c],s] .- node_voltage_phase_angle[shared_cand.nodeZone2[c],shared_cand.tNodeID2[c],s])))
         @constraint(Mod3, shared_cand_flow[c,s] .- ((1 ./ shared_cand.reacT[c]) .* (node_voltage_phase_angle[shared_cand.nodeZone1[c],shared_cand.tNodeID1[c],s] .- node_voltage_phase_angle[shared_cand.nodeZone2[c],shared_cand.tNodeID2[c],s])) .<= M .* (1 .- shared_line_decision_var[c]))
         #limiting the upper bound of power flow flowing within candidate shared lines
         @constraint(Mod3, shared_cand_flow[c,s] .<= shared_line_decision_var[c] .*shared_cand.ptMax[c]./100)
         #Limiting the lower bound of power flowing within the candidate shared lines
         @constraint(Mod3, -(shared_cand.ptMax[c]./100) .* shared_line_decision_var[c] .<= shared_cand_flow[c,s])
      end
   end
end

function shared_existing_line_constraint()
   #Shared existing lines
   for h in 1:nrow(shared_ex)
      for s in 1:nrow(scen_prob)
         @constraint(Mod3, shared_ex_flow[h,s] .== (1 ./ shared_ex.reacT[h]) .* (node_voltage_phase_angle[shared_ex.nodeZone1[h],shared_ex.tNodeID1[h],s] .- node_voltage_phase_angle[shared_ex.nodeZone2[h],shared_ex.tNodeID2[h],s]))
         @constraint(Mod3, shared_ex_flow[h,s] .<= shared_ex. ptMax[h]./100)
         @constraint(Mod3, -shared_ex.ptMax[h]./100 .<= shared_ex_flow[h,s])
      end
   end
end

function internal_candidate_line_constraint()
   #Zone 1 internal candidate lines
   for c in 1:nrow(int_c(1))
      for s in 1:nrow(scen_prob)
         @constraint(Mod3,-M .* (1 .- int_line_decision_var_1[c]) .<= int_cand_flow_1[c,s] .- ((1 ./ int_c(1).reacT[c]) .* (node_voltage_phase_angle[int_c(1).zoneNum[c],int_c(1).tNodeID1[c],s] .- node_voltage_phase_angle[int_c(1).zoneNum[c],int_c(1).tNodeID2[c],s])))
         @constraint(Mod3, int_cand_flow_1[c,s] .- ((1 ./ int_c(1).reacT[c]) .* (node_voltage_phase_angle[int_c(1).zoneNum[c],int_c(1).tNodeID1[c],s] .- node_voltage_phase_angle[int_c(1).zoneNum[c],int_c(1).tNodeID2[c],s])) .<= M .* (1 .- int_line_decision_var_1[c]))
         @constraint(Mod3, int_cand_flow_1[c,s] .<= int_line_decision_var_1[c] .* int_c(1).ptMax[c]./100)
         @constraint(Mod3, -(int_c(1).ptMax[c]./100) .*int_line_decision_var_1[c] .<= int_cand_flow_1[c,s])
      end
   end

   #Zone 2 internal candidate lines
   for c in 1:nrow(int_c(2))
      for s in 1:nrow(scen_prob)
         @constraint(Mod3,-M .* (1 .- int_line_decision_var_2[c]) .<= int_cand_flow_2[c,s] .- ((1 ./ int_c(2).reacT[c]) .* (node_voltage_phase_angle[int_c(2).zoneNum[c],int_c(2).tNodeID1[c],s] .- node_voltage_phase_angle[int_c(2).zoneNum[c],int_c(2).tNodeID2[c],s])))
         @constraint(Mod3, int_cand_flow_2[c,s] .- ((1 ./ int_c(2).reacT[c]) .* (node_voltage_phase_angle[int_c(2).zoneNum[c],int_c(2).tNodeID1[c],s] .- node_voltage_phase_angle[int_c(2).zoneNum[c],int_c(2).tNodeID2[c],s])) .<= M .* (1 .- int_line_decision_var_2[c]))
         @constraint(Mod3, int_cand_flow_2[c,s] .<= int_line_decision_var_2[c] .* int_c(2).ptMax[c]./100)
         @constraint(Mod3, -(int_c(2).ptMax[c]./100) .*int_line_decision_var_2[c] .<= int_cand_flow_2[c,s])
      end
   end

   #Zone 3 internal candidate lines
   for c in 1:nrow(int_c(3))
      for s in 1:nrow(scen_prob)
         @constraint(Mod3,-M .* (1 .- int_line_decision_var_3[c]) .<= int_cand_flow_3[c,s] .- ((1 ./ int_c(3).reacT[c]) .* (node_voltage_phase_angle[int_c(3).zoneNum[c],int_c(3).tNodeID1[c],s] .- node_voltage_phase_angle[int_c(3).zoneNum[c],int_c(3).tNodeID2[c],s])))
         @constraint(Mod3, int_cand_flow_3[c,s] .- ((1 ./ int_c(3).reacT[c]) .* (node_voltage_phase_angle[int_c(3).zoneNum[c],int_c(3).tNodeID1[c],s] .- node_voltage_phase_angle[int_c(3).zoneNum[c],int_c(3).tNodeID2[c],s])) .<= M .* (1 .- int_line_decision_var_3[c]))
         @constraint(Mod3, int_cand_flow_3[c,s] .<= int_line_decision_var_3[c] .* int_c(3).ptMax[c]./100)
         @constraint(Mod3, -(int_c(3).ptMax[c]./100) .*int_line_decision_var_3[c] .<= int_cand_flow_3[c,s])
      end
   end

end

function internal_existing_line_constraint()
   #Zone 1 internal existing lines
   for h in 1:nrow(int_e(1))
      for s in 1:nrow(scen_prob)
         @constraint(Mod3, int_ex_flow_1[h,s] .== (1 ./ int_e(1).reacT[h]) .* (node_voltage_phase_angle[int_e(1).zoneNum[h],int_e(1).tNodeID1[h],s] .- node_voltage_phase_angle[int_e(1).zoneNum[h],int_e(1).tNodeID2[h],s]))
         @constraint(Mod3, int_ex_flow_1[h,s] .<= int_e(1).ptMax[h]./100)
         @constraint(Mod3, -int_e(1).ptMax[h]./100 .<= int_ex_flow_1[h,s])
      end
   end

   #Zone 2 internal existing lines
   for h in 1:nrow(int_e(2))
      for s in 1:nrow(scen_prob)
         @constraint(Mod3, int_ex_flow_2[h,s] .== (1 ./ int_e(2).reacT[h]) .* (node_voltage_phase_angle[int_e(2).zoneNum[h],int_e(2).tNodeID1[h],s] .- node_voltage_phase_angle[int_e(2).zoneNum[h],int_e(2).tNodeID2[h],s]))
         @constraint(Mod3, int_ex_flow_2[h,s] .<= int_e(2).ptMax[h]./100)
         @constraint(Mod3, -int_e(2).ptMax[h]./100 .<= int_ex_flow_2[h,s])
      end
   end

   #Zone 3 internal existing lines
   for h in 1:nrow(int_e(3))
      for s in 1:nrow(scen_prob)
         @constraint(Mod3, int_ex_flow_3[h,s] .== (1 ./ int_e(3).reacT[h]) .* (node_voltage_phase_angle[int_e(3).zoneNum[h],int_e(3).tNodeID1[h],s] .- node_voltage_phase_angle[int_e(3).zoneNum[h],int_e(3).tNodeID2[h],s]))
         @constraint(Mod3, int_ex_flow_3[h,s] .<= int_e(3).ptMax[h]./100)
         @constraint(Mod3, -int_e(3).ptMax[h]./100 .<= int_ex_flow_3[h,s])
      end
   end

end

function solve_centralized(Mod3::Model)
   @objective(Mod3, Min, Mod3[:total_cost])
   optimize!(Mod3)
end

function write_ouput()
   if termination_status(Mod3) == MOI.OPTIMAL
      shared_line_decision = value.(shared_line_decision_var)
      shared_cand_power = (value.(shared_cand_flow)) .* 100
      node_voltage_phase_angle = value.(node_voltage_phase_angle)
      shared_ex_power = (value.(shared_ex_flow)) .* 100
      int_line_decision_1 = value.(int_line_decision_var_1)
      int_line_decision_2 = value.(int_line_decision_var_2)
      int_line_decision_3 = value.(int_line_decision_var_3)
      int_line_flow_1 = (value.(int_cand_flow_1)) .* 100
      int_line_flow_2 = (value.(int_cand_flow_2)) .* 100
      int_line_flow_3 = (value.(int_cand_flow_3)) .* 100
      generation_1 = (value.(gen_var_1)) .* 100
      generation_2 = (value.(gen_var_2)) .* 100
      generation_3 = (value.(gen_var_3)) .* 100
      internal_existing_line_flow_1 = (value.(int_ex_flow_1)) .* 100
      internal_existing_line_flow_2 = (value.(int_ex_flow_2)) .* 100
      internal_existing_line_flow_3 = (value.(int_ex_flow_3)) .* 100
      obj_value = objective_value(Mod3)
   end


