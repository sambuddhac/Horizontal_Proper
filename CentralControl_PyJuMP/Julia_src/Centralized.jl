using Pkg
Pkg.activate("GTheoryJulEnv")
#Defining sets
function milp_avg_hr_central(network::Dict, gen::Dict, cand_line::Dict, int_cand_line::Dict, tran::Dict, setup::Dict)
    C=network["sharedCLines"]
    I=network["internalCLines"] #The number of candidate Lines 
    S=network["CountofScenarios"] #Scenarios
    H=network["tranNumber"]+network["sharedELines"]    #Existing Lines
    G=network["genNumber"] #The number of generators
    T=network["Hours"] #Total number of hours considered (whole year)
    N=network["nodeNumber"]
    #Defining model
    if(setup["Solver"]=="Gurobi")
        # Set solver to use Gurobi for MIP or LP problems
        Gurobi.GurobiSolver
            # Optional setup parameters ############################################
            MyFeasibilityTol = 1e-6 # Constraint (primal) feasibility tolerances. See https://www.gurobi.com/documentation/8.1/refman/feasibilitytol.html
            if(haskey(setup, "Feasib_Tol")) MyFeasibilityTol = setup["Feasib_Tol"] end
            MyOptimalityTol = 1e-6 # Dual feasibility tolerances. See https://www.gurobi.com/documentation/8.1/refman/optimalitytol.html#parameter:OptimalityTol
            if(haskey(setup, "Optimal_Tol")) MyOptimalityTol = setup["Optimal_Tol"] end
            MyPresolve = -1 	# Controls presolve level. See https://www.gurobi.com/documentation/8.1/refman/presolve.html
            if(haskey(setup, "Pre_Solve")) MyPresolve = setup["Pre_Solve"] end
            MyAggFill = -1 		# Allowed fill during presolve aggregation. See https://www.gurobi.com/documentation/8.1/refman/aggfill.html#parameter:AggFill
            if(haskey(setup, "AggFill")) MyAggFill = setup["AggFill"] end
            MyPreDual = -1		# Presolve dualization. See https://www.gurobi.com/documentation/8.1/refman/predual.html#parameter:PreDual
            if(haskey(setup, "PreDual")) MyPreDual = setup["PreDual"] end
            MyTimeLimit = Inf	# Limits total time solver. See https://www.gurobi.com/documentation/8.1/refman/timelimit.html
            if(haskey(setup, "TimeLimit")) MyTimeLimit = setup["TimeLimit"] end
            MyMIPGap = 1e-4		# Relative (p.u. of optimal) mixed integer optimality tolerance for MIP problems (ignored otherwise). See https://www.gurobi.com/documentation/8.1/refman/mipgap2.html
            if(haskey(setup, "MIPGap")) MyMIPGap = setup["MIPGap"] end
            MyCrossover = -1 	# Barrier crossver strategy. See https://www.gurobi.com/documentation/8.1/refman/crossover.html#parameter:Crossover
            if(haskey(setup, "Crossover")) MyCrossover = setup["Crossover"] end
            MyMethod = -1		# Algorithm used to solve continuous models (including MIP root relaxation). See https://www.gurobi.com/documentation/8.1/refman/method.html
            if(haskey(setup, "Method")) MyMethod = setup["Method"] end
            MyBarConvTol = 1e-8 	# Barrier convergence tolerance (determines when barrier terminates). See https://www.gurobi.com/documentation/8.1/refman/barconvtol.html
            if(haskey(setup, "BarConvTol")) MyBarConvTol = setup["BarConvTol"] end
            MyNumericFocus = 0 	# Numerical precision emphasis. See https://www.gurobi.com/documentation/8.1/refman/numericfocus.html
            if(haskey(setup, "NumericFocus")) MyNumericFocus = setup["NumericFocus"] end
        ########################################################################
        CMod = Model(Gurobi.Optimizer)

        set_optimizer_attribute(CMod, "OptimalityTol", MyOptimalityTol)
        set_optimizer_attribute(CMod, "FeasibilityTol", MyFeasibilityTol),
        set_optimizer_attribute(CMod, "Presolve", MyPresolve)
        set_optimizer_attribute(CMod, "AggFill", MyAggFill)
        set_optimizer_attribute(CMod, "PreDual", MyPreDual)
        set_optimizer_attribute(CMod, "TimeLimit", MyTimeLimit)
        set_optimizer_attribute(CMod, "MIPGap", MyMIPGap)
        set_optimizer_attribute(CMod, "Method", MyMethod)
        set_optimizer_attribute(CMod, "BarConvTol", MyBarConvTol)
        set_optimizer_attribute(CMod, "NumericFocus", MyNumericFocus)

    elseif(setup["Solver"]=="CPLEX")
        # Set solve to use CPLEX for MIP or LP problems
            # Optional setup parameters ############################################
            MyFeasibilityTol = 1e-6 # Constraint (primal) feasibility tolerances. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpRHS.html
            if(haskey(setup, "Feasib_Tol")) MyFeasibilityTol = setup["Feasib_Tol"] end
            MyOptimalityTol = 1e-6 # Dual feasibility tolerances. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpOpt.html
            if(haskey(setup, "Optimal_Tol")) MyOptimalityTol = setup["Optimal_Tol"] end
            MyPresolve = 1 	# Decides if presolve is applied. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/PreInd.html
            if(haskey(setup, "Pre_Solve")) MyPresolve = setup["Pre_Solve"] end
            MyAggFill = 10 		# Allowed fill during presolve aggregation. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/AggFill.html
            if(haskey(setup, "AggFill")) MyAggFill = setup["AggFill"] end
            MyPreDual = 0		# Decides whether presolve should pass the primal or dual linear programming problem to the LP optimization algorithm. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/PreDual.html
            if(haskey(setup, "PreDual")) MyPreDual = setup["PreDual"] end
            MyTimeLimit = 1e+75	# Limits total time solver. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/TiLim.html
            if(haskey(setup, "TimeLimit")) MyTimeLimit = setup["TimeLimit"] end
            MyMIPGap = 1e-4		# Relative (p.u. of optimal) mixed integer optimality tolerance for MIP problems (ignored otherwise). See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
            if(haskey(setup, "MIPGap")) MyMIPGap = setup["MIPGap"] end
            MyCrossover = 0 	# Barrier crossver strategy. See https://www.ibm.com/support/knowledgecenter/hr/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarCrossAlg.html
            if(haskey(setup, "Crossover")) MyCrossover = setup["Crossover"] end
            MyMethod = 0		# Algorithm used to solve continuous models (including MIP root relaxation). See https://www.ibm.com/support/knowledgecenter/de/SSSA5P_12.7.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/LPMETHOD.html
            if(haskey(setup, "Method")) MyMethod = setup["Method"] end
            MyBarConvTol = 1e-8 	# Barrier convergence tolerance (determines when barrier terminates). See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarEpComp.html
            if(haskey(setup, "BarConvTol")) MyBarConvTol = setup["BarConvTol"] end
            MyNumericFocus = 0 	# Numerical precision emphasis. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/NumericalEmphasis.html
            if(haskey(setup, "NumericFocus")) MyNumericFocus = setup["NumericFocus"] end
            MyBarObjRng = 1e+75 	# Numerical precision emphasis. See https://www.ibm.com/support/knowledgecenter/en/SSSA5P_12.5.1/ilog.odms.cplex.help/CPLEX/Parameters/topics/NumericalEmphasis.html
            if(haskey(setup, "BarObjRng")) MyBarObjRng = setup["BarObjRng"] end
            MySolutionType = 2 	# Solution type for LP or QP. See https://www.ibm.com/support/knowledgecenter/hr/SSSA5P_12.8.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/SolutionType.html
            if(haskey(setup, "SolutionType")) MySolutionType = setup["SolutionType"] end
        ########################################################################
        CPLEX.Optimizer
        CMod=Model(CPLEX.Optimizer)
        set_optimizer_attribute(CMod, "CPX_PARAM_EPRHS", MyFeasibilityTol)
        set_optimizer_attribute(CMod, "CPX_PARAM_EPOPT", MyOptimalityTol)
        set_optimizer_attribute(CMod, "CPX_PARAM_AGGFILL", MyAggFill)
        set_optimizer_attribute(CMod, "CPX_PARAM_PREDUAL", MyPreDual)
        set_optimizer_attribute(CMod, "CPX_PARAM_TILIM", MyTimeLimit)
        set_optimizer_attribute(CMod, "CPX_PARAM_EPGAP", MyMIPGap)
        set_optimizer_attribute(CMod, "CPX_PARAM_BARCROSSALG", MyCrossover)
        set_optimizer_attribute(CMod, "CPX_PARAM_LPMETHOD", MyMethod)
        set_optimizer_attribute(CMod, "CPX_PARAM_BAREPCOMP", MyBarConvTol)
        set_optimizer_attribute(CMod, "CPX_PARAM_NUMERICALEMPHASIS", MyNumericFocus)
        set_optimizer_attribute(CMod, "CPX_PARAM_BAROBJRNG", MyBarObjRng)
        set_optimizer_attribute(CMod, "CPX_PARAM_SOLUTIONTYPE", MySolutionType)
    elseif(setup["Solver"]=="CBC")
        # Set solver to use Cbc for MIP problems
        Cbc.Optimizer
        # NEED TO TEST SETUP FOR CBC FOR MIP PROBLEMS
        CMod=Model(solver=CbcSolver(logLevel = 2))
    elseif(setup["Solver"]=="CLP")
        # Set solver to use Clp for LP problems
        Clp.Optimizer
        # NEED TO TES SETUP FOR CLP FOR LP PROBLEMS
        CMod=Model(solver=ClpSolver())
    elseif(setup["Solver"]=="GLPK")
        CMod = Model(GLPK.Optimizer)
    end


    #Defining variables
    @variable(CMod, dv_cand_line_decision[1:C], Bin) #Binary Integer Decision variable for building or not building shared candidate line
    @variable(CMod, dv_int_cand_line_decision[1:I], Bin) #Binary Integer Decision variable for building or not building internal candidate line 
    @variable(CMod, dv_cand_flow_mw[1:S, 1:C, 1:T])  #Power flowing on candidate candLineDecision 
    @variable(CMod, dv_int_cand_flow_mw[1:S, 1:I, 1:T])  #Power flowing on internal candidate candLineDecision
    @variable(CMod, 0 <= dv_phase_angle[1:S, 1:N, 1:T] <= 44/7) #Phase angle decision
    @variable(CMod, 0 <= dv_p_gen[1:S, 1:G, 1:T]) #Power of generator
    #@variable(CMod, F) #Objective function

    @expression(CMod, expTotalCost, sum(sum(prob[s,t].*sum(dv_p_gen[s,:,t].*gen["linCostCoeff"][:])for s in S)for t in T)
                                    .+sum(dv_cand_line_decision[:].*cand_line["costPerCap"][:].*cand_line["interestRate"][:]
                                    .*((ones(C).+cand_line["interestRate"][:]).^cand_line["lifeTime"][:])
                                    ./((ones(C).+cand_line["interestRate"][:]).^cand_line["lifeTime"][:].-ones(C)))
                                    .+sum(dv_int_cand_line_decision[:].*int_cand_line["costPerCap"][:].*int_cand_line["interestRate"][:]
                                    .*((ones(C).+int_cand_line["interestRate"][:]).^int_cand_line["lifeTime"][:])
                                    ./((ones(C).+int_cand_line["interestRate"][:]).^int_cand_line["lifeTime"][:].-ones(C))))
    for t in 1:T #The Hour/time loop
        for s in 1:S #The scenario loop
            for h in 1:H #Transmission lines loop
                @expression(CMod, expFlowMW[s,h,t], (dv_phase_angle[s, tran["fromNode"][h], t] - dv_phase_angle[s, tran["toNode"][h], t])/tran["Reactance"][h]) #Constraint regarding the power flowing on existing lines
                @constraint(CMod, expFlowMW[s,h,t]<= tran["lineLimit"][h])  # Line capacity constraints
                @constraint(CMod, -tran["lineLimit"][h]<= EFlowMW[s,h,t])  # Line capacity constraints
            end
            for k in 1:K
                if CandLine["tranZoneID"][h]==z
                    @constraint(CMod, -10000 * (1-candLineDecision[k,z]) .<= candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./CandLine["Reactance"][k] .- candPhaseAngleTo[s,k]./CandLine["Reactance"][k]]) #Constraint regarding the power flowing on shared existing lines
                    @constraint(CMod, candFlowMW[s,k] .- [candPhaseAngleFrom[s,k]./CandLine["Reactance"][k] .- candPhaseAngleTo[s,k]./CandLine["Reactance"][k]] .<= 10000 * (1-candLineDecision[k,z])) #Constraint regarding the power flowing on shared existing lines
                    @constraint(CMod, candFlowMW[s,k]<= candLineDecision[k,z].*CandLine["lineLimit"][k])
                    @constraint(CMod, -candLineDecision[k,z].*CandLine["lineLimit"][k]<=candFlowMW[s,k])
                end
            end
            for n in N
                if Node["nodeZoneID"][n]==z ####I think we should define a zoneID for each node
                    @constraint(CMod, sum(P_gen[s,g] for g in G if Gen["genNodeID"][g]==n).- P_demand[s,n].== sum(candFlowMW[s,k] for k in K if CandLine["fromnode"][k]==n)-sum(candFlowMW[s,k] for k in K if CandLine["tonode"][k]==n)+sum(EFlowMW[s,h] for h in H if Tran["fromnode"][h]==n)-sum(EFlowMW[s,h] for h in H if Tran["tonode"][h]==n))
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
    @objective(CMod, Min, F)
    status = optimize!(CMod)
    @show(status, value.(candLineDecision))
end