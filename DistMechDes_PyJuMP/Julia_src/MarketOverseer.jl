using Pkg
Pkg.activate("GTheoryJulEnv")
function milp_marketOverseer(network::Dict, sharedCLines::Dict, sharedELines::Dict, lagrangeMultPi::Dict, lagrangeMultXi::Dict, setup::Dict)
    K=network["sharedCLines"] #The number of shared candidate lines
    S=network["CountofScenarios"] #Scenarios
    H=network["sharedELines"]#Shared existing lines
    N=network["nodeNumber"] #Total number of nodes
    T=network["Hours"] #Total number of hours considered (whole year)
    zoneList = Array{Int8}(undef, Z) #Array of zones ####I think we should define an array of zones or for example create a dictionary file for ot
    Z=size(zoneList) #Total number of zones 
    flag=1 #If flag is 1 variable candLineDecision is binary, otherwise it is not (Determining relaxed LP or MIP)
    
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
        MOMod = Model(Gurobi.Optimizer)

        set_optimizer_attribute(MOMod, "OptimalityTol", MyOptimalityTol)
        set_optimizer_attribute(MOMod, "FeasibilityTol", MyFeasibilityTol),
        set_optimizer_attribute(MOMod, "Presolve", MyPresolve)
        set_optimizer_attribute(MOMod, "AggFill", MyAggFill)
        set_optimizer_attribute(MOMod, "PreDual", MyPreDual)
        set_optimizer_attribute(MOMod, "TimeLimit", MyTimeLimit)
        set_optimizer_attribute(MOMod, "MIPGap", MyMIPGap)
        set_optimizer_attribute(MOMod, "Method", MyMethod)
        set_optimizer_attribute(MOMod, "BarConvTol", MyBarConvTol)
        set_optimizer_attribute(MOMod, "NumericFocus", MyNumericFocus)

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
        MOMod=Model(CPLEX.Optimizer)
        set_optimizer_attribute(MOMod, "CPX_PARAM_EPRHS", MyFeasibilityTol)
        set_optimizer_attribute(MOMod, "CPX_PARAM_EPOPT", MyOptimalityTol)
        set_optimizer_attribute(MOMod, "CPX_PARAM_AGGFILL", MyAggFill)
        set_optimizer_attribute(MOMod, "CPX_PARAM_PREDUAL", MyPreDual)
        set_optimizer_attribute(MOMod, "CPX_PARAM_TILIM", MyTimeLimit)
        set_optimizer_attribute(MOMod, "CPX_PARAM_EPGAP", MyMIPGap)
        set_optimizer_attribute(MOMod, "CPX_PARAM_BARCROSSALG", MyCrossover)
        set_optimizer_attribute(MOMod, "CPX_PARAM_LPMETHOD", MyMethod)
        set_optimizer_attribute(MOMod, "CPX_PARAM_BAREPCOMP", MyBarConvTol)
        set_optimizer_attribute(MOMod, "CPX_PARAM_NUMERICALEMPHASIS", MyNumericFocus)
        set_optimizer_attribute(MOMod, "CPX_PARAM_BAROBJRNG", MyBarObjRng)
        set_optimizer_attribute(MOMod, "CPX_PARAM_SOLUTIONTYPE", MySolutionType)
    elseif(setup["Solver"]=="CBC")
        # Set solver to use Cbc for MIP problems
        Cbc.Optimizer
        # NEED TO TEST SETUP FOR CBC FOR MIP PROBLEMS
        MOMod=Model(solver=CbcSolver(logLevel = 2))
    elseif(setup["Solver"]=="CLP")
        # Set solver to use Clp for LP problems
        Clp.Optimizer
        # NEED TO TES SETUP FOR CLP FOR LP PROBLEMS
        CMod=Model(solver=ClpSolver())
    elseif(setup["Solver"]=="GLPK")
        MOMod = Model(GLPK.Optimizer)
    end

     
    @variable(MOMod, 0 <= candLineDecision[1:K] <= 1) #Decision variable 
    @variable(MOMod, 0 <= CSFlowMW[1:S, 1:K, 1:T])  #Power flowing on shared candidate candLineDecision 
    @variable(MOMod, 0 <= CSPhaseAngle[1:S, 1:K, 1:N, 1:T]) #Phase angle decision for  candidate shared lines
    @variable(MOMod, 0 <= ESFlowMW[1:S, 1:H, 1:T])  #Power flowing on existing shared lines
    @variable(MOMod, 0 <= ESPhaseAngle[1:S, 1:H, 1:N, 1:T]) #Phase angle decision for existing shared lines 
    @variable(MOMod, F[1:S,1:Z])
    for t in T
        for z in zoneList
            for s in 1:S
                for k in 1:K
                    for h in 1:H
                        @expression(MOMod, F[1:S,1:Z], -[sum(lagrangeMultPi[z,:].*candLineDecision[:]) 
                            .+ sum(lagrangeMultXi[z,n,s].*CSPhaseAngleFrom[s,k,n,t] for n in union(sharedCLines["fromNode"][k],sharedCLines["toNode"][k]))       ####I am not sure if langrange multipliers should have time index
                            .+ sum(lagrangeMultXi[z,n,s].*ESPhaseAngleFrom[s,h,n,t] for n in union(sharedELines["fromNode"][h],sharedELines["toNode"][h])) 
                            
                    end
                end
    
            
                for h in 1:H
                    @constraint(MOMod, ESFlowMW[s,h,t] .== (sum(ESPhaseAngle[s,h,n,t] for n in sharedELines["fromNode"][h]) 
                                        .- sum(ESPhaseAngle[s,h,n,t] for n in sharedELines["toNode"][h])) ./sharedELines["Reactance"][h]) #Constraint regarding the power flowing on shared existing lines
                    @constraint(MOMod, ESFlowMW[s,h,t]<= sharedELines["lineLimit"][h])  # Capacity constraints
                end
                for k in 1:K
                    if flag==1
                        @constraint(MOMod, candLineDecision[k] in MOI.Integer())
                    end
                    @constraint(MOMod, -10000 * (1-candLineDecision[k]) .<= CSFlowMW[s,k] 
                                            .- [(sum(CSPhaseAngle[s,k,n,t] for n in sharedCLines["fromNode"][k]) 
                                        .- sum(CSPhaseAngle[s,k,n,t] for n in sharedCLines["toNode"][k])) ./sharedCLines["Reactance"][k])] #Constraint regarding the power flowing on shared candidate lines
                    @constraint(MOMod, CSFlowMW[s,k] .- [(sum(CSPhaseAngle[s,k,n,t] for n in sharedCLines["fromNode"][k]) 
                                        .- sum(CSPhaseAngle[s,k,n,t] for n in sharedCLines["toNode"][k])) 
                                                    ./sharedCLines["Reactance"][k])] .<= 10000 * (1-candLineDecision[k])) #Constraint regarding the power flowing on shared candidate lines
                    @constraint(MOMod, CSFlowMW[s,k]<= sharedCLines["lineLimit"][k])
                end
            end
        end
    end
    @objective(MOMod, Min, sum(F[:,:]))
    status = optimize!(MOMod)
    @show(status, value.(candLineDecision))
    end
