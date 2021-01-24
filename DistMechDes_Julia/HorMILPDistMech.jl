#=
    #HorMILPDistMech() for zonal decision making for building cross-border and internal transmission lines and optimal operation

    #Author: Sambuddha Chakrabarti
    #This is the zonal decision making for building cross-border and internal transmission lines and optimal operation
=#

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

function HorMILPDistMech(coordInstanceRef, LagMultXi, LagMultPi, totalCandLineNum, totalSharedNodeNum)
	start_t = now()
    	if solChoice == 1
        	model = Model(with_optimizer(Gurobi.Optimizer, OUTPUTLOG=OUTPUTLOG, MAXTIME=-MAXTIME))
    	elseif solChoice == 2
        	model = Model(with_optimizer(GLPK.Optimizer, OUTPUTLOG=OUTPUTLOG, MAXTIME=-MAXTIME))
    	elseif solChoice == 3
        	model = Model(with_optimizer(MathOptInterfaceMosek.Optimizer, OUTPUTLOG=OUTPUTLOG, MAXTIME=-MAXTIME))
    	elseif solChoice == 4
        	model = Model(with_optimizer(Cbc.Optimizer, OUTPUTLOG=OUTPUTLOG, MAXTIME=-MAXTIME))
    	elseif solChoice == 5
        	model = Model(with_optimizer(Ipopt.Optimizer, OUTPUTLOG=OUTPUTLOG, MAXTIME=-MAXTIME))
    	else
        	error("Invalid Solver Choice:", solChoice)
    	end
	#C++ starts
        start_t = now()
	vector<int>::iterator diffZNIt; // Iterator for diffZoneNodeID
	vector<Powergenerator*>::iterator genIterator; // Iterator for Powergenerator objects
	vector<transmissionLine*>::iterator tranIterator; // Iterator for Transmission line objects
	vector<Load*>::iterator loadIterator; // Iterator for load objects
	vector<Node*>::iterator nodeIterator; // Iterator for node objects
	vector<candLine*>::iterator candIterator; // Iterator for candidate lines
	vector<SELine*>::iterator exsharedIterator; // Iterator for shared existing lines
        vector<intCandLine*>::iterator intCandIterator; // Iterator for candidate lines

	string outSummaryFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPK/OutSummaryGLPK" + to_string(zonalIndex) + ".txt";
	ofstream outPutFile(outSummaryFileName, ios::out); // Create Output File to output the Summary of Results
	if (!outPutFile){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	# C++ ends
        @variables model begin
                0 <= Pgen[1:countOfScenarios, 1:genNumber] #power generation MW outputs
                voltTheta[1:countOfScenarios, 1:nodeNumber] #voltage phase angles for internal zonal nodes
                sharedCandFlow[1:countOfScenarios, 1:sharedCLines] #power flow values for shared candidate lines
                sharedCandBinDec[1:sharedCLines], Bin #binary integer decision variable values for shared candidate lines
                otherTheta[1:countOfScenarios, 1:otherNodeCount] #voltage phase angles of other-zone nodes connected through shared existing and candidate lines
                internalCandFlow[1:countOfScenarios, 1:internalCLines] #power flow values for internal candidate lines
                internalCandBinDec[1:internalCLines], Bin #binary integer decision variable values for internal candidate lines
        end
	# Supply demand balances
    	@constraint(DCOPF, cBalance[i in N], 
        sum(GEN[g] for g in gens[gens.connnode .== i,:connnode]) 
            + sum(load for load in loads[loads.connnode .== i,:demand]) 
        == sum(FLOW[i,j] for j in lines[lines.fromnode .== i,:tonode])
    	)
        @constraints model begin
                #Constraints corresponding to supply-demand balance
		#C++ starts #### Dear Sam, I copied the following lines for marketoverseer optimization
                rCount = 1 #Initialize the row count
	        for scenCounter in 1:countOfScenarios
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			outPutFile << "\nGeneration\t" << rCount << "\n";
			int genListLength = (*nodeIterator)->getGenLength(); // get the number
			cout << "\nScenario Count: " << scenCounter << " node count: " << (*nodeIterator)->getNodeID() << " Conn Gen: " << genListLength << endl;
			for (int cCount = 1; cCount <= genListLength; ++cCount){
				cout << "\nSerial: " << cCount << " Generator Serial: " << (*nodeIterator)->getGenSer(cCount) << endl;
				ia.push_back(rCount), ja.push_back(scenCounter*genNumber+(*nodeIterator)->getGenSer(cCount)), ar.push_back(1.0);
				outPutFile << "\n" << rCount << "\t" << scenCounter*genNumber+(*nodeIterator)->getGenSer(cCount) << "\t" << 1.0 << endl;
				++index;
			}
			outPutFile << "\nIntrazonal Node Angles\t" << rCount << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+rCount), ar.push_back(((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+rCount << "\t" << ((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()) << "\t" << ((*nodeIterator)->getFromReact()) << ((*nodeIterator)->getToReact()) << endl;
			++index;
			outPutFile << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
			int connNodeListLength = (*nodeIterator)->getConNodeLength(); // get the number of intra-zonal nodes connected to this node
			for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
				if (((*nodeIterator)->getConnReact(cCount))<=0) {
					ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount))), ar.push_back(-((*nodeIterator)->getConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";
                                }
				else {
					ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount))), ar.push_back(((*nodeIterator)->getConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";
                                }
				++index;

			}
			outPutFile << "\nConnected Outer zonal Node Angles\t" << rCount << "\n";
			int connOutNodeLength = (*nodeIterator)->getExtraNodeLength(); // get the number of outer-zonal nodes connected to this node
			for (int cCount = 1; cCount <= connOutNodeLength; ++cCount){
                                if (((*nodeIterator)->getExtConnReact(cCount))<=0) {
					ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount)), ar.push_back(-((*nodeIterator)->getExtConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (-((*nodeIterator)->getExtConnReact(cCount))) << "\n";
				}
                                else {
					ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount)), ar.push_back(((*nodeIterator)->getExtConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (((*nodeIterator)->getExtConnReact(cCount))) << "\n";
				}
				++index;
			}
			outPutFile << "\nConnected Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connCandListLengthF = (*nodeIterator)->getCandLineLengthF(); // get the number of candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connCandListLengthF; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerF(cCount)), ar.push_back(-1);
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerF(cCount) << "\t" << -1.0 << "\n";
				++index;
			}
			outPutFile << "\nConnected Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connCandListLengthT = (*nodeIterator)->getCandLineLengthT(); // get the number of candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connCandListLengthT; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerT(cCount)), ar.push_back(1);
				outPutFile << "\n" << rCount << "\t" << genNumber+nodeNumber+(*nodeIterator)->getCandSerT(cCount) << "\t" << 1.0 << "\n";
				++index;
			}
                	outPutFile << "\nConnected Internal Candidate Lines for which this is the From node\t" << rCount << "\n";
                	int connintCandListLengthF = (*nodeIterator)->getIntCandLineLengthF(); // get the number of internal candidate lines connected to this from node 
                	for (int cCount = 1; cCount <= connintCandListLengthF; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerF(cCount)), ar.push_back(-1);			
                        	outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerF(cCount) << "\t" << -1.0 << "\n";
				++index;
               		}
               		outPutFile << "\nConnected Internal Candidate Lines for which this is the To node\t" << rCount << "\n";
                	int connintCandListLengthT = (*nodeIterator)->getIntCandLineLengthT(); // get the number of internal candidate lines connected to this to node 
                	for (int cCount = 1; cCount <= connintCandListLengthT; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerT(cCount)), ar.push_back(1);
                        	outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerT(cCount) << "\t" << 1.0 << "\n";
				++index;
                	}
			++rCount; // Increment the row count to point to the next node object
		}
	}

        end
        int dimRow = countOfScenarios*(2 * genNumber + 4 * sharedCLines + 2 * sharedELines + 2 * tranNumber + nodeNumber + 4*internalCLines); // Total number of rows of the A matrix (number of structural constraints of the LP) first term to account for lower and upper generating limits, second term for lower and upper line limits & lower and upper definition limits of candidate shared lines, third term for lower and upper line limits for shared existing lines, fourth term for lower and upper line limits for internal zonal lines, the fifth term to account for nodal power balance constraints, and sixth term to account for the internal candidate lines
        int dimCol = countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount+internalCLines)+sharedCLines+internalCLines

	//Row Definitions and Bounds Corresponding to Constraints/

	/*******************************************************************************************/

	// Constraints corresponding to supply-demand balance
	string outPGenFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputPowerGLPK/OutPowerGenGLPK" + to_string(zonalIndex) + ".txt"; 
	ofstream powerGenOut(outPGenFileName, ios::out);
	if (!powerGenOut){
		cerr << "\nCouldn't open the file" << endl;
		exit(1);
	}
	// Vectors for storing the output data
	vector<int> busCount; // vector for storing the node/bus serial
	outPutFile << "Constraints corresponding to Supply-Demand Balance right hand side" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (int rCount = 1; rCount <= nodeNumber; ++rCount){
			busCount.push_back(rCount);
			glp_set_row_name(milp, scenCounter*nodeNumber+rCount, NULL); // Specify the particular supply-demand balance constraint for the particular node in GLPK
			if (((nodeObject[rCount-1])->devpinitMessage(scenCounter))==0)
				glp_set_row_bnds(milp, scenCounter*nodeNumber+rCount, GLP_FX, ((nodeObject[rCount-1])->devpinitMessage(scenCounter)), 0.0); // Specify the right hand side which is the net demand for the particular node
			else 
				glp_set_row_bnds(milp, scenCounter*nodeNumber+rCount, GLP_FX, -((nodeObject[rCount-1])->devpinitMessage(scenCounter)), 0.0); // Specify the right hand side which is the net demand for the particular node
			outPutFile << "Connected load to node " << rCount << " in scenario " << scenCounter+1 << " is " << (nodeObject[rCount-1])->devpinitMessage(scenCounter)*100 << " MW" << endl;
			outPutFile << rCount << "\t";
			if (((nodeObject[rCount-1])->devpinitMessage(scenCounter))==0)
				outPutFile << ((nodeObject[rCount-1])->devpinitMessage(scenCounter))*100 << " MW" << endl;
			else
				outPutFile << -((nodeObject[rCount-1])->devpinitMessage(scenCounter))*100 << " MW" << endl;
		}
	}
	/*******************************************************************************************/

	// Constraints corresponding to Powergenerator Lower Bound
	outPutFile << "Constraints corresponding to Powergenerator Lower Bound" << endl;
	int genRun = 0;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator) {
			++genRun;
			int rCount = countOfScenarios*nodeNumber+genRun;
			glp_set_row_name(milp, rCount, NULL);
			glp_set_row_bnds(milp, rCount, GLP_LO, ((*genIterator)->getPMin()), 0.0);
			outPutFile << rCount << "\t";
			outPutFile << ((*genIterator)->getPMin())*100 << " MW" << endl;
		}
	}
	// Constraints corresponding to Powergenerator Upper Bound
	outPutFile << "Constraints corresponding to Powergenerator Upper Bound" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator) {
			++genRun;
			int rCount = countOfScenarios*nodeNumber+genRun;
			glp_set_row_name(milp, rCount, NULL);
			glp_set_row_bnds(milp, rCount, GLP_UP, 0.0, ((*genIterator)->getPMax()));
			outPutFile << rCount << "\t";
			outPutFile << ((*genIterator)->getPMax())*100 << " MW" << endl;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Constraints corresponding to Powergenerator Bounds: " << genRun << endl;
	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for internal zonal transmission lines
	outPutFile << "\nTesting of intrazonal transmission line forward flow limits" << endl; 
	int rowCount = countOfScenarios*(nodeNumber+2*genNumber)+1;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (*tranIterator)->getFlowLimit());
			outPutFile << rowCount << "\t";
			outPutFile << ((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for internal zonal transmission lines
	outPutFile << "\nTesting of intrazonal transmission line reverse flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -((*tranIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -((*tranIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for internal zonal transmission lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared existing lines
	outPutFile << "\nTesting of Shared Existing transmission line forward flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, (*exsharedIterator)->getFlowLimit());
			outPutFile << rowCount << "\t";
			outPutFile << ((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared existing lines
	outPutFile << "\nTesting of Shared Existing transmission line reverse flow limits" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -((*exsharedIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -((*exsharedIterator)->getFlowLimit())*100 << " MW" << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared existing lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to Line Forward Flow limits for shared candidate lines
	outPutFile << "\nTesting of Shared Candidate transmission line forward flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for shared candidate lines
	outPutFile << "\nTesting of Shared Candidate transmission line reverse flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for shared candidate lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to shared candidate lines flow definition upper bound
	outPutFile << "\nTesting of Definition of Flows on Shared Candidate transmission lines upper bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 2.5*((*candIterator)->getFlowLimit()));
			outPutFile << rowCount << "\t";
			outPutFile << BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to shared candidate lines flow definition lower bound
	outPutFile << "\nTesting of Definition of Flows on Shared Candidate transmission lines lower bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -2.5*((*candIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for shared candidate lines flow definition: " << rowCount << endl;
	/*******************************************************************************************/
	// Constraints corresponding to Line Forward Flow limits for internal candidate lines
	outPutFile << "\nTesting of Internal Candidate transmission line forward flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to Line Reverse Flow limits for internal candidate lines
	outPutFile << "\nTesting of Internal Candidate transmission line reverse flow limits" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, 0.0, 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << 0.0 << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for Line Flow limits for internal candidate lines: " << rowCount << endl;

	/*******************************************************************************************/

	// Constraints corresponding to internal candidate lines flow definition upper bound
	outPutFile << "\nTesting of Definition of Flows on Internal Candidate transmission lines upper bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_UP, 0.0, 2.5*((*intCandIterator)->getFlowLimit()));
			outPutFile << rowCount << "\t";
			outPutFile << BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << endl;
	// Constraints corresponding to internal candidate lines flow definition lower bound
	outPutFile << "\nTesting of Definition of Flows on Internal Candidate transmission lines lower bound" << endl; 
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_row_name(milp, rowCount, NULL);
			glp_set_row_bnds(milp, rowCount, GLP_LO, -2.5*((*intCandIterator)->getFlowLimit()), 0.0);
			outPutFile << rowCount << "\t";
			outPutFile << -BIGM << endl;
			++rowCount;
		}
	}
	outPutFile << "\nTotal number of rows after accounting for internal candidate lines flow definition: " << rowCount << endl;
	outPutFile << "\nConstraint bounds (rows) Specified" << endl;
	outPutFile << "\nTotal number of rows: " << rowCount - 1 << endl;
	/*******************************************************************************************/

	/*******************************************************************************************/

	/*Column Definitions, Bounds, and Objective Function Co-efficients*/
	glp_add_cols(milp, dimCol);
	int colCount = 1;
	//Columns corresponding to Power Generation continuous variables for different generators//
	outPutFile << "\nCoefficients of Power generator variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_LO, 0.0, 0.0);
			glp_set_obj_coef(milp, colCount, probability.at(scenCounter)*((*genIterator)->getLinCoeff()));
			outPutFile << colCount << "\t";
			outPutFile << ((*genIterator)->getLinCoeff())/100 << " $/MW" << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Power Generation continuous variables for different generators: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Voltage Phase Angles continuous variables for different intrazonal nodes//
	outPutFile << "\nCoefficients of Voltage Phase Angles continuous variables for different intrazonal nodes" << endl;
	outPutFile << "\nVariable Count\tShared Node\tGlobal Rank\tLagMultXiIndex\tLagMultXiValue" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {	
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_DB, 0, (44/7));
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
				glp_set_obj_coef(milp, colCount, LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank())]);
				outPutFile << colCount << "\tYes\t" << ((*nodeIterator)->getGlobalRank()) << "\t" << scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank()) << "\t" << LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+((*nodeIterator)->getGlobalRank())] << endl;	
			}
			else {
				glp_set_obj_coef(milp, colCount, 0.0);
				outPutFile << colCount << "\tNo\t" << "-" << "\t" << "-" << "\t" << "-" << endl;	
			}
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for different intrazonal nodes: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Shared Candidate Line Flows continuous variables//
	outPutFile << "\nCoefficients corresponding to Shared Candidate Line Flows continuous variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(milp, colCount, 0.0);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Flows continuous variables: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Shared Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nCoefficients corresponding to Shared Candidate Line Construction Decision Binary Integer variables" << endl;
	outPutFile << "\nVariable Count\tGlobal Rank\tLagMultPiIndex\tLagMultPiValue\tInvestment Cost" << endl;
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		glp_set_col_kind(milp, colCount, GLP_BV);
		glp_set_col_name(milp, colCount, NULL);
		glp_set_col_bnds(milp, colCount, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(milp, colCount, LagMultPi[(zonalIndex-1)*totalCandLineNum+((*candIterator)->getGlobalRank())]+((*candIterator)->returnOwnership())*0.5*((*candIterator)->getInvestCost()));
		outPutFile << colCount << "\t" << ((*candIterator)->getGlobalRank()) << "\t" << (zonalIndex-1)*totalCandLineNum+((*candIterator)->getGlobalRank()) << "\t" << LagMultPi[(zonalIndex-1)*totalCandLineNum+((*candIterator)->getGlobalRank())] << "\t" << ((*candIterator)->returnOwnership())*0.5*((*candIterator)->getInvestCost()) << endl;	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Shared Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	/*******************************************************************************************/
	/*int diffNodeTestCounter = 0;
	for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
		if (diffNodeTestCounter > 0) { // Skip the first element, since it's a dummy "0"	
			cout << " Column count after the shared node test is " << diffNodeTestCounter << endl;
		}
		++diffNodeTestCounter;
	}*/
	//Columns corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines//
	outPutFile << "\nCoefficients corresponding to Voltage Phase Angles continuous variables for other zone nodes for shared lines" << endl;
	outPutFile << "\nVariable Count\tGlobal Rank\tLagMultXiIndex\tLagMultXiValue" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffNodeCounter = 0; // counter flag to indicate the first element of the diffZoneNodeID list
		for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
			if (diffNodeCounter > 0) { // Skip the first element, since it's a dummy "0"
				glp_set_col_kind(milp, colCount, GLP_CV);
				glp_set_col_name(milp, colCount, NULL);
				glp_set_col_bnds(milp, colCount, GLP_DB, 0, (44/7));
				glp_set_obj_coef(milp, colCount, LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+(*globalIterator)]);
				outPutFile << colCount << (*globalIterator) << "\t" << scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+(*globalIterator) << "\t" << LagMultXi[scenCounter*zonalCount*totalSharedNodeNum+(zonalIndex-1)*totalSharedNodeNum+(*globalIterator)] << endl;		
				++colCount;
				//cout << " Column count after the shared node " << diffNodeCounter << " is " << colCount << endl;
			}
			++diffNodeCounter;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Voltage Phase Angles continuous variables for other zone nodes for shared lines: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Internal Candidate Line Flows continuous variables//
	outPutFile << "\nCoefficients corresponding to Internal Candidate Line Flows continuous variables" << endl;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			glp_set_col_kind(milp, colCount, GLP_CV);
			glp_set_col_name(milp, colCount, NULL);
			glp_set_col_bnds(milp, colCount, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(milp, colCount, 0.0);
			outPutFile << colCount << "\t";
			outPutFile << 0 << endl;
			++colCount;
		}
	}
	outPutFile << "\nTotal number of columns after accounting for Internal Candidate Line Flows continuous variables: " << colCount << endl;
	/*******************************************************************************************/

	//Columns corresponding to Internal Candidate Line Construction Decision Binary Integer variables//
	outPutFile << "\nCoefficients corresponding to Internal Candidate Line Construction Decision Binary Integer variables" << endl;
	outPutFile << "\nVariable Count\tInvestment Cost" << endl;
	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		glp_set_col_kind(milp, colCount, GLP_BV);
		glp_set_col_name(milp, colCount, NULL);
		glp_set_col_bnds(milp, colCount, GLP_DB, 0.0, 1.0);
		glp_set_obj_coef(milp, colCount, ((*intCandIterator)->getInvestCost()));
		outPutFile << colCount << "\t" << ((*intCandIterator)->getInvestCost()) << endl;	
		++colCount;
	}
	outPutFile << "\nTotal number of columns after accounting for Internal Candidate Line Construction Decision Binary Integer variables: " << colCount << endl;
	outPutFile << "\nTotal Number of columns for generation, angles, integer variables, and flows: " << colCount - 1 << endl;
	outPutFile << "\nDecision Variables and Objective Function defined" << endl;
	outPutFile << "\nTotal Number of columns: " << colCount - 1 << endl;
	/*******************************************************************************************/

	/*******************************************************************************************/

	//Non-Zero entries of A matrix (Constraint/Coefficient matrix entries)//
	int index = 1;
	ia.push_back(0), ja.push_back(0), ar.push_back(0.0);
	// Coefficients for the supply-demand balance constraints
	outPutFile << "\nNon-zero elements of A matrix" << endl;
	outPutFile << "\nRow Number\tColumn Number\tNon-zero Entry\tFrom Reactance\tToReactance" << endl;
	outPutFile << "\nCoefficients for the supply-demand balance constraints" << endl;
	int rCount = 1; // Initialize the row count
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			outPutFile << "\nGeneration\t" << rCount << "\n";
			int genListLength = (*nodeIterator)->getGenLength(); // get the number
			cout << "\nScenario Count: " << scenCounter << " node count: " << (*nodeIterator)->getNodeID() << " Conn Gen: " << genListLength << endl;
			for (int cCount = 1; cCount <= genListLength; ++cCount){
				cout << "\nSerial: " << cCount << " Generator Serial: " << (*nodeIterator)->getGenSer(cCount) << endl;
				ia.push_back(rCount), ja.push_back(scenCounter*genNumber+(*nodeIterator)->getGenSer(cCount)), ar.push_back(1.0);
				outPutFile << "\n" << rCount << "\t" << scenCounter*genNumber+(*nodeIterator)->getGenSer(cCount) << "\t" << 1.0 << endl;
				++index;
			}
			outPutFile << "\nIntrazonal Node Angles\t" << rCount << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+rCount), ar.push_back(((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+rCount << "\t" << ((*nodeIterator)->getToReact())-((*nodeIterator)->getFromReact()) << "\t" << ((*nodeIterator)->getFromReact()) << ((*nodeIterator)->getToReact()) << endl;
			++index;
			outPutFile << "\nConnected Intrazonal Node Angles\t" << rCount << "\n";
			int connNodeListLength = (*nodeIterator)->getConNodeLength(); // get the number of intra-zonal nodes connected to this node
			for (int cCount = 1; cCount <= connNodeListLength; ++cCount){
				if (((*nodeIterator)->getConnReact(cCount))<=0) {
					ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount))), ar.push_back(-((*nodeIterator)->getConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";
                                }
				else {
					ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount))), ar.push_back(((*nodeIterator)->getConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber+scenCounter*nodeNumber+((*nodeIterator)->getConnSer(cCount)) << "\t" <<  (((*nodeIterator)->getConnReact(cCount))) << "\n";
                                }
				++index;

			}
			outPutFile << "\nConnected Outer zonal Node Angles\t" << rCount << "\n";
			int connOutNodeLength = (*nodeIterator)->getExtraNodeLength(); // get the number of outer-zonal nodes connected to this node
			for (int cCount = 1; cCount <= connOutNodeLength; ++cCount){
                                if (((*nodeIterator)->getExtConnReact(cCount))<=0) {
					ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount)), ar.push_back(-((*nodeIterator)->getExtConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (-((*nodeIterator)->getExtConnReact(cCount))) << "\n";
				}
                                else {
					ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount)), ar.push_back(((*nodeIterator)->getExtConnReact(cCount)));
					outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*nodeIterator)->getExtConnSer(cCount) << "\t" <<  (((*nodeIterator)->getExtConnReact(cCount))) << "\n";
				}
				++index;
			}
			outPutFile << "\nConnected Candidate Lines for which this is the From node\t" << rCount << "\n";
			int connCandListLengthF = (*nodeIterator)->getCandLineLengthF(); // get the number of candidate lines connected to this from node 
			for (int cCount = 1; cCount <= connCandListLengthF; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerF(cCount)), ar.push_back(-1);
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerF(cCount) << "\t" << -1.0 << "\n";
				++index;
			}
			outPutFile << "\nConnected Candidate Lines for which this is the To node\t" << rCount << "\n";
			int connCandListLengthT = (*nodeIterator)->getCandLineLengthT(); // get the number of candidate lines connected to this to node 
			for (int cCount = 1; cCount <= connCandListLengthT; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber)+scenCounter*sharedCLines+(*nodeIterator)->getCandSerT(cCount)), ar.push_back(1);
				outPutFile << "\n" << rCount << "\t" << genNumber+nodeNumber+(*nodeIterator)->getCandSerT(cCount) << "\t" << 1.0 << "\n";
				++index;
			}
                	outPutFile << "\nConnected Internal Candidate Lines for which this is the From node\t" << rCount << "\n";
                	int connintCandListLengthF = (*nodeIterator)->getIntCandLineLengthF(); // get the number of internal candidate lines connected to this from node 
                	for (int cCount = 1; cCount <= connintCandListLengthF; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerF(cCount)), ar.push_back(-1);			
                        	outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerF(cCount) << "\t" << -1.0 << "\n";
				++index;
               		}
               		outPutFile << "\nConnected Internal Candidate Lines for which this is the To node\t" << rCount << "\n";
                	int connintCandListLengthT = (*nodeIterator)->getIntCandLineLengthT(); // get the number of internal candidate lines connected to this to node 
                	for (int cCount = 1; cCount <= connintCandListLengthT; ++cCount){
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerT(cCount)), ar.push_back(1);
                        	outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines+otherNodeCount)+sharedCLines+scenCounter*internalCLines+(*nodeIterator)->getIntCandSerT(cCount) << "\t" << 1.0 << "\n";
				++index;
                	}
			++rCount; // Increment the row count to point to the next node object
		}
	}
	/*******************************************************************************************/

	// Coefficients corresponding to lower generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			ia.push_back(rCount), ja.push_back(rCount - countOfScenarios*nodeNumber), ar.push_back(1.0);
			++index;
			outPutFile << rCount << "\t" << (rCount - countOfScenarios*nodeNumber) << "\t" << 1.0 << "\t" << (*genIterator)->getPMin() << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	// Coefficients corresponding to upper generation limits
	outPutFile << "\nCoefficients corresponding to lower generation limits\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			ia.push_back(rCount), ja.push_back(rCount - countOfScenarios*(genNumber + nodeNumber)), ar.push_back(1.0);
			++index;
			outPutFile << rCount << "\t" << (rCount - countOfScenarios*(genNumber + nodeNumber)) << "\t" << 1.0 << "\t" << ((*genIterator)->getPMax()) << endl;
			++rCount; // Increment the row count to point to the next generator object
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to intra-zone Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1()), ar.push_back(1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2()), ar.push_back(-1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
			++index;
			++rCount; // Increment the row count to point to the next transmission line object
		
		}
	}
	// Coefficients corresponding to intra-zone Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to intra-zone Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1()), ar.push_back(1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1() << "\t" << 1/((*tranIterator)->getReactance()) << "\t" << 1/((*tranIterator)->getReactance()) << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2()), ar.push_back(-1/((*tranIterator)->getReactance()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2() << "\t" << -1/((*tranIterator)->getReactance()) << "\t" << "-" << "\t" << -1/((*tranIterator)->getReactance()) << "\n";
			++index;
			++rCount; // Increment the row count to point to the next transmission line object
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared existing Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			if ((*exsharedIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	// Coefficients corresponding to shared existing Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared existing Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
			if ((*exsharedIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank()), ar.push_back(1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*exsharedIterator)->getExtNodeRank() << "\t" << 1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID()), ar.push_back(-1/((*exsharedIterator)->getReactance()));
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber + (*exsharedIterator)->getIntlNodeID() << "\t" << -1/((*exsharedIterator)->getReactance()) << "\n";
				++index;
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Forward Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared candidate Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines) << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*sharedCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines), ar.push_back(-((*candIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*sharedCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines << "\t" << -((*candIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
		}
	}
	// Coefficients corresponding to shared candidate Line Reverse Flow Limit Constraints
	outPutFile << "\nCoefficients corresponding to shared candidate Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines) << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines), ar.push_back(((*candIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines)-scenCounter*sharedCLines << "\t" << ((*candIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
		}
	}
	/*******************************************************************************************/
	
	// Coefficients corresponding to shared candidate Line Definition upper bound
	outPutFile << "\nCoefficients corresponding to shared candidate Line Definition upper bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			if ((*candIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines), ar.push_back(2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines << "\t" << BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines), ar.push_back(2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+sharedCLines)-scenCounter*sharedCLines << "\t" << BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	// Coefficients corresponding to shared candidate Line Definition lower bound
	outPutFile << "\nCoefficients corresponding to shared candidate Line Definition lower bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			if ((*candIterator)->getFlowDir() == 1) {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines), ar.push_back(-2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines << "\t" << -BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
			else {
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)), ar.push_back(1);
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines) << "\t" << 1 << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank()), ar.push_back(-1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount + (*candIterator)->getExtNodeRank() << "\t" << -1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID()), ar.push_back(1/((*candIterator)->getReactance()));
				++index;
				outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*candIterator)->getIntlNodeID() << "\t" << 1/((*candIterator)->getReactance()) << "\n";
				ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines), ar.push_back(-2.5*((*candIterator)->getFlowLimit()));
				++index;
				outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+2*sharedCLines)-scenCounter*sharedCLines << "\t" << -BIGM << "\n";
				++rCount; // Increment the row count to point to the next transmission line object
			}
		}
	}
	/*******************************************************************************************/

        // Coefficients corresponding to Internal candidate Line Forward Flow Limit Constraints
        outPutFile << "\nCoefficients corresponding to Internal candidate Line Forward Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(countOfScenarios*internalCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(-((*intCandIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*internalCLines+rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << -((*intCandIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
        	}
	}
        // Coefficients corresponding to Internal candidate Line Reverse Flow Limit Constraints\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			++index;
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(((*intCandIterator)->getFlowLimit()));
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << ((*intCandIterator)->getFlowLimit()) << "\n";
			++index;
			++rCount;
       		}
	}
	/*******************************************************************************************/

        // Coefficients corresponding to Internal candidate Line Definition upper bound
        outPutFile << "\nCoefficients corresponding to Internal candidate Line Definition upper bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1()), ar.push_back(-1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2()), ar.push_back(1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(2.5*((*intCandIterator)->getFlowLimit()));
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << BIGM << "\n";
			++rCount; // Increment the row count to point to the next transmission line object
        	}
	}
        // Coefficients corresponding to Internal candidate Line Definition lower bound
        outPutFile << "\nCoefficients corresponding to Internal candidate Line Definition lower bound\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+3*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount), ar.push_back(1);
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+3*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount << "\t" << 1 << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1()), ar.push_back(-1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID1() << "\t" << -1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2()), ar.push_back(1/((*intCandIterator)->getReactance()));
			++index;
			outPutFile << "\n" << rCount << "\t" << countOfScenarios*genNumber + scenCounter*nodeNumber+(*intCandIterator)->getTranslNodeID2() << "\t" << 1/((*intCandIterator)->getReactance()) << "\n";
			ia.push_back(rCount), ja.push_back(rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines), ar.push_back(-2.5*((*intCandIterator)->getFlowLimit()));
			++index;
			outPutFile << "\n" << rCount << "\t" << rCount-countOfScenarios*(genNumber+2*tranNumber+2*sharedELines+3*sharedCLines+2*internalCLines)+sharedCLines+countOfScenarios*otherNodeCount-scenCounter*internalCLines << "\t" << -BIGM << "\n";
			++rCount; // Increment the row count to point to the next transmission line object
        	}
	}
	/*******************************************************************************************/

	outPutFile << "\nCoefficient Matrix specified" << endl;
	clock_t end1 = clock(); // stop the timer
	double elapsed_secs1 = double(end1 - begin) / CLOCKS_PER_SEC; // Calculate the time required to populate the constraint matrix and objective coefficients
	outPutFile << "\nTotal time taken to define the rows, columns, objective and populate the coefficient matrix = " << elapsed_secs1 << " s " << endl;
	/* RUN THE OPTIMIZATION SIMULATION ALGORITHM */
	cout << "\nSimulation in Progress. Wait !!! ....." << endl;
	glp_load_matrix(milp, --index, &ia[0], &ja[0], &ar[0]); // Loads the Problem
	int lpMethodChoice; // Choice between Simplex or Interior Point Methods to solve the LP relaxation problem
	lpMethodChoice = lpSolveAlgo;
	lpMethodChoice = 1;
	glp_init_iocp(ipControlParam); // Initialize the Mixed Integer Programming Control Parameters with the default setting values
	ipControlParam->mip_gap = 1e-1; // Adjust the tolerance of the Integer Oprtimization part for faster convergence 
	cout << "\nTolerence for MIP is " << ipControlParam->mip_gap;
	switch (lpMethodChoice){
	case 1:
		glp_simplex(milp, NULL);
		break;
	case 2:
		glp_interior(milp, NULL);
		break;
	default:
		cout << "\nInvalid Choice" << endl;
		break;
	}
	string outLPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPK/OutLPLogGLPK" + to_string(zonalIndex) + ".txt";
	string outMILPLogFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputSummaryGLPK/OutMIPLogGLPK" + to_string(zonalIndex) + ".txt";
	glp_print_sol(milp, outLPLogFileName.c_str()); // write the solution log to the relaxed LP problem
	glp_intopt(milp, ipControlParam); // Calls the GLPK Branch & Bound-Cutting Plane and Simplex Method to solve the integer optimization problem
	glp_print_sol(milp, outMILPLogFileName.c_str()); // write the solution log to the integer optimization problem
	int stat = glp_mip_status(milp); // Outputs the solution status of the problem 

	/* DISPLAY THE SOLUTION DETAILS */
	switch (stat){
	case 1:
		outPutFile << "\nThe solution to the problem is UNDEFINED." << endl;
		cout << "\nThe solution to the problem is UNDEFINED." << endl;
		break;
	case 2:
		outPutFile << "\nThe solution to the problem is FEASIBLE." << endl;
		cout << "\nThe solution to the problem is FEASIBLE." << endl;
		break;
	case 3:
		outPutFile << "\nThe solution to the problem is INFEASIBLE." << endl;
		cout << "\nThe solution to the problem is INFEASIBLE." << endl;
		break;
	case 4:
		outPutFile << "\nNO FEASIBLE solution to the problem exists." << endl;
		cout << "\nNO FEASIBLE solution to the problem exists." << endl;
		break;
	case 5:
		outPutFile << "\nThe solution to the problem is OPTIMAL." << endl;
		cout << "\nThe solution to the problem is OPTIMAL." << endl;
		break;
	case 6:
		outPutFile << "\nThe solution to the problem is UNBOUNDED." << endl;
		cout << "\nThe solution to the problem is UNBOUNDED." << endl;
		break;
	default:
		outPutFile << "\nERROR in Solution Status." << endl;
		cout << "\nERROR in Solution Status." << endl;
		break;
	}

	//Get the Optimal Objective Value results//
	z = glp_mip_obj_val(milp);
	for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
		z += (*genIterator)->getNLCost();
	}

	// Open separate output files for writing results of different variables
	string outIntAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGLPK/internalAngleGLPK" + to_string(zonalIndex) + ".txt";
	string outCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPK/candFlowMWGLPK" + to_string(zonalIndex) + ".txt";
	string outCandDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/candidateLinesGLPK/candLineDecisionGLPK" + to_string(zonalIndex) + ".txt";
	string outExtAngFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/outputAnglesGLPK/externalAngleGLPK" + to_string(zonalIndex) + ".txt";
        string outintCandFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/intcandLinesGLPK/intcandFlowMWGLPK" + to_string(zonalIndex) + ".txt";
        string outintCandLineDecFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/intcandLinesGLPK/intcandLineDecisionGLPK" + to_string(zonalIndex) + ".txt";
        string outTranFlowFileName = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/tranLinesGLPK/tranLineFlowGLPK" + to_string(zonalIndex) + ".txt";
        string outSEFlowFileName1 = "/home/samie/code/Horizontal_Coordination_MILP_Dist_Mechanism_Design/Horizontal_Proper/output/SELinesGLPK/SELineFlowGLPK" + to_string(zonalIndex) + ".txt";
	ofstream internalAngleOut(outIntAngFileName, ios::out); //switchStateOut
	ofstream candFlowMWOut(outCandFlowFileName, ios::out); //switchOnOut
	ofstream candLineDecisionOut(outCandDecFileName, ios::out); //switchOffOut
	ofstream externalAngleOut(outExtAngFileName, ios::out);
        ofstream intCandFlowMWOut(outintCandFlowFileName, ios::out);
        ofstream intCandLineDecisionOut(outintCandLineDecFileName, ios::out);
        ofstream tranFlowOut(outTranFlowFileName, ios::out);
        ofstream SEFlowOut1(outSEFlowFileName1, ios::out);
	outPutFile << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	powerGenOut << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	cout << "\nThe Optimal Objective value (Generation Dispatch  and Line Building Decision cost) is: " << z << endl;
	x.push_back(0); // Initialize the decision Variable vector

	//Display Power Generation
	powerGenOut << "\n****************** GENERATORS' POWER GENERATION LEVELS (MW) *********************" << endl;
	powerGenOut << "GENERATOR ID" << "\t" << "GENERATOR MW" << "\n";
	int arrayInd = 1;
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (genIterator = genObject.begin(); genIterator != genObject.end(); ++genIterator){
			x.push_back(glp_mip_col_val(milp, arrayInd));
			powerGenOut << (*genIterator)->getGenID() << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	powerGenOut << "Finished writing Power Generation" << endl;

	// Display Internal node voltage phase angle variables
	internalAngleOut << "\n****************** INTERNAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	internalAngleOut << "NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (nodeIterator = nodeObject.begin(); nodeIterator != nodeObject.end(); ++nodeIterator){
			if ( (((*nodeIterator)->getSharedFlag()) == 1) || (((*nodeIterator)->getCandFlag()) == 1) ) {
				x.push_back(glp_mip_col_val(milp, arrayInd));
				internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << glp_mip_col_val(milp, arrayInd) << endl;
				coordInstanceRef.populateAngleDec(glp_mip_col_val(milp, arrayInd), (zonalIndex-1), scenCounter, ((*nodeIterator)->getGlobalRank())); // Passing on the shared node angle decision message to the MO		
				++arrayInd;
			}
			else {
				x.push_back(glp_mip_col_val(milp, arrayInd));
				internalAngleOut << (*nodeIterator)->getNodeID() << "\t" << glp_mip_col_val(milp, arrayInd) << endl;		
				++arrayInd;	
			}		
		}
	}
	internalAngleOut << "Finished writing Internal Node Voltage Phase Angles" << endl;

	// Display Shared Candidate lines' Power Flow variables
	candFlowMWOut << "\n****************** SHARED CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
	candFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
			x.push_back(glp_mip_col_val(milp, arrayInd));
			candFlowMWOut << (*candIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
			++arrayInd;
		}
	}
	candFlowMWOut << "Finished writing Shared Candidate lines' Power Flow variables" << endl;

	// Display Shared Candidate lines' Construction Decisions
	candLineDecisionOut << "\n****************** SHARED CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
	candLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
	for (candIterator = candLineObject.begin(); candIterator != candLineObject.end(); ++candIterator){
		x.push_back(glp_mip_col_val(milp, arrayInd));
		candLineDecisionOut << (*candIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
		coordInstanceRef.populateLineDec(glp_mip_col_val(milp, arrayInd), (zonalIndex-1), ((*candIterator)->getGlobalRank())); // Passing on the line building decision message to the MO
		++arrayInd;
	}
	candLineDecisionOut << "Finished writing Shared Candidate lines' Construction decisions" << endl;

	// Display Outer Zonal node angles
	externalAngleOut << "\n****************** OUTER ZONAL NODE VOLTAGE PHASE ANGLE VALUES *********************" << endl;
	externalAngleOut << "EXTERNAL NODE ID" << "\t" << "VOLTAGE PHASE ANGLE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
		int diffNodeC = 0;
		for (globalIterator = globalRankDiffNode.begin(); globalIterator != globalRankDiffNode.end(); ++globalIterator){
			if (diffNodeC > 0) { // Skip the dummy element "0" at the beginning
				x.push_back(glp_mip_col_val(milp, arrayInd));
				externalAngleOut << (*globalIterator) << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
				coordInstanceRef.populateAngleDec(glp_mip_col_val(milp, arrayInd), (zonalIndex-1), scenCounter, (*globalIterator)); // Passing on the shared node angle decision message to the MO
				++arrayInd;
			}
			++diffNodeC;
		}
	}
	externalAngleOut << "Finished writing outer zonal node voltage phase angle values" << endl;

        // Display Internal Candidate lines' Power Flow variables
        intCandFlowMWOut << "\n****************** INTERNAL CANDIDATE LINES POWER FLOW VALUES *********************" << endl;
        intCandFlowMWOut << "CANDIDATE LINE ID" << "\t" << "MW FLOW VALUE" << "\n";
	for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
        	for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
			x.push_back(glp_mip_col_val(milp, arrayInd));
			intCandFlowMWOut << (*intCandIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd))*100 << " MW" << endl;
			++arrayInd;
        	}
	}
        intCandFlowMWOut << "Finished writing Internal Candidate lines' Power Flow variables" << endl;

        // Display Internal Candidate lines' Construction Decisions
        intCandLineDecisionOut << "\n****************** INTERNAL CANDIDATE LINES CONSTRUCTION DECISIONS *********************" << endl;
        intCandLineDecisionOut << "CANDIDATE LINE ID" << "\t" << "CONSTRUCTION DECISIONS" << "\n";
        for (intCandIterator = intCandLineObject.begin(); intCandIterator != intCandLineObject.end(); ++intCandIterator){
		x.push_back(glp_mip_col_val(milp, arrayInd));
		intCandLineDecisionOut << (*intCandIterator)->getTranslID() << "\t" << (glp_mip_col_val(milp, arrayInd)) << endl;
                ++arrayInd;
        }
        intCandLineDecisionOut << "Finished writing Internal Candidate lines' Construction decisions" << endl;
        // Display Internal Transmission lines' Flows
        tranFlowOut << "\n****************** INTERNAL TRANSMISSION LINES FLOWS *********************" << endl;
        tranFlowOut << "SCENARIO ID" << "\t" << "TRANSMISSION LINE ID" << "\t" << "MW FLOW" << "\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
                for (tranIterator = translObject.begin(); tranIterator != translObject.end(); ++tranIterator){
                        tranFlowOut << scenCounter << "\t" << (*tranIterator)->getTranslID() << "\t" << (1/((*tranIterator)->getReactance()))*((glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID1())))-(glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*tranIterator)->getTranslNodeID2()))))*100 << " MW" << endl;
                }
        }
        tranFlowOut << "Finished writing Internal Transmission lines' MW Flows" << endl;

        // Display SE lines' Flows
        SEFlowOut1 << "\n****************** SHARED EXISTING TRANSMISSION LINES FLOWS *********************" << endl;
        SEFlowOut1 << "SCENARIO ID" << "\t" << "SHARED EXISTING TRANSMISSION LINE ID" << "\t" << "LINE REACTANCE" << "\t" << "FROM ANGLE" << "\t" << "TO ANGLE" << "\t" << "MW FLOW" <<"\n";
        for (int scenCounter=0; scenCounter < countOfScenarios; ++scenCounter) {
                for (exsharedIterator = SELineObject.begin(); exsharedIterator != SELineObject.end(); ++exsharedIterator){
                        if ((*exsharedIterator)->getFlowDir() == 1) {
                                SEFlowOut1 << scenCounter << "\t" << (*exsharedIterator)->getTranslID() << "\t" << (1/((*exsharedIterator)->getReactance())) << "\t" << (glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID()))) << "\t" << (glp_mip_col_val(milp, (countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank()))) << "\t" << (1/((*exsharedIterator)->getReactance()))*((glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID())))-(glp_mip_col_val(milp, (countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank()))))*100 << " MW" << endl;
                        }
                        else {
                                SEFlowOut1 << scenCounter << "\t" << (*exsharedIterator)->getTranslID() << "\t" << (1/((*exsharedIterator)->getReactance())) << "\t" << (glp_mip_col_val(milp, (countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank()))) << "\t" << (glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID()))) << "\t" << (1/((*exsharedIterator)->getReactance()))*((glp_mip_col_val(milp, (countOfScenarios*(genNumber+nodeNumber+sharedCLines)+sharedCLines+scenCounter*otherNodeCount+(*exsharedIterator)->getExtNodeRank())))-(glp_mip_col_val(milp, (countOfScenarios*genNumber + scenCounter*nodeNumber+(*exsharedIterator)->getIntlNodeID()))))*100 << " MW" << endl;
                        }
                }
        }
        SEFlowOut1 << "Finished writing Shared Existing Transmission lines' MW Flows" << endl;
	delete ipControlParam; // free the memory of the Integer Programming Control Parameter struct
	glp_delete_prob(milp); // Free the memory of the GLPK Problem Object
	clock_t end2 = clock(); // stop the timer
	double elapsed_secs2 = double(end2 - begin) / CLOCKS_PER_SEC; // Calculate the Total Time
	outPutFile << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;
	cout << "\nTotal time taken to solve the MILP Line Construction Decision Making Problem instance and retrieve the results = " << elapsed_secs2 << " s " << endl;

	// Close the different output files
	outPutFile.close();
	powerGenOut.close();
	internalAngleOut.close();
	candFlowMWOut.close();
	candLineDecisionOut.close();
	externalAngleOut.close();
        intCandFlowMWOut.close();
        intCandLineDecisionOut.close();
        tranFlowOut.close();
        SEFlowOut1.close();
	cout << "\nSimulation Completed.\nResults written on the different output files" << endl; 
	#C++ ends
	return z;
end
#Sample Julia/JuMP code starts
datadir = joinpath("ieee_test_cases") 
gens = CSV.read(joinpath(datadir,"Gen14.csv"), DataFrame);
lines = CSV.read(joinpath(datadir,"Tran14_b.csv"), DataFrame);
loads = CSV.read(joinpath(datadir,"Load14.csv"), DataFrame);

# Rename all columns to lowercase (by convention)
for f in [gens, lines, loads]
    rename!(f,lowercase.(names(f)))
end

# create generator ids 
gens.id = 1:nrow(gens);

# create line ids 
lines.id = 1:nrow(lines);
# add set of rows for reverse direction with same parameters
lines2 = copy(lines)
lines2.f = lines2.fromnode
lines2.fromnode = lines.tonode
lines2.tonode = lines2.f
lines2 = lines2[:,names(lines)]
append!(lines,lines2)

# calculate simple susceptance, ignoring resistance as earlier 
lines.b = 1 ./ lines.reactance

# keep only a single time period
loads = loads[:,["connnode","interval-1_load"]]
rename!(loads,"interval-1_load" => "demand");

lines
#=
Function to solve DC OPF problem using IEEE test cases
Inputs:
    gen_info -- dataframe with generator info
    line_info -- dataframe with transmission lines info
    loads  -- dataframe with load info
=#
function dcopf_ieee(gens, lines, loads)
    DCOPF = Model(GLPK.Optimizer) # You could use Clp as well, with Clp.Optimizer
    
    # Define sets based on data
      # Set of generator buses
    G = gens.connnode
    
      # Set of all nodes
    N = sort(union(unique(lines.fromnode), 
            unique(lines.tonode)))
    
      # sets J_i and G_i will be described using dataframe indexing below

    # Define per unit base units for the system 
    # used to convert from per unit values to standard unit
    # values (e.g. p.u. power flows to MW/MVA)
    baseMVA = 100 # base MVA is 100 MVA for this system
    
    # Decision variables   
    @variables(DCOPF, begin
        GEN[N]  >= 0     # generation        
        # Note: we assume Pmin = 0 for all resources for simplicty here
        THETA[N]         # voltage phase angle of bus
        FLOW[N,N]        # flows between all pairs of nodes
    end)
    
    # Create slack bus with reference angle = 0; use bus 1 with generator
    fix(THETA[1],0)
                
    # Objective function
    @objective(DCOPF, Min, 
        sum( gens[g,:c1] * GEN[g] for g in G)
    )
    
    # Supply demand balances
    @constraint(DCOPF, cBalance[i in N], 
        sum(GEN[g] for g in gens[gens.connnode .== i,:connnode]) 
            + sum(load for load in loads[loads.connnode .== i,:demand]) 
        == sum(FLOW[i,j] for j in lines[lines.fromnode .== i,:tonode])
    )

    # Max generation constraint
    @constraint(DCOPF, cMaxGen[g in G],
                    GEN[g] <= gens[g,:pgmax])

    # Flow constraints on each branch; 
    # In DCOPF, line flow is a function of voltage angles
       # Create an array of references to the line constraints, 
       # which we "fill" below in loop
    cLineFlows = JuMP.Containers.DenseAxisArray{Any}(undef, 1:nrow(lines)) 
    for l in 1:nrow(lines)
        cLineFlows[l] = @constraint(DCOPF, 
            FLOW[lines[l,:fromnode],lines[l,:tonode]] == 
            baseMVA * lines[l,:b] * 
            (THETA[lines[l,:fromnode]] - THETA[lines[l,:tonode]])
        )
    end
    
    # Max line flow limits
       # Create an array of references to the line constraints, 
       # which we "fill" below in loop
    cLineLimits = JuMP.Containers.DenseAxisArray{Any}(undef, 1:nrow(lines)) 
    for l in 1:nrow(lines)
        cLineLimits[l] = @constraint(DCOPF,
            FLOW[lines[l,:fromnode],lines[l,:tonode]] <=
            lines[l,:capacity]
        ) 
    end

    # Solve statement (! indicates runs in place)
    optimize!(DCOPF)

    # Output variables
    generation = DataFrame(
        node = gens.connnode,
        gen = value.(GEN).data[gens.connnode]
        )
    
    angles = value.(THETA).data
    
    flows = DataFrame(
        fbus = lines.fromnode,
        tbus = lines.tonode,
        flow = baseMVA * lines.b .* (angles[lines.fromnode] .- 
                        angles[lines.tonode]))
    
    # We output the marginal values of the demand constraints, 
    # which will in fact be the prices to deliver power at a given bus.
    prices = DataFrame(
        node = N,
        value = dual.(cBalance).data)
    
    # Return the solution and objective as named tuple
    return (
        generation = generation, 
        angles,
        flows,
        prices,
        cost = objective_value(DCOPF),
        status = termination_status(DCOPF)
    )
end
#Sample Julia/JuMP code ends
