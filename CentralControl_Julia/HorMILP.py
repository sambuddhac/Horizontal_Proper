# HorMILP.py : Defines the entry point for the centralized/merged regions/Single control area Horizontal Investment Coordination MILP Simulation application.
# Main Method for running the Horizontal Investment Coordination MILP Simulation; Parts of the code follows the code design philosophy of Nick Laws of NREL
# Authors of the code: Sambuddha Chakrabarti & Hosna Khajeh under the guidance and supervision of Dr. Mohammad Reza Hesamzadeh & Tom Nudell, Ph.D.
import julia
import os
import subprocess
import pandas as pd
import numpy as np
import json
import sys
import traceback
from Python_src.log import log
from Python_src.profiler import Profiler
import gurobipy as gp
from gurobipy import GRB
from Python_src.nettran import Nettran

profiler = Profiler()

if sys.platform in ["darwin", "linux"]:
	log.info("Using Julia executable in {}".format(str(subprocess.check_output(["which", "julia"]), 'utf-8').strip('\n')))
elif sys.platform in ["win32", "win64", "cygwin"]:
	log.info("Using Julia executable in {}".format(str(subprocess.check_output(["which", "julia"]), 'utf-8').strip('\n')))

log.info("Loading Julia...")
profiler.start()
julSol = julia.Julia()
julSol.using("Pkg")
julSol.eval('Pkg.activate(".")')
julSol.include(os.path.join("JuMP_src", "HorMILPCentralized.jl")) # definition of Gensolver class for base case scenario first interval
log.info(("Julia took {:..2f} seconds to start and include Horizontal Investment Coordination centralized models.".format(profiler.get_interval())))

def HorMILPCentral(): # Main method begins program execution
	'''Future Work
	#Choose the type of objective function
    	'''
	systemChoice = int(input("Choose the type of System to be simulated: 1 for Simple two bus/two region, 2 for system combined of IEEE 14, 30, and 5 node systems"))
	curveChoice = 1 # Number to indicate the type of Objective function among average heat rate, piecewise linear, or polynomial; Assume Average Heat Rate for now
	# Read the master zones file, for deciding upon which other files to read for building the model
	if systemChoice==1:
		zoneSummaryFile = open(os.path.join("data", "masterZonesSummary.json")) #opens the master zones summary file
	else:
		zoneSummaryFile = open(os.path.join("data", "masterZonesSummaryRevised.json")) #opens the master zones summary file
	
	matrixFirstFile = json.load(zoneSummaryFile) #opens the file

	numberOfZones = int(input("\nEnter the number of zones")) #Number of zones between which horizontal investment coordination for transmission lines to be built is considered
	log.info("\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n")
	nettranInstance = Nettran(matrixFirstFile, numberOfZones, curveChoice) #create the network instances for the different zones
	log.info("\n*** NETWORK INITIALIZATION STAGE ENDS: ZONAL SUB-NETWORKS CREATED ***\n")
	log.info("\n*** SOLUTION OF SINGLE AREA MILP HORIZONTAL COORDINATION BEGINS ***\n")
	if curveChoice == 1:
		log.info("\nSOLVING MILP")
		solution = nettranInstance.MILPAvgHR() #Perform unit commitment for average heat rate objective
		log.info("\nMILP SOLVED")
	elif curveChoice == 2:
		log.info("\nSOLVING MILPPiecewiseLin")
		nettranInstance.MILPPiecewiseLin() #Perform unit commitment for piecewise linear objective
		log.info("\nMILPPiecewiseLin SOLVED")
	elif curveChoice == 2:
		log.info("\nSOLVING MILPPolynomial")
		nettranInstance.MILPPolynomial() #Perform unit commitment for polynomial objective
		log.info("\nMILPPolynomial SOLVED")
	else:
		log.info("\nInvalid choice of Objective function")

	log.info("\n*** SOLUTION OF SINGLE AREA MILP HORIZONTAL COORDINATION ENDS ***\n")

print("\nThis is the simulation program for the centralized/merged regions/Single control area Horizontal Investment Coordination MILP Simulation application.\n")

try:
    if __name__ == '__main__': HorMILPCentral()
except:
    log.warning("Simulation FAILED !!!!")
