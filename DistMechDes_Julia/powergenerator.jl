import numpy as np # typically this is what you want, and then instantiate numpy arrays wherever needed as np.array(). However, not sure, in this particular class, if we will need one
#import powergenerator #I am not sure abot this line #You don't need this line, since in Python (unlike C++), we don't create a separate header file for the class declaration
# Include definition of Node class
from Python_src.node import node
from math import *

# constructor definition for the linear objective function
class Powergenerator(object):
	def __init__(self, ID, nodeConng, incCost, noLoad, Max, Min): #Since Python does dynamic type checking, we don't typically need to specify the type of the variables
 		self.genID=ID # Initializer list to initialize data members that don't need validity check
 		self.connNodegPtr= nodeConng
 		self.noLoadCost=noLoad
 		self.incrementalCost = incCost # Initialize this data member to zero (unused for piecewise linear objective)
		self.setGenParamsSimple(Max, Min) # call the set function to perform validity check on parameter value ranges and assign the values
		self.connNodegPtr.setgConn(self.genID) # increments the generation connection variable to node
 	# end of constructor

	####### Make sure you take care of proper indentation. I have taken care of indentation here, if functions belong in the same module or class, you need to indent all those properlyarly
 	def setGenParamsSimple(self, Max, Min): # set function to set the Powergenerator class data members min max limits
		self.pMax = Max
		self.pMin = Min
 	# end of setGenParams function
	
	def getGenID(self): # returns the Powergenerator ID number
		return self.genID
 	# end of getGenID function

	def getPMax(self): # function getPMax begins ####### change this function name to "getPMax" instead of "PowerggetPMax"
		return self.pMax 
	# getPMax ends

	def getPMin(self): # function getPMin begins
		return self.pMin
	# getPMax ends

	def getLinCoeff(self): # Gets the linear coefficient (Incremental production cost
		return self.incrementalCost
	# function getLinCoeff ends

	def getNLCost(self): # Gets the no load cost 
		return self.noLoadCost
	# function getNLCost ends 

