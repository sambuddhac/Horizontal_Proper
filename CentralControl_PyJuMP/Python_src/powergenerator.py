#Include definition of Node class 
from Python_src.node import Node
from Python_src.node import Node

#constructor definition for the piecewise linear objective function
class Powergenerator(object):
	def __init__(self, ID, nodeConng, incCost, noLoad, max, Min):
		self.genID = ID #Initializer list to initialize data members that don't need validity check
		self.connNodegPtr = nodeConng
		self.noLoadCost = noLoad
		self.incrementalCost = incCost #Initialize this data member to zero (unused for piecewise linear objective)
		self.setGenParamsSimple(max, Min) #call the set function to perform validity check on parameter value ranges and assign the values 
		self.connNodegPtr.setgConn(self.genID) #increments the generation connection variable to node
		#end of constructor


	def setGenParamsSimple(self, max, Min): #set function to set the Powergenerator class data members min max limits
		self.pMax = max
		self.pMin = Min 
		#end of setGenParams function

	def __del__(): #destructor definition
		#log.info("\nGenerator object {}  destroyed".format(self.genID))
		#end of destructor

	def getGenID(self): #returns the Powergenerator ID number
		return self.genID
		#end of getGenID function

	def getPMax(self): #function getPMax begins
		return self.pMax
		#getPMax ends

	def getPMin(self): #function getPMin begins
		return self.pMin
		#getPMax ends

	def getLinCoeff(self): #Gets the linear coefficient (Incremental production cost
		return self.incrementalCost
		#function getLinCoeff ends

	def getNLCost(self): #Gets the no load cost
		return self.noLoadCost
		#function getNLCost ends 

