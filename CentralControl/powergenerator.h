#ifndef POWERGENERATOR_H
#define POWERGENERATOR_H
#include <vector>
using namespace std;
class Node; // forward declaration
// Powergenerator class definition
class Powergenerator
{
private:
	double incrementalCost; // Average Incremental Cost
	double noLoadCost; // no-load cost
	double pMax, pMin; //Maximum and minimum generating capability respectively
	int genID; // Generator object id number
	Node *connNodegPtr; // connection node object

public:
	Powergenerator(int, Node*, double, double, double, double); // constructor for both average HR & piecewise linear case
	~Powergenerator(); // destructor
	void setGenParamsSimple(double, double); // set function to assign values and check validity range of parameter values for average heat rate case
	int getGenID(); // returns the ID number of the generator object
	double getLinCoeff(); // returns the linear coefficient
	double getNLCost(); // returns the no load cost
	double getPMax(); // returns the maximum power output for the generator
	double getPMin(); // returns the minimum power output for the generator
}; // end of class definition

#endif // POWERGENERATOR_H
