import java.util.*;
import java.lang.Math;
import java.io.*;
import gurobi.*;
package thread_class;
public class HorMILPCent {

	public static void main(String[] args) // Main method begins program execution
	{
		int curveChoice = 0; // Number to indicate the type of Objective function among average heat rate, piecewise linear, or polynomial
		int systemChoice = 0; // Choice of the system to be simulated
		System.out.println("\nChoose the type of System to be simulated: 1 for Simple two bus/two region, 2 for system combined of IEEE 14, 30, and 5 node systems");
		Scanner systemChoiceInput = new Scanner(System.in);
		systemChoice = systemChoiceInput.nextInt();
		curveChoice = 1; //Assume Average Heat Rate for now
		// Read the master zones file, for deciding upon which other files to read for building the model
		int numberOfZones = 0; // Number of zones between which horizontal investment coordination for transmission lines to be built is considered
		int numberOfFields = 0; // Number of rows or individual file types for each of the zones
		string inputMasterFile;
		if (systemChoice == 1)
			inputMasterFile = "masterZonesSummary.txt";
		else
			inputMasterFile = "masterZonesSummaryRevised.txt";
		ifstream zoneSummaryFile (inputMasterFile, ios::in )
		; // ifstream constructor opens the master zones summary file
		stringstream buffer; // stringstream object to store the read information from the summary file
		// exit program if ifstream could not open file
		if (!zoneSummaryFile) {
			cerr << "\nMaster file for Zones Summary could not be opened\n" << endl;
			exit(1);
		} // end if

		//zoneSummaryFile >> numberOfZones >> numberOfFields; // get the number of zones and the number of fields: Future expansion
		System.out.println("\nEnter the number of zones");
		Scanner zoneNumInput = new Scanner(System.in);
		numberOfZones = zoneNumInput.nextInt(); // User input the number of zones/regions
		GRBEnv environmentGUROBI = new GRBEnv("GUROBILogFile.log"); // GUROBI Environment object for storing the different optimization models
		numberOfFields = 7; // Number of fields
		buffer << zoneSummaryFile.rdbuf(); // reads the data in the summary file
		string test = buffer.str(); // Extract the strings from the buffer to "test"

		//create variables that will act as "cursors". we'll take everything between them.
		size_t pos1 = 0;
		size_t pos2;
		//create the array to store the strings.
		string str[ numberOfFields * numberOfZones];
		//Read the summary input file
		for (int i = 0; i < numberOfFields; ++i) {
			for (int j = 0; j < numberOfZones; ++j) {
				if (j == numberOfZones - 1) {
					pos2 = test.find("\n", pos1); //search for the bar "\n". pos2 will be where the bar was found.
					str[i * numberOfZones + j] = test.substr(pos1, (pos2 - pos1)); //make a substring, wich is nothing more
					//than a copy of a fragment of the big string.
					pos1 = pos2 + 1; // sets pos1 to the next character after pos2.
				} else {
					pos2 = test.find(" ", pos1); //search for the bar " ". pos2 will be where the bar was found.
					str[i * numberOfZones + j] = test.substr(pos1, (pos2 - pos1)); //make a substring, wich is nothing more
					//than a copy of a fragment of the big string.
					pos1 = pos2 + 1; // sets pos1 to the next character after pos2.
				}
			}
		}
		System.out.println("\n*** NETWORK INITIALIZATION STAGE BEGINS ***\n");
		Nettran nettranInstance = new Nettran(str, numberOfZones, curveChoice); // create the network instances for the different zones
		System.out.println("\n*** NETWORK INITIALIZATION STAGE ENDS: ZONAL SUB-NETWORKS CREATED ***\n");
		System.out.println("\n*** SOLUTION OF SINGLE AREA MILP HORIZONTAL COORDINATION BEGINS ***\n");
		switch (curveChoice) {
			case 1 -> {
				System.out.println("\nSOLVING MILP");
				double solution = nettranInstance.MILPAvgHRGUROBI(environmentGUROBI); // Perform unit commitment for average heat rate objective
				System.out.println("\nMILP SOLVED");
			}
			case 2 -> nettranInstance.MILPPiecewiseLin(); // Perform unit commitment for piecewise linear objective
			case 3 -> nettranInstance.MILPPolynomial(); // Perform unit commitment for polynomial objective
			default -> System.out.println("\nInvalid choice of Objective function");
		}
		System.out.println("\n*** SOLUTION OF SINGLE AREA MILP HORIZONTAL COORDINATION ENDS ***\n");
		delete nettranInstance; // Free the memory of the Nettran class object
		delete environmentGUROBI; // Free the memory of the GUROBI environment object
	} // End of main method
}

