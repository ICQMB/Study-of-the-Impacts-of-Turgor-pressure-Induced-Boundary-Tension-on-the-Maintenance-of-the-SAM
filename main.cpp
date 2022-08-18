#ifndef _MAIN_INCLUDE_
#define _MAIN_INCLUDE_

// Main.cpp
//============================

//===========================
// Include Dependencies
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <stdio.h>
#include <memory>
#include <random>
#include <functional>
#include <chrono>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//#include "rand.h"
//==============================

using namespace std;

//FREQUENTLY CHANGED VALUES
//Flags
bool OUT_OF_PLANE_GROWTH = true; //./batchGenerator -par -OOP_flag <1 or 0>
bool WUS_LEVEL = false;
//EXPERIMENTAL PARAMTERS
bool UNIFORM_OOP_DIST = false; // Flag -Unif <1 or 0>
double OOP_PROBABILITY = 0.3; // UNIFORM ONLY - Defaults to 0.3
double MECH_DIV_PROB = 0.5;
int DIV_MECHANISM = 5; //./batchGenerator -par -div <int>
//1 - Errera, 2 - Chem, 3 - Mech, 4 - Merged, 5 - Layer-Specific
double WUS_RAD_CONTRACTION_FACTOR = 1;//./batchGenerator -par -WR <double>
double CK_RAD_CONTRACTION_FACTOR = 1; //./batchGenerator -par -CKR <double>
int PRINT_VTKS = true; // ./batchGenerator -par -PRINT <1 or 0>
int TENSILE_CALC = 4; //./batchGenerator -par TC <int> 
int NUM_STEPS_PER_FRAME = 2500;
int VTK_PER_DATA_POINT = 5;
int RECENT_DIV_NUM_FRAMES = 10;
double THETA_ABC = 0.78539816339; // Default to pi / 4
bool CHEMICAL_GD = true; //./batchGenerator -par -Chem_GD <1 or 0>
bool BOUNDARY_PULL = true; //./batchGenerator -par -BP <1 or 0>
int BOUNDARY_PULL_TYPE = 2; // ./batchGenerator -par -BPT
//0 is off;
//1 is boundary pulled until theta formed is >= THETA_ABC (Anchored Boundary Condition / ABC)
//2 is continuous pulling by magnitude BOUNDARY_FORCE_MAGNITUDE along curve (Forced Boundary Condition / FBC)
double BOUNDARY_FORCE_MAGNITUDE = 233.0091; // ./batchGenerator -par -BFM : Experimental lower bound is 139.9599, upper bound is 313.5101, default is mean.
bool L1_L2_FORCED_ANTICLINAL_DIV = true;
int Weird_WUS = 0;
//Must be declared in externs.h
//For clarity, listed as comments in phys.h

int main(int argc, char* argv[]) {

	// reads in name of folder that 
	// stores vtk files, given in run.sh
	//This is "Animate_Cyt"
	string anim_folder = argv[1];
	//Reads in the name of folder that stores VTK files for
	//visualization without cytoplasm nodes.
	//This is "Animate_No_Cyt"
	string no_cyt_folder = argv[4];
	//reads in name of folder that 
	//stores data output, given in run.sh
	//This is "Nematic"
	string locations_cyt_folder = argv[3];
	//this is a folder that holds data output
	//without cytoplasm node info
	//This is "Locations"
	string locations_no_cyt_folder = argv[2];
	//Folders to hold Cell and Tissue level data
	string cell_data_folder = argv[5];
	string tissue_data_folder = argv[6];
	//keep track of time
	for (int i = 1; i < argc; i++) { 
		if (!strcmp(argv[i], "-WR")) { 
			WUS_RAD_CONTRACTION_FACTOR = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-CKR")) { 
			CK_RAD_CONTRACTION_FACTOR = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-div")) { 
			DIV_MECHANISM = stoi(argv[i+1]);
		} else if (!strcmp(argv[i], "-TC")) { 
			TENSILE_CALC = stoi(argv[i+1]);
		} else if (!strcmp(argv[i], "-OOP_flag")) { 
			OUT_OF_PLANE_GROWTH = stoi(argv[i+1]) ? true : false;
		} else if (!strcmp(argv[i], "-Unif")) { 
			UNIFORM_OOP_DIST = stoi(argv[i+1]) ? true : false;
		} else if (!strcmp(argv[i], "-OOP_P")) { 
			OOP_PROBABILITY = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-Chem_GD")) { 
			CHEMICAL_GD = stoi(argv[i+1]) ? true : false;
		} else if (!strcmp(argv[i], "-L1_Anti")) { 
			//If this is 1 or true, then L1 mechanism is geometric
			L1_L2_FORCED_ANTICLINAL_DIV = stoi(argv[i+1]) ? true : false;
		} else if (!strcmp(argv[i], "-PRINT")) { 
			PRINT_VTKS = stoi(argv[i+1]) ? true : false;
		} else if(!strcmp(argv[i], "-WUS_loc")) {
			Weird_WUS = stoi(argv[i+1]);
		} else if(!strcmp(argv[i], "-WUS_change")) {
			WUS_LEVEL = stoi(argv[i+1]) ? true : false;
		} else if (!strcmp(argv[i], "-BP")) { 
			BOUNDARY_PULL = stoi(argv[i+1]) ? true : false;
		} else if (!strcmp(argv[i], "-BPT")) { 
			BOUNDARY_PULL_TYPE = stoi(argv[i+1]);
		} else if (!strcmp(argv[i], "-BFM")) { 
			BOUNDARY_FORCE_MAGNITUDE = stod(argv[i+1]);
		} else if (!strcmp(argv[i], "-theta")) { 
			THETA_ABC = stod(argv[i+1]);
		}
	}
	if (DIV_MECHANISM == 0) { 
		//cout << "DIV_MECHANISM not set.  Exiting..." << endl;
		exit(1);
	} else if (DIV_MECHANISM > 5 || DIV_MECHANISM < 1) { 
		//cout << "DIV_MECHANISM Input failed. Exiting..." << endl;
		exit(1);
	}
	int start = clock();	
	std::random_device rd;
	std::mt19937 gen(rd());
	srand(time(0));
	
	//distrubution for different cell cycle lengths
	/*auto normal_rand1 = bind(normal_distribution<double> (15800,2300),mt19937(seed));
	auto normal_rand2 = bind(normal_distribution<double> (19400,1800),mt19937(seed));
	auto normal_rand3 = bind(normal_distribution<double> (24800,3600), mt19937(seed));
	auto normal_rand4 = bind(normal_distribution<double> (44600,18200), mt19937(seed));
	vector<int> dist1;
	vector<int> dist2;
	vector<int> dist3;
	vector<int> dist4;
	for (unsigned int i=0;i<600;i++){
		dist1.push_back((int) normal_rand1());
		dist2.push_back((int) normal_rand2());
		dist3.push_back((int) normal_rand3());
		dist4.push_back((int) normal_rand4());
	}*/
	//.txt file that tells initial
	//cell configuration 
	//cout << "before cell file is read in" << endl;
	string init_tissue = "initial_conditions/staggered_generated.txt";
	//string init_tissue = "one_cell.txt";
	//cout << "Read in cell starter" << endl;	


	//instantiate tissue
	//new cell and node objects
	//are made in this call
	Tissue growing_Tissue(init_tissue,gen);
	//Tissue growing_Tissue_experiment(init_tissue);
	//growing_Tissue.assign_dist_vecs(dist1, dist2, dist3, dist4);
	//cout << "Finished creating Cells" << endl;
	growing_Tissue.update_Signal(true);
	growing_Tissue.update_growth_direction();
	//cout << "Signal" << endl;
	//growing_Tissue.update_growth_direction();
	//cout << "growth direction" << endl;
	//parameters for time step
	//double numSteps = 500;

	// Variable for dataoutput
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "/Plant_Cell_";
	string Filename;
	ofstream ofs_anim;
	int out = 0;

	string noCyt_Filename;
	ofstream ofs_noCyt;
	string no_cyt_initial = "/Plant_Cell_NC_";

	//int digits2;
	string Number2;
	string Locations_no_cyt;
	ofstream ofs_loc_no_cyt;
	int out2 = 0;

	//int digits3;
	string Number3;
	string Locations_cyt;
	ofstream ofs_loc_cyt;
	string locations_initial = "/Locations_";
	int out3 = 0;
	//int digits4;
	string Number4;
	string cell_data;
	ofstream cell_data_file;
	string cell_data_initial = "/cell_data_";
	int out4 = 0;
	//int digits5;
	string Number5;
	string tissue_data;
	ofstream tissue_data_file;
	string tissue_data_initial = "/tissue_data_";
	int out5 = 0;

	//(No longer use this loop, instead flag for terminal.)
	//loop for time steps
	//which matlab file tells you how many
	//seconds each time step reprsents (2.5?)
	//for(int Ti = 0; Ti*dt < numSteps; Ti++) {

	int Ti = 0;
	int terminal_timer = 0;
	bool is_terminal;
	//is_terminal = false;
	//For hard-set time limit

	if (BOUNDARY_PULL_TYPE == 1) { 
		is_terminal = false;
		growing_Tissue.update_Boundary_Directions();
		growing_Tissue.set_Theta_Flag(false);
	} else { 
		is_terminal = true;
		growing_Tissue.set_Theta_Flag(false);
	}
	cout << "Initialized with Theta_Flag: " << growing_Tissue.get_Theta_Flag() << endl; 

	int terminal_timeout = 362500; //Plant stops 40.28 hours (exactly 145 vtks) after simulation begins
	//cout << "Setup complete" << endl;

	//Delta t is approximately 0.4s
	//int terminal_timeout = 182000; //Plant stops at 20.2 hours after all cells divide
	//int terminal_timeout = 12500; // Plant stops at 1.5 hours after all cells divide
	while (terminal_timer < terminal_timeout) {
		//keep track of simulation runs
		if (!is_terminal) {
			if (Ti%1000 == 0) {
				//is_terminal = growing_Tissue.terminal_Tissue();
				growing_Tissue.update_Theta_Flag();
				is_terminal = growing_Tissue.get_Theta_Flag();
				cout << "Timestep: " << Ti << "\n TF: " << is_terminal << endl;
				if (is_terminal) { 
					growing_Tissue.set_Time_Frozen(Ti);
				}
			}
		} else { 
			terminal_timer++;
		}

		// Tissue Growth

		//fills vector of neighbor cells for each cell
		//in tissue class this goes through each cell and calls
		//updated neighbor on each cell
		//growing_Tissue.BAD_CATCH(1,Ti);
		if(Ti%5000==0) {
			//cout << "Find Neighbors" << endl;
			growing_Tissue.update_Neighbor_Cells();
		}	

		//This if statement seems redundant
		
		
		if(Ti == 10000) {
			growing_Tissue.update_Signal(false);
			growing_Tissue.update_growth_direction();
		}
		if(Ti % 30000 == 0) {
			//cout << "update signal" << endl;
			growing_Tissue.update_Signal(false);
			growing_Tissue.update_growth_direction();
		}
		//adds one new cell wall node per cell everytime it is called
		//dont call it right away to give cell time to find initial configuration
		if(Ti >= 10000){
			if(Ti%1000 == 0){	
				//cout << "add new cell wall nodes if needed" << endl;
				growing_Tissue.add_Wall(Ti);
			}
		}
		//was used previously to help stability of cells
		//deletes a cell wall node if too close together
		if(Ti > 10000){
			//if(Ti%1000 == 500){}
			if (Ti % 1000 == 500) { 
				//cout << "delete wall" << endl;
				growing_Tissue.delete_Wall(Ti);
				//growing_Tissue.delete_Wall_Node_Check(Ti);
			}
		}

		//make adhesion pairs for each cell
		if(Ti < 10000) {
			if(Ti%1000 == 0) {
				//cout << "adhesion early" << endl;
				growing_Tissue.update_Adhesion();
			}
		} else {	
			if(Ti%50000 == 0) {
				//cout << "adhesion"<< endl;
				growing_Tissue.update_Adhesion();
			}
		}

		//Defines boundary nodes for pulling after the first
		//adhesion connections are formed.
		if(Ti == 0) { 
			growing_Tissue.identify_Boundaries();
			growing_Tissue.update_Num_Boundary_Nodes();
		}

		//adds internal node according to 
		//individual cell growth rate
		if (Ti >= 10000){
			//cout << "cell cycle" << endl;
			growing_Tissue.update_Cell_Cycle(Ti);
		}
		//will divide cell if time
		//cout << "divide necessary cells" << endl;
		if (Ti >= 10000) {
			growing_Tissue.division_check();
		}

		//Calculate new forces on cells and nodes
		//cout << "forces" << endl;
		growing_Tissue.calc_New_Forces(Ti);

		//Update node positions
		//cout << "locations" << endl;
		growing_Tissue.update_Cell_Locations(Ti);	

		//cout << "Finished" << endl;

		// print to dataOutput and vtk files

		if (Ti % NUM_STEPS_PER_FRAME == 0) { 
			growing_Tissue.refresh_Wall_Nodes();
		}
		if (PRINT_VTKS) { 

			if(Ti % NUM_STEPS_PER_FRAME == 0) {
				digits = ceil(log10(out + 1));
				if (digits == 1 || digits == 0) {
					Number = "0000" + to_string(out);
				}	
				else if (digits == 2) {
					Number = "000" + to_string(out);
				}	
				else if (digits == 3) {
					Number = "00" + to_string(out);
				}
				else if (digits == 4) {
					Number = "0" + to_string(out);
				}

				Filename = anim_folder + initial + Number + format;

				//cout << "Opening file for cyt" << endl;
				ofs_anim.open(Filename.c_str());
				//true is in reference to printing cytoplasm
				//THIS MUST BE RUN FIRST: Tensile stress is calculated
				//once for this, and NOT before the next to prevent
				//double-calculating.
				growing_Tissue.one_To_One_Check();
				growing_Tissue.print_VTK_File(ofs_anim,true);
				ofs_anim.close();	
				//cout << "Opening file for no cyt" << endl;

				noCyt_Filename = no_cyt_folder + no_cyt_initial + Number + format;

				ofs_noCyt.open(noCyt_Filename.c_str());
				//false is in reference to printing cytoplasm
				growing_Tissue.print_VTK_File(ofs_noCyt,false);
				ofs_noCyt.close();	

				out++;
			}
		}
		//cout << "VTK printed" << endl;
		//growing_Tissue.BAD_CATCH(11,Ti);
		/*if(Ti%1000==0) {
			digits2 = ceil(log10(out2 + 1));
			if (digits2 == 1 || digits2 == 0) {
				Number2 = "0000" + to_string(out2);
			}	
			else if (digits2 == 2) {
				Number2 = "000" + to_string(out2);
			}	
			else if (digits2 == 3) {
				Number2 = "00" + to_string(out2);
			}
			else if (digits2 == 4) {
				Number2 = "0" + to_string(out2);
			}

			nem_Filename = nem_folder+ initial2 + Number2 + format;

			ofs_nem.open(nem_Filename.c_str());
			growing_Tissue.print_VTK_Direction_File(ofs_nem);
			ofs_nem.close();	
			out2++;
		}*/
		//data output from simulations
		//for cell center etc
		/*if(Ti%10000 == 1){
		  nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
		  ofs_nem.open(nem_Filename.c_str());
		  growing_Tissue.nematic_output(ofs_nem);
		  ofs_nem.close();
		  Number2++;
		  }
		  if(Ti%1000 == 1){
		  locations_Filename = locations_folder + locations_initial + to_string(Number3) + ".txt";
		  ofs_loc.open(locations_Filename.c_str());
		  growing_Tissue.locations_output(ofs_loc);
		  ofs_loc.close();
		  Number3++;
		  }

			nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
			ofs_nem.open(nem_Filename.c_str());
			growing_Tissue.nematic_output(ofs_nem);
			ofs_nem.close();
			Number2++;
		}*/
		//locations with cyt nodes
		if(Ti % (NUM_STEPS_PER_FRAME * VTK_PER_DATA_POINT) == 0) {
			Locations_no_cyt = locations_no_cyt_folder + locations_initial + to_string(out2) + ".txt";
			ofs_loc_no_cyt.open(Locations_no_cyt.c_str());
			//cout << "Opened locations file..." << endl;
			growing_Tissue.locations_output(ofs_loc_no_cyt,false,Ti);
			ofs_loc_no_cyt.close();
			out2 += VTK_PER_DATA_POINT;
			//Locations_cyt = locations_cyt_folder + locations_initial + to_string(out3) + ".txt";
			//ofs_loc_cyt.open(Locations_cyt.c_str());
			//growing_Tissue.locations_output(ofs_loc_cyt,true,Ti);
			//ofs_loc_cyt.close();
			out3 += VTK_PER_DATA_POINT;

			//TISSUE MUST BE RUN FIRST.  It updates
			//various tissue characteristics for cell
			//data to be computed:
			//Cell Depth, (maybe more in future)
			tissue_data = tissue_data_folder + tissue_data_initial + to_string(out5) + ".txt";
			tissue_data_file.open(tissue_data.c_str());
			//cout << "Tissue output file opened..." << endl;
			growing_Tissue.tissue_Data_Output(tissue_data_file,Ti);
			tissue_data_file.close();
			out5 += VTK_PER_DATA_POINT;

			cell_data = cell_data_folder + cell_data_initial + to_string(out4) + ".txt";
			cell_data_file.open(cell_data.c_str());
			//cout << "Cell data output opened..." << endl;
			growing_Tissue.cell_Data_Output(cell_data_file,Ti);
			cell_data_file.close();
			out4 += VTK_PER_DATA_POINT;
		}
		//growing_Tissue.BAD_CATCH(12,Ti);

		Ti++;
	}
	/*nem_Filename = nematic_folder + nem_initial + to_string(Number2) + ".txt";
	  ofs_nem.open(nem_Filename.c_str());
	  growing_Tissue.nematic_output(ofs_nem);
	  ofs_nem.close();*/

	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;

	return 0;

}










#endif
