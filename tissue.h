//tissue.h
//=================
//Include Guards
#ifndef _TISSUE_H_INCLUDED_
#define _TISSUE_H_INCLUDED_
//=========================
//forward declarations

//=======================
//Include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <random>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "externs.h"
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
//=========================
// Tissue Class Declaration

//class Tissue: public enable_shared_from_this<Tissue> {
class Tissue {
	private:
		// We'll need a better data structure for later
		vector<shared_ptr<Cell>> cells;
		int num_cells;
		int num_deleted;
		int num_IP_divs;
		int num_divs;
		int num_boundary_nodes;
		int time_frozen;
		double avg_cell_diam;
		bool theta_flag; //For Anchored boundary, identifies whether 
		// the simulation has reached its anchoring theta.
		Coord top_cell_center;
		Coord left_boundary_dir;
		Coord right_boundary_dir;
		mt19937 gen;
		vector<int> dist1;
		vector<int> dist2;
		vector<int> dist3;
		vector<int> dist4;
		vector<int> counts;	
		//vector of theta values of div planes,
		//angle compared to (1,0) (positive horiz. unit vector).
		vector<double> divplanes;
		vector<int> divplanes_layers;
	public:
		Tissue(string filename, mt19937 gen);
		void get_Cells(vector<shared_ptr<Cell>>& cells);
		//set/get the number of cells in the tissue
		void update_Num_Cells(shared_ptr<Cell>& new_Cell);
		int  get_num_cells() {return num_cells;}
		Coord Compute_L1_AVG();
		void set_Theta_Flag(bool tf);
		void update_Theta_Flag();
		bool get_Theta_Flag() {return theta_flag;}
		void update_Avg_Cell_Diameter();
		void set_Time_Frozen(int Ti);
		int get_Time_Frozen() {return time_frozen;}		
		double get_Avg_Cell_Diameter() {return avg_cell_diam;}
		int return_counts(int index);
		void set_up_counts();
		void refresh_Wall_Nodes();
		void identify_Boundaries();
		void set_counts(int index);
		void assign_dist_vecs(vector<int> dist1, vector<int>dist2, vector<int> dist3, vector<int> dist4);
		int get_next_random(int dist, int count);
		void update_Signal(bool initial_update);
		void update_growth_direction();
		void update_Neighbor_Cells();
		void add_Wall(int Ti);
		void delete_Wall(int Ti);
		void update_Adhesion();
		double get_normal_number(double mean, double sigma); //mt19937 get_Random_Generator(){return gen;}
		
		// Boundary Functions

		void update_Boundary_Directions();
		void update_Num_Boundary_Nodes();
		int get_Num_Boundary_Nodes(){return num_boundary_nodes;}
		Coord get_Right_Boundary_Direction(){return right_boundary_dir;}
		Coord get_Left_Boundary_Direction(){return left_boundary_dir;}
		// Not in use --------
		void update_Linear_Bending_Springs();
		// ------
		
		void update_Cell_Cycle(int Ti);
		void division_check();
		void update_Divs();
		void update_IP_Divs();
		void update_Divplane_Vector(Coord plane,int layer_of_div);
		void calc_New_Forces(int Ti);
		void update_Cell_Locations(int Ti);
		
		//For tissue halting purposes

		bool terminal_Tissue();
		//stuff for data output
		void plot_direction_vec(ofstream& ofs);
		void print_Data_Output(ofstream& ofs);

		void locations_output(ofstream& ofs,bool cytoplasm, int Ti);
		//int update_VTK_Indices();
		//void print_VTK_File(ofstream& ofs);

		//void locations_output(ofstream& ofs);
		int update_VTK_Indices(bool cytoplasm);
		void print_VTK_File(ofstream& ofs, bool cytoplasm);
		void print_VTK_Direction_File(ofstream& ofs);	
		void inc_Num_Deleted();
		Coord get_Top_Cell_Center() {return top_cell_center;}
		//Data output functions
		void cell_Data_Output(ofstream& ofs, int Ti);
		void tissue_Data_Output(ofstream& ofs, int Ti);
		vector<double> calc_Width();
		double calc_Height();
		//Random functions
		int unifRandInt(int a, int b);
		double unifRand();
		double unifRand(double a, double b);
		//Destructor
		~Tissue();
		//Debugging
		void one_To_One_Check();
		void NAN_CATCH(int Ti);
		void BAD_CATCH(int call, int Ti);
};


//===========================
//End of file

#endif
