//cell.h
//===================
// Inlcude Guards
//
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_
//===================
// forward declarations
class Tissue;
//===================
// include dependencies
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "externs.h"
#include <boost/random.hpp>
//===================
// Cell Class Declaration

class Cell: public enable_shared_from_this<Cell> {

	private:
		Tissue* my_tissue;
		int rank;
		int layer;
		int lineage;
		int boundary;
		int stem;
		double damping;
		int life_length;
		int num_cyt_nodes;
		vector<shared_ptr<Cyt_Node>> cyt_nodes;
		int num_wall_nodes;
		vector<shared_ptr<Wall_Node>> wall_nodes;
		double Cell_Progress;
		Coord cell_center;
		double cytokinin;
		double wuschel;
		int growth_rate;
		double init_Num_Nodes;
		double perimeter;
		bool growing_this_cycle;
		Coord growth_direction;
		vector<shared_ptr<Cell>> neigh_cells;
		vector<shared_ptr<Cell>> adh_neighbors;
		shared_ptr<Wall_Node> left_Corner;
		bool recent_div;
		//0 is interphase, 1 is Mother, 2 is Daughter
		bool terminal; //When all cells are marked as terminal, simulation stops.
		int recent_div_MD;
	public:

		Cell(Tissue* tissue);
		Cell(int rank, Coord center, double radius, Tissue* tiss, int layer, int boundary, int stem);
		void make_nodes(double radius);	
		// Destructor
		~Cell();

		// Getters and Setters
		//set tissue
		Tissue* get_Tissue() {return my_tissue;}
		//set/get rank
		void set_Rank(const int id);
		int get_Rank() {return rank;}
		bool get_recent_div() {return recent_div;}
		//set/get layer
		void set_Layer(int layer);
		int get_Layer() {return layer;}
		//set/get damping
		void set_Damping(double new_damping);
		double get_Damping() {return damping;}
		// get boundary
		int get_Boundary() {return boundary;}
		int get_Num_Boundary_Nodes();
		//set/get life length
		void update_Life_Length();
		void rescale_Life_Length(int old_growth_rate,bool init_phase);
		int get_life_length() {return life_length;} 
		void reset_Life_Length();
		//get total node count
		int get_Node_Count();
		void set_perimeter(double new_perimeter);
		double get_perimeter(){return perimeter;}
		double get_curr_perimeter();
		//get counts of individual nodes
		int get_Wall_count() {return num_wall_nodes;}
		int get_cyt_count() {return num_cyt_nodes;}
		//get wall nodes
		void get_Wall_Nodes_Vec(vector<shared_ptr<Wall_Node>>& walls);
		//add new wall node
		void add_Wall_Node_Vec(shared_ptr<Wall_Node> curr);
		//refresh list of wall nodes - ensures no fake nodes
		void refresh_Wall_Nodes_Vec();
		//get cyt nodes
		void get_Cyt_Nodes_Vec(vector<shared_ptr<Cyt_Node>>& cyts);
		//add new cyt node
		void update_cyt_node_vec(shared_ptr<Cyt_Node> new_node);
		//reset cell_Progress
		void reset_Cell_Progress();
		//set/get cell progress
		void update_Cell_Progress(int& Ti);
		double get_Cell_Progress() {return Cell_Progress;}
		double get_Init_Num_Nodes() {return init_Num_Nodes;}
		void set_Init_Num_Nodes(double inn);
		//Identify boundaries
		void identify_Boundaries(); 
		//get cell center
		void update_Cell_Center();
		Coord get_Cell_Center() {return cell_center;}
		//set/get WUS conc
		void calc_WUS(Coord L1_AVG, double WUS_dropdown);
		double get_WUS_concentration() {return wuschel;}
		void calc_CK(Coord L1_AVG, double CK_dropdown);
		double get_CYT_concentration() {return cytokinin;}
		double getRandomDoubleUsingNormalDistribution(double mean, double sigma);
		//set growth rate based on WUS
		void set_growth_rate(bool first_growth_rate);
		int get_growth_rate(){return growth_rate;}
		//set/get growth direction
		void update_growth_direction();
		void update_node_parameters_for_growth_direction();
		void set_growth_direction(Coord gd);
		Coord get_growth_direction(){return growth_direction;}
		double hill_Prob();
		//get current neighbor cells		
		void get_Neighbor_Cells(vector<shared_ptr<Cell>>& cells);
		//set/get growing_this_cycle
		void set_Growing_This_Cycle(bool gtc);
		bool is_Growing_This_Cycle(){return growing_this_cycle;}
		double calc_My_Oop_Prob();
		void get_ADH_Neighbors_Vec(vector<shared_ptr<Cell>>& adh_neighbs);
		//set/get left corner
		void set_Left_Corner(shared_ptr<Wall_Node> new_left_corner);
		shared_ptr<Wall_Node> get_Left_Corner() {return left_Corner;}			      //is this necessary?
		void set_Wall_Count(int number_nodes);
		//set/get cell terminal state
		void set_Terminal(bool t);
		bool is_Terminal(){return terminal;}

		//Computations
		double compute_membr_thresh(shared_ptr<Wall_Node> current);
		double compute_k_lin(shared_ptr<Wall_Node> current);
		double compute_k_bend(shared_ptr<Wall_Node> current);
		double compute_k_bend_div(shared_ptr<Wall_Node> current);
		void update_Wall_Angles();
		void update_Wall_Equi_Angles();
		void update_Wall_Equi_Angles_Div();
		//this function is not in use
		void update_Linear_Bending_Springs();	

		//Keep track of neighbor cells
		void update_Neighbor_Cells();
		void update_Neighbor_Cells(vector<shared_ptr<Cell>>& cell_vec, shared_ptr<Cell>);

		//adhesion
		void clear_adhesion_vectors();
		void update_adhesion_springs();
		void update_Adh_Neighbors();

		//Forces and Positionsing
		void calc_New_Forces(int Ti);
		void update_Node_Locations(int Ti);

		//Growth of a cell
		//void update_Cell_Progress(int& Ti);
		double calc_Cell_Maturity(bool cross_section_check);
		void division_check();
		double calc_Area();
		void add_Wall_Node_Check(int Ti);
		void delete_Wall_Node_Check(int Ti);
		void refresh_Walls();
		void add_Wall_Node(int Ti);
		void delete_Wall_Node(int Ti);
		void delete_Specific_Wall_Node(int Ti, shared_ptr<Wall_Node> wall);

		void find_Smallest_Length(shared_ptr<Wall_Node>& right);
		void find_Largest_Length(shared_ptr<Wall_Node>& node);//vector<pair<shared_ptr<Wall_Node>,double>>& nodes);
		Coord compute_direction_of_highest_tensile_stress();
		Coord compute_point_on_line(double t);;
		void add_Cyt_Node();

		//Functions for Division
		void find_nodes_for_div_plane(Coord& orientation, vector<shared_ptr<Wall_Node>>& nodes, int search_amount);
		void find_nodes_for_div_plane_mechanical(vector<shared_ptr<Wall_Node>>& nodes);
		void Errera_div(vector<shared_ptr<Wall_Node>>& nodes);
		void move_start_end_points(shared_ptr<Wall_Node> first, shared_ptr<Wall_Node> second, vector<shared_ptr<Wall_Node>>& daughter_ends);
		void move_cyt_nodes(Coord center_pt);
		Coord produce_random_vec();	
		Coord find_orientation_circ();
		//Output Functions
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices(int& id,bool cytoplasm);
		int num_Neighbors() {return adh_neighbors.size();} 
		void print_VTK_Adh(ofstream& ofs);
		Coord average_coordinates();
		void print_direction_vec(ofstream& ofs);

		void print_locations(ofstream& ofs,bool cytoplasm, int Ti);
		/*void print_VTK_Points(ofstream& ofs, int& count);
		  void print_VTK_Scalars_Wall_Pressure(ofstream& ofs);
		  void print_VTK_Scalars_Average_Pressure(ofstream& ofs);
		  void print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs);
		  void print_VTK_Scalars_WUS(ofstream& ofs);
		  void print_VTK_Scalars_WUS_cell(ofstream& ofs);
		  void print_VTK_Scalars_CK(ofstream& ofs);
		  void print_VTK_Scalars_Total(ofstream& ofs);
		  void print_VTK_Vectors(ofstream& ofs);
		  void print_VTK_Scalars_Node(ofstream& ofs);	
		  void print_VTK_Tensile_Stress(ofstream& ofs);
		  void print_VTK_Shear_Stress(ofstream& ofs);*/
		void print_locations(ofstream& ofs);
		void print_Cell_Data(ofstream& ofs, int Ti);
		vector<double> calc_Orientation_Stats();
		double calc_Long_Length();
		double calc_Depth();
		int get_Lineage() {return lineage;} 
		void set_Lineage(int parent_lineage);
		void print_VTK_Points(ofstream& ofs, int& count, bool cytoplasm);
		void print_VTK_Scalars_Wall_Pressure(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_Average_Pressure(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_WUS(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_WUS_cell(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_CK(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_Total(ofstream& ofs, bool cytoplasm);
		void print_VTK_Vectors(ofstream& ofs, bool cytoplasm);
		void print_VTK_Scalars_Node(ofstream& ofs, bool cytoplasm);	
		void print_VTK_Tensile_Stress(ofstream& ofs, bool cytoplasm);
		void print_VTK_Shear_Stress(ofstream& ofs, bool cytoplasm);
		void print_VTK_Cell_Progress(ofstream& ofs, bool cytoplasm);
		void print_VTK_Neighbors(ofstream& ofs, bool cytoplasm); 
		void print_VTK_Curved(ofstream& ofs, bool cytoplasm); 
		//void print_VTK_Corners(ofstream& ofs, bool cytoplasm);
		void print_VTK_Boundary(ofstream& ofs, bool cytoplasm); 
		void print_VTK_Growth_Dir(ofstream& ofs, bool cytoplasm); 
		void print_VTK_MD(ofstream& ofs, bool cytoplasm);
		void print_VTK_OOP(ofstream& ofs, bool cytoplasm);
		void print_VTK_Lineage(ofstream& ofs, bool cytoplasm);

		vector<pair<double,shared_ptr<Wall_Node>>> get_Angle_Wall_Sorted();
		//vector<shared_ptr<Wall_Node>> get_Corner_Nodes();
		void set_MD(int n){ recent_div_MD = n; return;}

		//Division 
		shared_ptr<Cell> division();
		//Debugging
		void one_To_One_Check();
		//void ensure_Only_Real_Walls();
		void NAN_CATCH(int Ti);
};


// End Cell Class
//===================

#endif

