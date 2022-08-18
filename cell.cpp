//cell.cpp:
//===================
// Forward Declarations
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <random>
#include <memory>
#include "phys.h"
#include "coord.h"
#include "node.h"
//#include "rand.h"
#include "cell.h"
#include "tissue.h"
#include <boost/random.hpp>
//===================
//

// Cell Class Member functions

// Constructors
// this constructor is used in
// the divsion function
Cell::Cell(Tissue* tissue) {
	my_tissue = tissue;
	//rank assigned in division function
	//layer inherited from parent	
	//boundary cells don't divide
	boundary = 0;
	//stem cells don't divide 
	stem = 0;
	//damping assigned in div function
	//just divided so reset life length
	life_length = 0;
	recent_div = true;
	//cyt nodes reassigned in division function
	//start at zero
	num_cyt_nodes = 0;
	//wall nodes renumbered in division function
	//start at zero
	num_wall_nodes = 0;
	//cell progress decided in div function
	Cell_Progress = 0;
	//center calculate in division function
	//will calculate signals in div function
	wuschel = 0;
	cytokinin = 0;
	//growth_rate assigned in div function
	//growth direction assigned in division
	//neighbors assigned in div function
	//left corner assigned in division
	//recent_div_MD = 0 means that no cells are recently divided
	//same with recent_div = false
	recent_div = false;
	set_Terminal(true);
	recent_div_MD = 0;
}
//this constructor is used to initialize first set of cells
//calls set_growth_rate which detemrines growth rate based on WUS CONC
Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer, int boundary, int stem)    {
	this->my_tissue = tiss;
	this->rank = rank;
	this->layer = layer;
	set_Lineage(rank);
	recent_div = false;
	//Determines if this tissue is growing out of plane initially
	if (my_tissue->unifRand()  < calc_My_Oop_Prob()) { 
		set_Growing_This_Cycle(false);
	} else { 
		set_Growing_This_Cycle(true);
	}
	//if boundary is equal to one 
	//then the cell will have higher damping
	//which is assigned below
	//boundary conditions are read in from initial text file
	this->boundary = boundary;
	this->stem = stem;
	//set damping for cells that act as anchor points
	if(this->stem == 1) {
		this->damping = STEM_DAMP;
		set_Terminal(true);
	} else if((this->boundary == 1)) {
		this->damping =  BOUNDARY_DAMP;
		set_Terminal(true);
	} else {
		this->damping = REG_DAMP;
		set_Terminal(false);
	}
	//life_length = 0;
	//cyt nodes initialized in tissue constructor which
	//calls the makes nodes function on each new cell
	num_cyt_nodes = 0;
	//wall nodes initialized in tissue constructor which 
	//calls the make nodes function on each new cell
	num_wall_nodes = 0;

	set_Init_Num_Nodes(static_cast<double>(INIT_NUM_CYT_NODES));


	/*if((this->layer == 1)||(this->layer == 2)){
		Cell_Progress = my_tissue->unifRand(0.25,0.5);
	} else {
		Cell_Progress = my_tissue->unifRand(0.5,0.75);
	}*/
	Cell_Progress = my_tissue->unifRand(0.15,0.85);
	
	//cout << "CELL PROGRESS CONSTRUCTOR: " << Cell_Progress << endl;
	//Cell_Progress = my_tissue->unifRandInt(0,10);
	this->cell_center = center;
	//this gets reupdated after singal is assigned
	//in tissue constructor
	growth_direction = Coord(0,0);
	//This is a placeholder.  This is updated in main after signal calculation.
	recent_div = false;
	recent_div_MD = 0;



	//cout << "layer" << this->layer << endl;
	//cout << "stem" << this->stem << endl;
	//cout << "boundary" << this-> boundary << endl;
	//cout << "gd" << this->growth_direction << endl;
	//cout << "damping" << this->damping << endl;
	//neighbors update function called after initialization
	//left corner set in make nodes function called by tissue constuctor
}
//calls compute membrane equi length for each node
//calls compute linear spring constant for each node
//calls compute bending spring constant for each node
//calls update wall equi angles for each node
//calls update wall angles to get initial angle of each node
/*void Cell::make_nodes_experimental(string filename){
  ifstream ifs(filename.c_str());
  if(!ifs){
  cout << filename << "is not available" << endl;
  return;
  }

  stringstream ss;
  string line;
  string temp;
  char trash;
  int rank;
  double x, y;
  while (getline(ifs,line)) {
  ss.str(line);

  getline(ss,temp,',');

  if (temp == "CellRank") {
  ss >> rank;
  }
  else if (temp == "Center") {
  ss >> x >> trash >> y;
  Coord loc(x,y);
  center = loc;
  }
  else if (temp == "Radius") {
  ss >> radius;
  }
  else if (temp == "Layer") {
  ss >> layer;
  }
  else if (temp == "Boundary"){
  ss >> boundary;
  }
  else if(temp == "Stem"){
  ss >> stem;
  }
  else if (temp == "End_Cell") {
//create new cell with collected data 
//and push onto vector that holds all cells in tissue 
//cout<< "making a cell" << endl;
shared_ptr<Cell> curr= make_shared<Cell>(rank, center, radius, my_tissue, layer,boundary, stem);
//give that cell wall nodes and internal nodes
//curr->make_nodes(radius);
curr->make_nodes_experimental("experimental_nodes.txt");
//cout<< "make nodes" << endl;
num_cells++;
cells.push_back(curr);
}
ss.clear();
}
ifs.close();
}*/

void Cell::make_nodes(double radius){

	//assemble the membrane
	int num_Init_Wall_Nodes = INIT_WALL_NODES + 2*(floor(calc_Cell_Maturity(growing_this_cycle)));
	double angle_increment = (2*pi)/num_Init_Wall_Nodes;

	//make all wall nodes
	double curr_X;
	double curr_Y;
	Coord location = this->cell_center;;
	double curr_theta = 0;
	curr_X = cell_center.get_X() + radius*cos(curr_theta);
	curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
	location = Coord(curr_X,curr_Y);
	//make the first node
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Wall_Node> prevW = make_shared<Wall_Node>(location,this_cell);
	shared_ptr<Wall_Node> currW = prevW;
	wall_nodes.push_back(prevW);
	num_wall_nodes++;
	shared_ptr<Wall_Node> orig(prevW);
	//this will be the "starter" node
	this->left_Corner = orig;
	//make successive nodes
	for (int i = 0; i < num_Init_Wall_Nodes - 1; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = cell_center.get_X() + radius*cos(curr_theta);
		curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		shared_ptr<Wall_Node> new_node =make_shared<Wall_Node>(location,this_cell);
		currW = new_node;
		wall_nodes.push_back(currW);
		num_wall_nodes++;
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}

	//connect last node to starter node
	currW->set_Left_Neighbor(orig);
	orig->set_Right_Neighbor(currW);
	this->perimeter = this->get_curr_perimeter();
	//where is where most private member variables are set
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//damping inherited from cell
	double new_damping = this->get_Damping();
	//variable to hold membr_equ_len
	double l_thresh = 0;
	double k_lin = 0;
	double k_bend = 0;
	for(unsigned int i = 0; i < walls.size();i++) {	
		walls.at(i)->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(walls.at(i));
		k_lin = compute_k_lin(walls.at(i));
		k_bend = compute_k_bend(walls.at(i));
		walls.at(i)->set_membr_len(l_thresh);
		walls.at(i)->set_K_LINEAR(k_lin);
		walls.at(i)->set_K_BEND(k_bend);
	}	
	//double new_damping = this->get_Damping();
	//insert cytoplasm nodes
	//int num_init_cyt_nodes = Init_Num_Cyt_Nodes + Cell_Progress;
	int num_init_cyt_nodes = floor(
				calc_Cell_Maturity(is_Growing_This_Cycle()) 
			);
	//this->Cell_Progress = num_init_cyt_nodes;
	double scal_x_offset = 0.8;
	//Coord location;
	double x;
	double y;
	for (int i = 0; i < num_init_cyt_nodes; i++) {
		// USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		double rand_radius = (static_cast<double>(rand()) / RAND_MAX)*scal_x_offset*radius;
		double rand_angle = (static_cast<double>(rand()) / RAND_MAX)*2*pi;
		x = cell_center.get_X()+ rand_radius*cos(rand_angle);
		y = cell_center.get_Y()+ rand_radius*sin(rand_angle);
		location = Coord(x,y);

		shared_ptr<Cyt_Node> cyt = make_shared<Cyt_Node>(location,this_cell);
		cyt->set_Damping(new_damping);
		cyt_nodes.push_back(cyt);
		num_cyt_nodes++;
	}

	//update equilibrium angle
	update_Wall_Equi_Angles();
	//update wall angles
	update_Wall_Angles();	
	//is_divided = false;
	return;
}

// Destructor
Cell::~Cell() {
	//not needed using smartpointers
}
//=============================================================
//========================================
// Getters and Setters
//========================================
//=============================================================
void Cell::set_Rank(const int id) {
	this->rank = id;
	return;
}
void Cell::set_Layer(int layer) {
	this->layer = layer;
	return;
}
void Cell::set_Damping(double new_damping) {
	this->damping = new_damping;
	return;
}

int Cell::get_Num_Boundary_Nodes() {
	if(get_Boundary() != 1) { 
		return 0;
	}
	int num = 0;
	shared_ptr<Wall_Node> curr = left_Corner;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;
	do { 
		next = curr->get_Left_Neighbor();
		if (curr->is_Boundary()) {
			num++;
		}
		curr = next;
	} while(next != orig);

	return num;
}

//Life length is the number of time steps this cell has lived 
//in order to be at its current cell progress (Not actually the time since
//creation, since changing this changes growth rate.
void Cell::update_Life_Length() {
	life_length++;
	return;
}
void Cell::reset_Life_Length(){
	life_length = 0;
	return;
}
int Cell::get_Node_Count() {
	return num_wall_nodes + num_cyt_nodes;
}
void Cell::get_Wall_Nodes_Vec(vector<shared_ptr<Wall_Node>>& walls) {
	walls = this->wall_nodes;
	return;
}
void Cell::add_Wall_Node_Vec(shared_ptr<Wall_Node> curr) {
	this->wall_nodes.push_back(curr);
	this->num_wall_nodes++;
	return;
}
void Cell::get_Cyt_Nodes_Vec(vector<shared_ptr<Cyt_Node>>& cyts) {
	cyts = cyt_nodes;
	return;
}
void Cell::update_cyt_node_vec(shared_ptr<Cyt_Node> new_node){
	this->cyt_nodes.push_back(new_node);
	this->num_cyt_nodes++;
	return;
}
void Cell::reset_Cell_Progress(){
	this->Cell_Progress = 0;
	return;
}
/*
void Cell::update_Cell_Progress() {
	this->Cell_Progress++;
	return;
}*/
void Cell::calc_WUS(Coord L1_AVG, double WUS_dropdown) {
	
	//new data from eric
	//CZ ~5 cells wide
	//layer 1
	//from 2018 paper
	if (Weird_WUS == 0) {
		double distance = (cell_center-(L1_AVG-Coord(0,WUS_dropdown))).length();
		distance = distance * WUS_RAD_CONTRACTION_FACTOR; 
		this->wuschel = 84.6*exp(-0.01573*(distance));
	} else if (Weird_WUS == 1) {
		double distance = (cell_center-(L1_AVG-Coord(0,WUS_dropdown))).length();
		distance = distance * WUS_RAD_CONTRACTION_FACTOR; 
		if(layer == 3){
			this->wuschel = 84.6*exp(-0.01573*(distance));
		} else {
			this->wuschel = 65;
		}
	} else if (Weird_WUS == 2) {
		double distance = (cell_center-(L1_AVG-Coord(0,WUS_dropdown))-Coord(-40,0)).length();
		distance = distance * WUS_RAD_CONTRACTION_FACTOR; 
		if((layer == 1)||(layer == 2)||(layer ==3)){
			this->wuschel = 84.6*exp(-0.01573*(distance));
		} else {
			this->wuschel = 65;
		}
	}
	if (WUS_LEVEL) {
		double distance = (cell_center-(L1_AVG-Coord(0,WUS_dropdown))).length();
		distance = distance * WUS_RAD_CONTRACTION_FACTOR; 
		this->wuschel = 64.6*exp(-0.01573*(distance));
	}


	return;
}
void Cell::calc_CK(Coord L1_AVG, double CK_dropdown) {
	double distance = (cell_center-(L1_AVG-Coord(0,CK_dropdown))).length();
	distance = distance * CK_RAD_CONTRACTION_FACTOR;
	if((this->get_Layer() == 1)||(this->get_Layer() == 2)) {
		this->cytokinin = 0;
	} else {
		// if((this->get_Layer() >2) && (this->get_Layer() < 6))
		this->cytokinin = 110*exp(-0.01637*distance);
	}


	return;
}
void Cell::set_growth_rate(bool first_growth_rate) {
	//first quartile is 15-21 hours
	//second quartile is 21-27 hours
	//third qurtile is 27-39 hours
	//fourth quartile is 39-93 hours
	//.1 min per timestep
	//first distribution mean/sigma 10800/1800
	//second distribution mean/sigma 14400/1800
	//third distribution mean/sigma 19800/3600
	//fourth distribution mean/sigma 39600/16200
	double mean;
	double sigma;
	int old_growth_rate;
	if (!first_growth_rate) {
		old_growth_rate = growth_rate;
	}
	/*

	if(this->wuschel <55 ) {
		//cout << "Count 1: " << my_tissue->return_counts(0) << endl;
		this->growth_rate = this->my_tissue->get_next_random(1,this->my_tissue->return_counts(0));
		this->my_tissue->set_counts(0);
	} else if(this->wuschel <65 ) {

		//cout << "Count 2: " << my_tissue->return_counts(1) << endl;	
		this->growth_rate = this->my_tissue->get_next_random(2,this->my_tissue->return_counts(1));
		this->my_tissue->set_counts(1);
	} else if(this->wuschel <75 ) {

		//cout << "Count 3: " << my_tissue->return_counts(2) << endl;
		this->growth_rate = this->my_tissue->get_next_random(3,this->my_tissue->return_counts(2));
		this->my_tissue->set_counts(2);
	} else {

		//cout << "Count 4: " << my_tissue->return_counts(3) << endl;
		this->growth_rate = this->my_tissue->get_next_random(4,this->my_tissue->return_counts(3));
		this->my_tissue->set_counts(3);
	}*/

	//this->growth_rate = my_tissue->unifRandInt(5000,30000);
	//mt19937 gen = this->get_Tissue()->get_random_generator();
	if(this->wuschel < 55) {
		//WUS less than 55
		mean = 10800; // 18 hr
		sigma = 1800/2; // 1.5 hr
		this->growth_rate = getRandomDoubleUsingNormalDistribution(mean,sigma);
		//this->growth_rate = my_tissue->unifRandInt(2000,10000);
		//cout << "growth rate:" << growth_rate << endl;
	} else if(this->wuschel < 65) {
		// WUS betseen 55 and 65
		mean = 14400; // 24 hr
		sigma = 1800/2; // 1.5 hr
		this->growth_rate = getRandomDoubleUsingNormalDistribution(mean,sigma);
		//this->growth_rate = my_tissue->unifRandInt(10000,12510);
		//cout << "growth rate" << growth_rate << endl;
	} else if(this->wuschel < 75) {
		//Wus between 65 and 75
		mean = 19800; // 33hr
		sigma = 3600/2; // 6 hr
		this->growth_rate = getRandomDoubleUsingNormalDistribution(mean,sigma);
	} else {
		//If WUS > 75
		mean = 39600; //66 hr
		sigma = 16200/2; //13.5 hr
		this->growth_rate = getRandomDoubleUsingNormalDistribution(mean,sigma);
	}

	//else if((this->wuschel >=78) && (this->wuschel <81.8)){
	//	this->growth_rate = my_tissue->unifRandInt(30024,35526);
	//}
	//else{
	//	this->growth_rate = my_tissue->unifRandInt(37530,40032);
		//}
	//if((this->cytokinin >= 80)){

	//this->growth_rate = growth_rate*.7;
	//}

	// TEST
	//r
	//this->growth_rate = growth_rate*.7; 
	//h
	//this->growth_rate = growth_rate*0.35;

	//this->growth_rate = 5000;
	//2018 paper	
	/*if(this->wuschel < 11){
	  this->growth_rate = my_tissue->unifRandInt(5000,8000);
	  }
	  else if((this->wuschel >= 12) &&(this->wuschel <24)) {
	  this->growth_rate = my_tissue->unifRandInt(8000,10000);
	  }
	  else if((this->wuschel >= 24) && (this->wuschel <36)){
	  this->growth_rate = my_tissue->unifRandInt(10000,12000);
	  }
	  else if ((this->wuschel >= 36) && (this->wuschel <48)){
	  this->growth_rate = my_tissue->unifRandInt(12000,14000);
	  }
	  else if ((this->wuschel >= 48) && (this->wuschel < 60)){
	  this->growth_rate = my_tissue->unifRandInt(14000,16000);
	  }	
	  else if ((this->wuschel >= 60) && (this->wuschel <72)){
	  this->growth_rate = my_tissue->unifRandInt(16000,17000);
	  }
	  else if ((this->wuschel >= 72) && (this->wuschel < 84)){
	  this->growth_rate = my_tissue->unifRandInt(17000,18000);
	  }
	  else if ((this->wuschel >= 84) && (this->wuschel < 96)){
	  this->growth_rate = my_tissue->unifRandInt(18000,20000);
	  }
	  else if((this->wuschel >=96)&&(this->wuschel < 108)) {
	  this->growth_rate = my_tissue->unifRandInt(200000,22000);
	  }
	  else if((this->wuschel >=108)&&(this->wuschel < 120)) {
	  this->growth_rate = my_tissue->unifRandInt(22000,24000);
	  }
	  else if((this->wuschel >=120)&&(this->wuschel < 132)) {
	  this->growth_rate = my_tissue->unifRandInt(24000,27000);
	  }
	  else if(this->wuschel>= 132) {
	  this->growth_rate = my_tissue->unifRandInt(27000,30000);
	  }*/


	if(!first_growth_rate) { 
		rescale_Life_Length(old_growth_rate,false);
	} else {
		rescale_Life_Length(this->growth_rate,true);
	}

	return;

}

double Cell::getRandomDoubleUsingNormalDistribution(double mean, double sigma){
	double gr;
	gr = this->my_tissue->get_normal_number(mean,sigma);
	return gr;
}

void Cell::rescale_Life_Length(int old_growth_rate, bool init_phase) { 
	if (!init_phase) { 
		//Sanity check
		if (old_growth_rate == this->growth_rate) { 
			return;
		} 
		//Translate LL into CP from old life length
		//Cell_Progress = (double)life_length / (double)(30*old_growth_rate);
		Cell_Progress = (double)life_length / (double)(15*old_growth_rate);
	//cout << "CELL PROGRESS RLL: " << Cell_Progress << endl;
	} else { 
		//Do nothing
	}
	//Calculate new LL from CP
	this->life_length = floor(Cell_Progress * 
			static_cast<double>(15*this->growth_rate));
	//static_cast<double>(30*this->growth_rate));

	return;
}

void Cell::update_growth_direction(){
	//signaling stuff
	//Isotropic growth if you're not growing this cycle or if
	//you're in the boundary
	
	if (CHEMICAL_GD) { 
		if (HILL_PROB) { 

			if(this->stem == 1) { //Stem grows vertically
				this->growth_direction = Coord(0,1);
			} else if(this->boundary == 1 || !growing_this_cycle) {
				this->growth_direction = Coord(0,0);
			} else if(my_tissue->unifRand() <  hill_Prob()) {
				//hill_prob is probability of periclinal
				//growth as a function of CK/WUS
				this->growth_direction = Coord(0,1);
			} else { 
				this->growth_direction = Coord(1,0);
			}

		} else { 
			//Deprecated block - unused
			//
			if(this->boundary == 1 || !growing_this_cycle){
				this->growth_direction = Coord(0,0);
			} else if(this->stem == 1) { //Stem grows vertically
				this->growth_direction = Coord(0,1);
			} else if((this->layer == 1)||(this->layer == 2)) { //Top layers grow horizontally
				this->growth_direction = Coord(1,0);
			} else if(this->wuschel > this->cytokinin) { //Wus > CK means horizontal
				this->growth_direction = Coord(1,0);
			} else { //CK > Wus means vertical
				this->growth_direction = Coord(0,1);
			}
		}
	} else {  
		if (this->boundary == 1 || !growing_this_cycle) {
			this->growth_direction = Coord(0,0);
		} else if(this->stem == 1) {
			this->growth_direction = Coord(0,1);
		} else if((this->layer == 1)||(this->layer == 2)) {
			this->growth_direction = Coord(1,0);
		} else {
			if(rand()% 100 < 36) {
				this->growth_direction = Coord(1,0);
			} else {
				this->growth_direction = Coord(0,1);
			}
		}
	}
	this->update_node_parameters_for_growth_direction();
	return;
}
double Cell::hill_Prob() {
	double lambda = this->cytokinin / this->wuschel;
	if (HILL_K <= 0) { 
		cout << "Hill K must be >= 0!" << endl;
		return 0;
	}
	if (lambda <= HILL_K) { 
		double lambdaN = pow(lambda,HILL_N);
		double HILL_KN = pow(HILL_K,HILL_N);
		return lambdaN / (HILL_KN + lambdaN);
	} else { 
		return 1.0 / (1.0 + pow(HILL_K / lambda, HILL_N));
	}
}
void Cell::update_node_parameters_for_growth_direction(){
	vector<shared_ptr<Wall_Node>> walls;
	double k_bend;
	this->get_Wall_Nodes_Vec(walls);
	for(unsigned int i = 0; i < walls.size();i++) {	
		k_bend = compute_k_bend(walls.at(i));
		walls.at(i)->set_K_BEND(k_bend);
	}
	this->update_Wall_Equi_Angles();
	return;
}
void Cell::set_growth_direction(Coord gd){
	this->growth_direction = gd;
	return;
}
void Cell::get_Neighbor_Cells(vector<shared_ptr<Cell>>& cells) {
	cells = this->neigh_cells;
	return;
}
void Cell::set_Left_Corner(shared_ptr<Wall_Node> new_left_corner) {
	this->left_Corner = new_left_corner;
	return;
}
void Cell::set_Wall_Count(int number_nodes) {
	this->num_wall_nodes = number_nodes;
	return;
}
void Cell::set_Terminal(bool t) { 
	this->terminal = t;
	return;
}
double Cell::compute_membr_thresh(shared_ptr<Wall_Node> current) {
	//ran simulations changing equilibrium length to 
	//introduce a growth bias and it did not have a 
	//big effect 
	//could do more sensitivity analysis with this but 
	//for now all nodes will have same equilibrium length
	/*double l_thresh = 0;
	  double theta = 0;
	  double costheta = 0;
	  double curr_len = 0;
	  double growth_len = 0;
	  Coord curr_vec;	
	  curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	  curr_len = curr_vec.length();	
	  growth_len = 1;
	  costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	  theta = acos( min( max(costheta,-1.0), 1.0) );
	  if(this->growth_direction == Coord(0,1)){
	  if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
	  l_thresh = Membr_Equi_Len_Long;
	  }
	  else { 
	  l_thresh = Membr_Equi_Len_Short;
	  }
	  }*/
	double l_thresh = Membr_Equi_Len_Short;

	return l_thresh;
}

double Cell::compute_k_lin(shared_ptr<Wall_Node> current) {
	//had the idea that changin linear springs would create a
	//growth direction bias
	//but linear springs did not have the larges effect
	//for now all spring have the same linear spring constant
	/*double k_lin = 0;
	  double theta = 0;
	  double costheta = 0;
	  double curr_len = 0;
	  double growth_len = 0;
	  Coord curr_vec;	
	  curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
	  curr_len = curr_vec.length();	
	  growth_len = 1;
	  costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
	  theta = acos( min( max(costheta,-1.0), 1.0) );
	  k_lin = K_LINEAR_LOOSE + K_LINEAR_STIFF*(1-pow(costheta,2));*/

	double k_lin;
	
	if (this->boundary == 0) {
		k_lin = K_LINEAR_LOOSE;
	} else { 
		k_lin = K_LINEAR_STIFF;
	}

	return k_lin;
}
double Cell::compute_k_bend(shared_ptr<Wall_Node> current) {
	//coefficient of bending spring is very important
	//nodes that are parrallel to growth direction have 
	//high bending coefficient
	//nodes that are perpendicular to growth direction have 
	//low bending coefficient
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))||(growth_direction == Coord(0,0))) {
		//fine
	} else {
	//cout << "No growth direction assigned" << endl;
		exit(1);
	}
	double k_bend = 0;

	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))){
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		growth_len = 1;
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "Theta: " << theta << endl;
		if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
			//End caps
			k_bend = K_BEND_STIFF;
		} else { 
			//Walls
			k_bend = K_BEND_LOOSE;
		}
	} else {
		//cout << "elsesssseses" << endl;
		k_bend = K_BEND_UNIFORM;
	}
	//if((layer == 1)||(layer == 2)){

	//	k_bend = K_BEND_UNIFORM;
	//}
	//cout << "K bend: " << k_bend << endl;
	return k_bend;
}
double Cell::compute_k_bend_div(shared_ptr<Wall_Node> current) {
	//coefficient of bending spring is very important
	//nodes that are parrallel to growth direction have 
	//high bending coefficient
	//nodes that are perpendicular to growth direction have 
	//low bending coefficient
	//cout << "k bend div" << endl;
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0)) || (growth_direction == Coord(0,0))) {
		//fine
	} else {
	//cout << "No growth direction assigned" << endl;
		exit(1);
	}
	double k_bend = 0;
	if((growth_direction == Coord(0,1)) || (growth_direction == Coord(1,0))){
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		curr_vec = current->get_Left_Neighbor()->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		growth_len = 1;
		costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		theta = acos( min( max(costheta,-1.0), 1.0) );
		//cout << "Theta: " << theta << endl;
		if((theta < ANGLE_FIRST_QUAD_Div) || (theta > ANGLE_SECOND_QUAD_Div)){
			k_bend = K_BEND_STIFF;
		} else { 
			k_bend = K_BEND_LOOSE;
		}
	} else {
		k_bend = K_BEND_UNIFORM;
	}
	//cout << "K bend: " << k_bend << endl;
	return k_bend;
}

void Cell::update_Wall_Angles() {
	//cout << "wall angles" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);

#pragma omp parallel for schedule(static,1)
	for(unsigned int i=0; i< walls.size();i++) {
		//cout<< "updating" <<endl;
		walls.at(i)->update_Angle();
	}
	//cout << "Success" << endl;
	return;
}
void Cell::update_Wall_Equi_Angles() {
	//cout << "equi angles" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		int counter = 0;
		double new_equi_angle = 0; 
		double circle_angle  = (this->num_wall_nodes-2)*pi/(this->num_wall_nodes);
#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < walls.size();i++) {	
			if(this->growth_direction != Coord(0,0)) {
				curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
				curr_len = curr_vec.length();	
				growth_len = this->growth_direction.length();
				costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
				theta = acos( min( max(costheta,-1.0), 1.0) );
				if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
					new_equi_angle = pi;
				} else {
					counter++;
					new_equi_angle = circle_angle;

				}
			} else {
				//cout << "elseeeeee" << endl;
				new_equi_angle = circle_angle;
			}
			//if((layer == 1)||(layer == 2)){
			//	new_equi_angle = circle_angle;
			//}


			//this was an idea to make the round part of the cell
			//smaller in radius but not necessary
			//if(new_equi_angle != pi) {
			//	new_equi_angle =  (counter*2-2)*pi/(counter*2);
			//}

			walls.at(i)->update_Equi_Angle(new_equi_angle);
		}
	}
	return;
}
void Cell::update_Wall_Equi_Angles_Div() {
	//cout << "equi angles div" << endl;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{
		double theta = 0;
		double costheta = 0;
		double curr_len = 0;
		double growth_len = 0;
		Coord curr_vec;	
		//int counter = 0;
		double new_equi_angle = 0; 
		double circle_angle  = (this->num_wall_nodes-2)*pi/(this->num_wall_nodes);
#pragma omp parallel for schedule(static,1)
		for(unsigned int i = 0; i < walls.size();i++) {	
			if(this->growth_direction != Coord(0,0)){
				curr_vec = walls.at(i)->get_Left_Neighbor()->get_Location() - walls.at(i)->get_Location();
				curr_len = curr_vec.length();	
				growth_len = this->growth_direction.length();
				costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
				theta = acos( min( max(costheta,-1.0), 1.0) );
				if((theta < ANGLE_FIRST_QUAD_Div) || (theta > ANGLE_SECOND_QUAD_Div)){
					new_equi_angle = circle_angle;
				} else {
					//counter++;
					new_equi_angle = circle_angle;
				}
			} else {
				new_equi_angle = circle_angle;
			}
			//this was an idea to make the round part of the cell
			//smaller in radius but not necessary
			//if(new_equi_angle != pi) {
			//	new_equi_angle =  (counter*2-2)*pi/(counter*2);
			//}

			walls.at(i)->update_Equi_Angle(new_equi_angle);
		}
	}
	return;
}
void Cell::update_Cell_Center() {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);

	Coord total_location = Coord();
#pragma omp parallel
	{
		Coord curr_loc;
#pragma omp declare reduction(+:Coord:omp_out+=omp_in) initializer(omp_priv(omp_orig))
#pragma omp for reduction(+:total_location) schedule(static,1)
		for (unsigned int i = 0; i < walls.size(); i++) {
			curr_loc = walls.at(i)->get_Location();
			total_location += curr_loc;
		}
	}
	this->cell_center = total_location*((1.0/static_cast<double>(num_wall_nodes)));
	//cout << cell_center<<"center"<< endl;
	return;
}
void Cell::update_Linear_Bending_Springs(){
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//double new_damping = this->get_Damping();
	//double k_lin = 0;
	double k_bend = 0;
	//double l_thresh = 0;

	for(unsigned int i = 0; i < walls.size();i++) {	
		//walls.at(i)->set_Damping(new_damping);
		//walls.at(i)->set_membr_len(MembrEquLen);
		//k_lin = compute_k_lin(walls.at(i));
		k_bend = compute_k_bend(walls.at(i));
		//l_thresh = compute_k_bend(walls.at(i));
		//walls.at(i)->set_K_LINEAR(k_lin);
		walls.at(i)->set_K_BEND(k_bend);
		//walls.at(i)->set_membr_len(l_thresh);
	}
	update_Wall_Equi_Angles();
	return;
}
//=============================================================
//=========================================
// Keep Track of neighbor cells and Adhesion springs
//=========================================
//=============================================================
void Cell::update_Neighbor_Cells() {
	//clear prev vector of neigh cells
	//cout << "clear" << endl;
	neigh_cells.clear();
	//grab all cells from tissue
	//cout << "cleared" << endl;
	
	vector<shared_ptr<Cell>> curr_ADH_Neighbors;
	get_ADH_Neighbors_Vec(curr_ADH_Neighbors);
	vector<shared_ptr<Cell>> all_Cells;
	my_tissue->get_Cells(all_Cells);

	// Empty variables for holding info about other cells
	double prelim_threshold = 20;
	shared_ptr<Cell> sp_this = shared_from_this();
	//cout << "made pointer to cell" << endl;
	// iterate through all cells
#pragma omp parallel
	{
		Coord curr_Cent;
		Coord distance;
#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < all_Cells.size(); i++) {
			shared_ptr<Cell> curr = all_Cells.at(i);
			//cout << "made pointer to current neighbor" << endl;
			if (curr != sp_this) {
				curr_Cent = curr->get_Cell_Center();
				//			//cout << "got center" << endl;
				// Check if cell centers are close enough together
				distance = sp_this->cell_center - curr_Cent;
				//cout << "Distance = " << distance << endl;
				if ( distance.length() < prelim_threshold ) {
#pragma omp critical
					neigh_cells.push_back(curr);
					//cout << rank << "has neighbor" << curr->get_Rank() << endl;
				}

			}
			//else you're pointing at yourself and shouldnt do anything

		}
	}

	for (unsigned int i = 0; i < curr_ADH_Neighbors.size(); i++) { 
		bool already_neigh = false;
		for (unsigned int j = 0; j < neigh_cells.size(); j++) { 
			already_neigh |= curr_ADH_Neighbors.at(i) == neigh_cells.at(j);
			if (already_neigh) break;
		}
		if (!already_neigh) { 
			neigh_cells.push_back(curr_ADH_Neighbors.at(i));
		}
	}

	//cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;

}
void Cell::update_Neighbor_Cells(vector<shared_ptr<Cell>>& cell_vec, shared_ptr<Cell> sister_cell){
	neigh_cells.clear();
	this->neigh_cells = cell_vec;
	neigh_cells.push_back(sister_cell);
	return;
}

//Updates a vector of cells that are connected to this via adhesion.
void Cell::get_ADH_Neighbors_Vec(vector<shared_ptr<Cell>>& adh_vec){
	adh_vec = this->adh_neighbors;
	return;
}
void Cell::update_Adh_Neighbors() { 
	vector<shared_ptr<Cell>> temp;
	shared_ptr<Wall_Node> curr = left_Corner;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;
	shared_ptr<Wall_Node> partner = NULL;
	shared_ptr<Cell> partnerCell = NULL; 
	bool unlisted_neighbor;
	do { 
		next = curr->get_Left_Neighbor();
		for (unsigned int i = 0; i < curr->get_adh_vec().size(); i++) { 
			unlisted_neighbor = true;
			partner = curr->get_adh_vec().at(i);
			partnerCell = partner->get_My_Cell();
			for (unsigned int j = 0; j < temp.size(); j++) { 
				if (temp.at(j) == partnerCell) {
					unlisted_neighbor = false;
					break;
				}
			}
			if (unlisted_neighbor) { 
				temp.push_back(partnerCell);
			}
		}
		curr = next;
	} while(next != orig);

	adh_neighbors = temp;

	return;
}
//each cell wall node holds a vector of adhesion 
//connections and this function clears that for
//all cell wall nodes in the cell
void Cell::clear_adhesion_vectors() {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel
	{	

#pragma omp for schedule(static,1)	
		for(unsigned int i=0; i< walls.size();i++) {
			walls.at(i)->clear_adhesion_vec();
		}
	}
	return;
}
//for each cell wall node on current cell 
//this function searches through all the cell wall nodes on neighboring cells
//if a cell wall node on a neighboring cell is within the ADHthresh
//this function updates adhesion vector which is a private
//member variable for each cell wall node on the current cell
//this function pushes the current wall node on the neighboring cell
//onto adhesion vector
void Cell::update_adhesion_springs() {
	//get wall nodes for this cell
	vector<shared_ptr<Wall_Node>> current_cell_walls;
	this->get_Wall_Nodes_Vec(current_cell_walls);
	vector<shared_ptr<Wall_Node>> nghbr_walls_total;
	vector<shared_ptr<Wall_Node>> nghbr_walls_current;
	//int counter;
	//get all neighboring cells to this cell
	vector<shared_ptr<Cell>> neighbors;
	this->get_Neighbor_Cells(neighbors);
	for (unsigned int i = 0; i < neighbors.size(); i++) {
		neighbors.at(i)->get_Wall_Nodes_Vec(nghbr_walls_current);
		nghbr_walls_total.insert(nghbr_walls_total.end(), nghbr_walls_current.begin(), nghbr_walls_current.end());
	}	
	//	#pragma omp parallel 
	//	{
	//		#pragma omp for schedule(static,1)
	for (unsigned int i = 0; i < current_cell_walls.size(); i++) {
		//counter++;
		//cout<< counter << endl;
		//cout << "Wall node" << current_cell_walls.at(i) << endl;
		current_cell_walls.at(i)->make_connection(nghbr_walls_total);
		//cout << "connection made" << endl;
	}
	//	}
	//for all cell wall nodes
	//look at adh vector, is a nodes is in the curr cell wall
	//nodes adh vector make sure the current cell wall node
	//is in that nodes adh vector
	for (unsigned int i = 0; i < current_cell_walls.size(); i++) {
		current_cell_walls.at(i)->one_to_one_check();
	}
	return;
}
//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
//calls calc_forces
void Cell::calc_New_Forces(int Ti) {
	//cout << "cyts forces" << endl;	
	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);

#pragma omp parallel for
	for (unsigned int i = 0; i < cyts.size(); i++) {
		cyts.at(i)->calc_Forces(Ti);
	}
	//cout << "cyts done" << endl;
	//calc forces on wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);

	//cout<< "walls  forces" << endl;
#pragma omp parallel
	{
#pragma omp for schedule(static,1)
		for(unsigned int i=0; i < walls.size(); i++) {
			walls.at(i)->calc_Forces(my_tissue->get_Num_Boundary_Nodes(),Ti);
		}	
	}

	return;
}
//calls update location
//calls update cell center
//calls update wall angles
void Cell::update_Node_Locations(int Ti) {
	//update cyt nodes
	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);
#pragma omp parallel 
	{
#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < cyts.size(); i++) {
			cyts.at(i)->update_Location(Ti);
		}	
	}

	//update wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{	
#pragma omp for schedule(static,1)
		for(unsigned int i = 0; i < walls.size(); i++) {
			//cout << "update locaation" << endl;
			walls.at(i)->update_Location(Ti);
		}
	}
	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	/*if((this->life_length == 2000)) {
	 * update wall equi angles();
	}*/
	update_Wall_Equi_Angles();
	update_Wall_Angles();
	//cout << "done" << endl;
	return;
}
//=====================================================================
//==========================================
// Growth of Cell
//==========================================
//=====================================================================

void Cell::update_Cell_Progress(int& Ti) {
	//update life length of the current cell
	this->update_Life_Length();
	int cutoff_time = RECENT_DIV_NUM_FRAMES * NUM_STEPS_PER_FRAME;
	if(life_length > cutoff_time) {
		this->recent_div = false;
		this->set_MD(0);
	}

	this->Cell_Progress = static_cast<double>(this->life_length) /
		static_cast<double>(15*this->growth_rate);
		//static_cast<double>(30*this->growth_rate);

	//if (Cell_Progress > 0.6) //cout << "CELL PROGRESS INCREASING" << endl;
	bool cross_section_check = this->growing_this_cycle || !OUT_OF_PLANE_GROWTH;
	double maturity = calc_Cell_Maturity(cross_section_check);
	double max_maturity = (cross_section_check) ? 31 : 23;
	if (maturity >= num_cyt_nodes + 1 && maturity < max_maturity) {
		//cout << "cyt node added "<< endl;
		//if (cross_section_check) this->add_Cyt_Node();
		this->add_Cyt_Node();
	} 
	else { 
	////cout << "Not adding node. Mat: " << maturity << "C_P" << this->Cell_Progress << endl;
	}
	return;
}

double Cell::calc_Cell_Maturity(bool cross_section_check) { 
	double maturity;
	double finish;
	if (cross_section_check) { 
		finish = 30;
	} else { 
		finish = 22;
	}
	double exponent = (NONLINEAR_GROWTH) ? 2.0/3.0 : 1.0;
	double start = this->get_Init_Num_Nodes();
	//if (start > 25 || start < 10) { 
	////cout << "Init_Num_Nodes WRONG: " <<  start << endl;
	//}
	if (start >= finish) finish = start + 1;
	double maturity_normalized = pow(this->Cell_Progress, exponent);
	maturity = maturity_normalized * (finish-start) + start;

	return maturity;
}

void Cell::division_check() {
	vector<shared_ptr<Cell>> adh_cells;
	vector<shared_ptr<Cell>> neighbor_curr_adhesions;
	shared_ptr<Cell> this_cell = shared_from_this();
	//cout <<"Before div progress" << Cell_Progress << endl;	
	bool cross_section_check = 
		this->growing_this_cycle || !OUT_OF_PLANE_GROWTH;

	bool boundary_check = 
		(this->boundary == 0 && this->stem == 0) || BOUNDARY_DIVISION;

	bool cell_cycle_check = (this->Cell_Progress >= 1); 

	//Case where the cell divides.
	if (cross_section_check && boundary_check && cell_cycle_check) { 
		my_tissue->update_IP_Divs();
		my_tissue->update_Divs();

		//cout << "dividing cell" << this->rank <<  endl;
		//orientation of division should be 
		//fed to the division  function here
		shared_ptr<Cell> new_Cell= this->division();
		//cout << "division success" << endl;
		//cout << "Parent cell prog" << Cell_Progress << endl;
		//cout << "Sister cell prog" << new_Cell->get_Cell_Progress() << endl;
		this->my_tissue->update_Num_Cells(new_Cell);
		//setting info about new cell
		//cout << "Num cells" << this->my_tissue->get_num_cells() << endl;
		new_Cell->set_Rank(this->my_tissue->get_num_cells()-1);
		//cout << "set rank" << endl;
		//cout << "Parent rank: " << this->rank << endl;
		//cout << "Sister rank: " << new_Cell->get_Rank() << endl;
		//cout << new_Cell->get_Wall_count() << endl;
		//cout << new_Cell->get_cyt_count() << endl;
		//cout << this->get_Wall_count() << endl;
		//cout << this->get_cyt_count() << endl;
		//cout << "Parent: " << this << endl;
	////cout << "Parent progress: " << this->get_Cell_Progress() << endl;
	////cout << "New cell: " << new_Cell << endl;
	////cout << "New progress: " << new_Cell->get_Cell_Progress() << endl;
	//	//for printing mother/daughter
		this->recent_div = true;

		//layer in division function		
		//damping in division function
		//boundary needs to be figured out
		//life length set to 0 in constructor
		//life length of parent cell reset in 
		//division function
		//cyt nodes in divison function
		//wall nodes in division function
		//all cell progress set to 0 in division 
		//cell center in division function
		//cyt and wus in division function
		//growth rate in div function
		//growth direction inherited in div function
		//left corner in divison function  
		//cout << "adhesion division" << endl;
		new_Cell->update_Neighbor_Cells(this->adh_neighbors,this_cell);
		new_Cell->clear_adhesion_vectors();
		new_Cell->update_adhesion_springs();
		this->get_ADH_Neighbors_Vec(adh_cells);
		for (unsigned int i = 0; i < adh_cells.size(); i++) {
			adh_cells.at(i)->get_ADH_Neighbors_Vec(neighbor_curr_adhesions);
			adh_cells.at(i)->update_Neighbor_Cells(neighbor_curr_adhesions,new_Cell);
			adh_cells.at(i)->clear_adhesion_vectors();
			adh_cells.at(i)->update_adhesion_springs();
		}
		this->update_Neighbor_Cells(this->adh_neighbors,new_Cell);
		this->clear_adhesion_vectors();
		this->update_adhesion_springs();	
		set_Terminal(true);
	} else if (!cross_section_check && cell_cycle_check && boundary_check) { 
		my_tissue->update_Divs();
		//Case where the cell "divides out of plane"
		this->reset_Cell_Progress();
		this->reset_Life_Length();
		if (my_tissue->unifRand() < calc_My_Oop_Prob()) { 
			set_Growing_This_Cycle(false);
		} else { 
			set_Growing_This_Cycle(true);
		}

		set_Terminal(true);
	} else if (cell_cycle_check && !boundary_check) {
		this->reset_Cell_Progress();
		this->reset_Life_Length();
	}
	return;
}

double Cell::calc_Area() {
	shared_ptr<Wall_Node> curr = left_Corner;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;
	Coord a_i;
	Coord a_j;
	double area = 0;
	double curr_area = 0;
	do {
		next = curr->get_Left_Neighbor();
		a_i = curr->get_Location() - cell_center;
		a_j = next->get_Location() - cell_center;
		curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
		area += curr_area;
		curr = next;
	} while(next != orig);
	//cout << "Area: " << area << endl;
	return area;
}
double Cell::get_curr_perimeter() {
	//measure perimeter
	vector<shared_ptr<Wall_Node>> walls; 
	this->get_Wall_Nodes_Vec(walls);
	shared_ptr<Wall_Node> current;
	Coord curr_vec;
	double curr_len;
	double curr_perimeter = 0;
	shared_ptr<Wall_Node> left_neighbor;
	//#pragma omp parallel for reduction(+:curr_perimeter)
	current = walls.at(0);
	shared_ptr<Wall_Node> start = current;
	do{
		left_neighbor = current->get_Left_Neighbor();
		curr_vec = left_neighbor->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		curr_perimeter += curr_len;
		current = left_neighbor;
	}while(current != start);

	return curr_perimeter;
}
void Cell::set_perimeter(double new_perimeter){
	this->perimeter = new_perimeter;
	return;
}
void Cell::set_Growing_This_Cycle(bool gtc) {
	if(!OUT_OF_PLANE_GROWTH) {gtc = true;}
	this->growing_this_cycle = gtc;
	return;
}
double Cell::calc_My_Oop_Prob() { 
	if (UNIFORM_OOP_DIST) {
		//If we're using uniform OoP
		return OOP_PROBABILITY;
	}
	if (abs(get_Cell_Center().get_X()) < CZ_RADIUS) {
		//If we're in the center of the SAM
		if (layer == 1) return CZ_OOP_PROB_L1;
		else if (layer == 2) return CZ_OOP_PROB_L2;
		else if (layer == 3 || layer == 4) return CZ_APICAL_CORPUS_OOP_PROB;
		else return CZ_BASAL_CORPUS_OOP_PROB;
	}  
	//If we're in the periphery
	if (layer == 1) return PZ_OOP_PROB_L1;
	else if (layer == 2) return PZ_OOP_PROB_L2;
	else if (layer == 3 || layer == 4) return PZ_APICAL_CORPUS_OOP_PROB;
	else return PZ_BASAL_CORPUS_OOP_PROB;

}
void Cell::set_Init_Num_Nodes(double inn) { 
	this->init_Num_Nodes = inn;
	return;
}
void Cell::identify_Boundaries() {
	vector<shared_ptr<Wall_Node>> walls;
	double x2 = BOUNDARY_X2;
	double x1 = BOUNDARY_X1;
	double y2 = BOUNDARY_Y2;
	double y1 = BOUNDARY_Y1;
	double positive_m = (y2-y1)/(x2-x1);
	//Y intercept is the same due to evenness of initial condition.
	double y_int = positive_m*(-x1)+y1;
	double my_x, my_y;
	this->get_Wall_Nodes_Vec(walls);
	shared_ptr<Wall_Node> current;
	shared_ptr<Wall_Node> left_neighbor;
	current = walls.at(0);
	shared_ptr<Wall_Node> start = current;
	do {
		if (this->boundary == 0 || this->layer == STEM_LAYER) { 
			current->set_Is_Boundary(false);
		} else if (current->get_adh_vec().size() > 0)  {
			current->set_Is_Boundary(false);
		} else {
			//Non-adhered boundary cells beneath line
			my_x = current->get_Location().get_X();
			my_y = current->get_Location().get_Y();
			if ((my_y < (positive_m*(my_x) + y_int)) ||
					(my_y < ((-1)*positive_m*(my_x) + y_int) )) { 
				current->set_Is_Boundary(true);
			}

			
		}
		left_neighbor = current->get_Left_Neighbor();
		current = left_neighbor;
	} while(current != start);

	return;
}
void Cell::add_Wall_Node_Check(int Ti) {
	//cout << "adding a wall node" << endl;
	//#pragma omp for schedule(static,1)
	if(this->Cell_Progress < 0.05){
		//do nothing
	} else {
		shared_ptr<Wall_Node> lc;
		shared_ptr<Wall_Node> curr;
		shared_ptr<Wall_Node> temp;

		double curr_perim = this->get_curr_perimeter();
		double increase = curr_perim - this->get_perimeter();
		//cout << "curr perim " << curr_perim << endl;
		//cout << "old perim" << this->get_perimeter()<< endl;
		//cout << "increase" << increase << endl;
		this->set_perimeter(curr_perim);
		if(increase > PERIM_INCREASE){
			add_Wall_Node(Ti);
		}
	}
	return;
}
void Cell::delete_Wall_Node_Check(int Ti){
	bool repeat;
	double currAngle;
	shared_ptr<Wall_Node> currW;
	do {
		repeat = false;
		vector<pair<double,shared_ptr<Wall_Node>>> angle_wall_pairs = get_Angle_Wall_Sorted();
		for (unsigned int i = 0; i < angle_wall_pairs.size(); i++) { 
			currAngle = angle_wall_pairs.at(i).first;
			currW = angle_wall_pairs.at(i).second;
			if (currAngle < pi / 2 || currAngle > 3*pi / 2) {  
				currW->one_to_one_check();
				delete_Specific_Wall_Node(Ti, currW);
				repeat = true;
				this->get_Tissue()->inc_Num_Deleted();
			} else if ( currW->get_Updated_Tensile_Stress() < 0) { 
				//This accidentally worked nicely for 0.3.
				currW->one_to_one_check();
				delete_Specific_Wall_Node(Ti, currW);
				repeat = true;
				this->get_Tissue()->inc_Num_Deleted();
			}
			if (repeat) {

				break;
			}
		}
	} while (repeat);
	refresh_Walls();
	return;
}

void Cell::refresh_Walls() { 

	this->wall_nodes.clear();
	this->num_wall_nodes = 0;
	
	shared_ptr<Wall_Node> curr = this->left_Corner;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = curr;

	do {
		this->wall_nodes.push_back(curr);
		next = curr->get_Left_Neighbor();
		num_wall_nodes++;
		curr = next;
	} while(next != orig);
	//cout << "update equi angles" << endl;

	return;
}

void Cell::add_Wall_Node(int Ti) {

	//find node to the right of largest spring
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Wall_Node> right = NULL;
	//vector<pair<shared_ptr<Wall_Node>,double>> nodes;
	//cout  << "Find largest length" << endl;
	find_Largest_Length(right);
	//cout << "Largest found" << endl;
	shared_ptr<Wall_Node> left;
	Coord location;
	double l_thresh;
	double k_lin;
	double k_bend;
	//if(nodes.size() >0) {
	if(right != NULL){
		//find location and set neighbors for new node
		//for(int i = 0; i< nodes.size(); i++){
		//right = nodes[i].first;
		//cout << "adding node" << endl;
		left = right->get_Left_Neighbor();
		location  = (right->get_Location() + left->get_Location())*0.5;
		shared_ptr<Wall_Node> added_node = make_shared<Wall_Node>(location, this_cell, left, right);
		this->add_Wall_Node_Vec(added_node);
		double new_damping = this->get_Damping();
		right->set_Left_Neighbor(added_node);
		left->set_Right_Neighbor(added_node);
		//set the variables for the new node
		l_thresh = compute_membr_thresh(added_node);
		k_lin = compute_k_lin(added_node);
		k_bend = compute_k_bend(added_node);
		added_node->set_Damping(new_damping);
		added_node->set_K_LINEAR(k_lin);
		added_node->set_K_BEND(k_bend);
		added_node->set_membr_len(l_thresh);	
		added_node->set_added(1);
		//adhesion for the new node
		vector<shared_ptr<Cell>> neighbors;
		this->get_Neighbor_Cells(neighbors);
		//this->get_ADH_Neighbors_Vec(neighbors);
		//cout << "adh added node find closest" << endl;
		vector<shared_ptr<Wall_Node>> nghbr_walls_total;
		vector<shared_ptr<Wall_Node>> nghbr_walls_current;
		for(unsigned int i = 0; i < neighbors.size(); i++) {
			neighbors.at(i)->get_Wall_Nodes_Vec(nghbr_walls_current);
			nghbr_walls_total.insert(nghbr_walls_total.end(), nghbr_walls_current.begin(), nghbr_walls_current.end());
		}
		added_node->make_connection(nghbr_walls_total);
		//cout << "adh added node success" << endl;
		//update angles
		//should i update k_bennnnndddddd?????
		//update_Linear_Bending_Springs();
		update_Wall_Equi_Angles();
		update_Wall_Angles();
		//}
		//

		if (right->is_Boundary() || left->is_Boundary()) { 
			added_node->set_Is_Boundary(true);
			my_tissue->update_Num_Boundary_Nodes();
		} else { 
			added_node->set_Is_Boundary(false);
		}
	}
	return;
}

void Cell::refresh_Wall_Nodes_Vec() { 
	vector<shared_ptr<Wall_Node>> new_wall_nodes;
	shared_ptr<Wall_Node> curr = left_Corner;
	int new_num_wall_nodes = 0;

	do { 
		new_wall_nodes.push_back(curr);
		new_num_wall_nodes++;
		curr = curr->get_Left_Neighbor();
	} while (curr !=left_Corner);

	wall_nodes = new_wall_nodes;
	num_wall_nodes = new_num_wall_nodes;

	return;
}

void Cell::delete_Wall_Node(int Ti) {
	//shared_ptr<Wall_Node> left = NULL;
	//shared_ptr<Wall_Node> right = NULL;
	shared_ptr<Wall_Node> small = NULL;
	//vector<Cell*>neighbors;

	this->find_Smallest_Length(small);

	if (small) small->one_to_one_check();

	this->delete_Specific_Wall_Node(Ti,small);

	if (small) small->one_to_one_check();
	/*
	if(small != NULL) {
		//cout << "delete initiated" << endl;
		left = small->get_Left_Neighbor();
		right = small->get_Right_Neighbor();

		//if small is the left corner cell reassign
		if(this->left_Corner == small) {
			//cout << " set left corner" << endl;
			this->set_Left_Corner(left);
		}
		//need to make sure all nodes connected to wall
		//via adhesion are erased
		small->remove_from_adh_vecs();
		small->clear_adhesion_vec();

		//set new neighbors so nothing points at  wall
		left->set_Right_Neighbor(right);
		right->set_Left_Neighbor(left);
		//reindex
		this->wall_nodes.clear();
		this->num_wall_nodes = 0;

		shared_ptr<Wall_Node> curr = this->left_Corner;
		shared_ptr<Wall_Node> next = NULL;
		shared_ptr<Wall_Node> orig = curr;

		do {
			this->wall_nodes.push_back(curr);
			next = curr->get_Left_Neighbor();
			num_wall_nodes++;
			curr = next;
		} while (next != orig);
		//cout << "update equi angles" << endl;

		update_Wall_Equi_Angles();

		//cout << "update angles" << endl;
		update_Wall_Angles();
	}*/
	return;
}

void Cell::delete_Specific_Wall_Node(int Ti, shared_ptr<Wall_Node> wall) {
	shared_ptr<Wall_Node> left = NULL;
	shared_ptr<Wall_Node> right = NULL;
	//vector<Cell*>neighbors;
	if(wall != NULL) {
		//cout << "delete initiated" << endl;
		left = wall->get_Left_Neighbor();
		right = wall->get_Right_Neighbor();

		//if small is the left corner cell reassign
		if(this->left_Corner == wall) {
			//cout << " set left corner" << endl;
			this->set_Left_Corner(left);
		}
		//need to make sure all nodes connected to small
		//via adhesion are erased
		wall->remove_from_adh_vecs();
		wall->clear_adhesion_vec();

		//set new neighbors so nothing points at small
		left->set_Right_Neighbor(right);
		right->set_Left_Neighbor(left);
		//reindex
		this->wall_nodes.clear();
		this->num_wall_nodes = 0;

		shared_ptr<Wall_Node> curr = this->left_Corner;
		shared_ptr<Wall_Node> next = NULL;
		shared_ptr<Wall_Node> orig = curr;

		do {
			this->wall_nodes.push_back(curr);
			next = curr->get_Left_Neighbor();
			num_wall_nodes++;
			curr = next;
		} while(next != orig);
		//cout << "update equi angles" << endl;

		update_Wall_Equi_Angles();

		//cout << "update angles" << endl;
		update_Wall_Angles();
	}
	return;
}
//finds right neighbor node of smallest length on membrane
//double Cell::find_Smallest_Length(shared_ptr<Wall_Node>& right) {}
void Cell::find_Smallest_Length(shared_ptr<Wall_Node>& right) {
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 100;
	//#pragma omp parallel
	//{
	shared_ptr<Wall_Node> left_neighbor;
	double curr_len = 0;
	//#pragma omp for schedule(static,1)
	for (unsigned int i = 0; i < walls.size();i++) {
		left_neighbor = walls.at(i)->get_Left_Neighbor();
		curr_len = (walls.at(i)->get_Location()-left_neighbor->get_Location()).length();
		if(curr_len < .05){
			if(curr_len < max_len) {
				//#pragma omp critical
				max_len = curr_len;
				right = walls.at(i);
			}
		}
	}
	//}
	return;
}
//finds right neighbor node of largest length on membrane
void Cell::find_Largest_Length(shared_ptr<Wall_Node>& node) {
	vector<shared_ptr<Wall_Node>> walls; 
	this->get_Wall_Nodes_Vec(walls);
	double max_len = 0;
	//double second_max_len = 0;
	shared_ptr<Wall_Node> biggest;
	//shared_ptr<Wall_Node> second_biggest;
	//double temp_len;
	//shared_ptr<Wall_Node> temp_pointer;
	// UNUSED double theta = 0;
	// UNUSED double costheta = 0;
	double curr_len = 0;
	// UNUSED double growth_len = 0;
	Coord curr_vec;	

	//#pragma omp parallel 
	//{

	int start = my_tissue->unifRandInt(0,num_wall_nodes-1); 
	shared_ptr<Wall_Node> starter = walls.at(start);
	shared_ptr<Wall_Node> left_neighbor;
	shared_ptr<Wall_Node> current = starter;
	//double average_length;
	//#pragma omp for schedule(static,1)
	do {
		left_neighbor = current->get_Left_Neighbor();
		curr_vec = left_neighbor->get_Location() - current->get_Location();
		curr_len = curr_vec.length();	
		// UNUSED growth_len = 1;
		// UNUSED costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
		// UNUSED theta = acos( min( max(costheta,-1.0), 1.0) );
		if(curr_len > max_len){
			max_len = curr_len;
			biggest = current;
		}
		//other stuff for vector of nodes etc
		//average_length = average_length + curr_len;
		//if((theta < ADD_WALL_NODE_ANGLE_FIRST_QUAD) ||(theta > ADD_WALL_NODE_ANGLE_SECOND_QUAD)){
		//if(curr_len > MEMBR_THRESH_LENGTH) {			
		//cout << "got in here" << endl;
		//if(curr_len > second_max_len){
		//	if(curr_len > max_len){
		//temp_len = max_len;
		//temp_pointer = biggest;
		//max_len = curr_len;
		//biggest = current;
		//second_max_len = temp_len;
		//second_biggest = temp_pointer;
		//}
		//else {
		//second_max_len = curr_len;
		//second_biggest = current;
		//}
		//cout << "max_len" << max_len << endl;
		//cout << "second_max_len"<< second_max_len << endl;
		//}
		//}	

		//}
		current = left_neighbor;
	}while (left_neighbor != starter);
	//cout << "have max lens" << endl;
	//nodes.push_back(make_pair(biggest,max_len));
	//nodes.push_back(make_pair(second_biggest, second_max_len));
	//cout << "average length" << average_length/num_wall_nodes<< endl;
	node = biggest;

	return;

}
Coord Cell::compute_direction_of_highest_tensile_stress(){
	//average position of all cell wall nodes
	vector<shared_ptr<Wall_Node>> wall_nodes;
	this->get_Wall_Nodes_Vec(wall_nodes);
	Coord next_coord;
	Coord curr_coord;
	Coord direction_vec;
	shared_ptr<Wall_Node> curr = wall_nodes.at(0);
	shared_ptr<Wall_Node> orig = curr;
	shared_ptr<Wall_Node> next;
	double delta_x = 0;
	double delta_y = 0;
	double x = 0;
	double y = 0;
	double average_x;
	double average_y;
	int counter = 0;
	double curr_length;
	double strain;
	do{
		next = curr->get_Left_Neighbor();
		curr_coord = curr->get_Location();
		next_coord = next->get_Location();
		curr_length = (next->get_Location() - curr->get_Location()).length();
		delta_x = (next_coord.get_X() - curr_coord.get_X())/curr_length;
		delta_y = (next_coord.get_Y() - curr_coord.get_Y())/curr_length;
		if(delta_x < 0) {
			delta_x = delta_x*-1;
		}
		if(delta_y < 0) {
			delta_y = delta_y*-1;
		}

		strain = (curr_length - Membr_Equi_Len_Long)/Membr_Equi_Len_Long; 
		//cout << "strain" << strain << endl;
		x = x+ strain*delta_x;
		y = y+ strain*delta_y;
		counter++;
		curr = next;
	} while(next != orig);
	average_x = x/counter;
	average_y = y/counter;

	direction_vec = Coord(-average_y,average_x);
	return direction_vec;
}
/*Coord Cell::compute_point_on_line(double t){
  Coord r_0 = this->get_Cell_Center();
  Coord v = this->compute_longest_axis();
  Coord p = Coord(r_0.get_X() + v.get_X()*t, r_0.get_Y() + v.get_Y()*t);
  return p;
  }*/
void Cell::add_Cyt_Node() {
	//cout << "cyt" << endl;
	double new_damping = this->get_Damping();
	shared_ptr<Cell> this_cell = shared_from_this();
	shared_ptr<Cyt_Node> cyt = make_shared<Cyt_Node>(cell_center, this_cell);
	cyt_nodes.push_back(cyt);
	cyt->set_Damping(new_damping);

	num_cyt_nodes++;
	return;
}
//===========================================================
//==================================
// Output Functions
//==================================
//===========================================================

void Cell::print_Data_Output(ofstream& ofs) {

	ofs << "This is where data output goes" << endl;

	return;
}

int Cell::update_VTK_Indices(int& id, bool cytoplasm) {
	//cout << "ID before: " << id << endl;
	int rel_cnt = 0;

	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do { 
		curr_wall->update_VTK_Id(id);
		id++;
		for(unsigned int i = 0; i < curr_wall->get_adh_vec().size(); i++){
			rel_cnt++;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			cyt_nodes.at(i)->update_VTK_Id(id);
			id++;
		}
	}
	//cout << "ID after: " << id << endl;
	return rel_cnt;
}
void Cell::print_VTK_Adh(ofstream& ofs) {

	int my_id, nei_id;
	shared_ptr<Wall_Node> neighbor = NULL;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	vector<shared_ptr<Wall_Node>> nodes;
	do {
		for(unsigned int i = 0; i < curr_wall->get_adh_vec().size(); i++) {
			nodes = curr_wall->get_adh_vec();
			neighbor = nodes.at(i);
			if (neighbor != NULL) {
				my_id = curr_wall->get_VTK_Id();
				nei_id = neighbor->get_VTK_Id();
				ofs.flush();	
				ofs << 2 << ' ' << my_id << ' ' << nei_id << endl;
				if (abs(nei_id) > 500000 || abs(my_id) > 500000) { 
					cout << "ERR7! Cell: " << rank << ", node stats:" << endl;
					cout << "Num adhesions: " << curr_wall->get_adh_vec().size() << endl;
					cout << "Current value of i: " << i << endl;
					cout << "My address: " << curr_wall << ", VTK ID: " << my_id << endl;
					cout << "Neighbor address: " << neighbor << ", VTK ID: " << nei_id << endl;
					cout << "Neighbor belongs to: " << neighbor->get_My_Cell() << " AkA cell rank " << neighbor->get_My_Cell()->get_Rank() << endl;
				} else if (ofs.bad()) { 
					cout << "ERR8!  Bad Bit in ofs! " << endl;
					cout << "Num adhesions: " << curr_wall->get_adh_vec().size() << endl;
					cout << "Current value of i: " << i << endl;
					cout << "My address: " << curr_wall << ", VTK ID: " << my_id << endl;
					cout << "Neighbor address: " << neighbor << ", VTK ID: " << nei_id << endl;
					cout << "Neighbor belongs to: " << neighbor->get_My_Cell() << " AkA cell rank " << neighbor->get_My_Cell()->get_Rank() << endl;
				} else if (ofs.fail()) { 
					cout << "ERR9! Fail bit in ofs! " << endl;
					cout << "Num adhesions: " << curr_wall->get_adh_vec().size() << endl;
					cout << "Current value of i: " << i << endl;
					cout << "My address: " << curr_wall << ", VTK ID: " << my_id << endl;
					cout << "Neighbor address: " << neighbor << ", VTK ID: " << nei_id << endl;
					cout << "Neighbor belongs to: " << neighbor->get_My_Cell() << " AkA cell rank " << neighbor->get_My_Cell()->get_Rank() << endl;
				}
			}
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while(curr_wall != left_Corner);
	return;
}
Coord Cell::average_coordinates(){
	Coord direction = Coord(0,0);
	for(unsigned int i = 0; i <wall_nodes.size() ;i++){
		direction = direction + wall_nodes.at(i)->get_Location();
	}
	double num_nodes = (double) wall_nodes.size();
	direction = direction/num_nodes;
	return direction;
}
void Cell::print_direction_vec(ofstream& ofs){
	//average coordinates of all nodes
	Coord sum = this->average_coordinates();
	//point 1 is center + 1*direction vector
	Coord point1 = this->cell_center + sum*.05;
	Coord point2 = this->cell_center - sum*.05;
	ofs << point1.get_X() << ' ' << point1.get_Y() << ' ' << 0 << endl;
	ofs << point2.get_X() << ' ' << point2.get_Y() << ' ' << 0 << endl;

	return;
}
//LOCATIONS
void Cell::print_locations(ofstream& ofs,bool cytoplasm, int Ti) {
	///ofs << this->get_Rank() << ' ' << this->get_Layer();
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	shared_ptr<Wall_Node> orig = curr_wall;
	//int num_neighbors = num_Neighbors();
	////cout << "knows left corner" << endl;

	do {


		Coord loc = curr_wall->get_Location();
		ofs << this->get_Rank() << ' ' << this->get_Layer() << ' ' << loc.get_X() << ' ' << loc.get_Y() << ' ' << curr_wall->get_Tensile_Stress() << ' ' << Ti << endl;
		//Add 
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();

		//cout << "did it  " << count << endl;
	} while (curr_wall != orig);

	//cout << "walls worked" << endl;

	/*
	if(cytoplasm) {
		for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
			Coord loc = cyt_nodes.at(i)->get_Location();
			ofs << this->get_Rank() << ' '<< this->get_Layer() << ' ' << loc.get_X() << ' ' << loc.get_Y() << " 0 " << num_neighbors << endl;
		}
	}
	*/
	return;
}
void Cell::print_Cell_Data(ofstream& ofs, int Ti) {
	//Rank Layer Area Orientation AR(Actual,prescribed) Depth Lineage M/D N_Neigh WUS CK Ti
	
	ofs << this->get_Rank() << " " << this->get_Layer() << " ";
	ofs << calc_Area() << " ";
	vector<double> orientation_stats = calc_Orientation_Stats();
	for (unsigned int i = 0; i < 3; i++) { 
		ofs << orientation_stats.at(i) << " ";
	}
	ofs << calc_Depth() << " ";
	ofs << get_Lineage() << " ";
	ofs << recent_div_MD << " ";
	ofs << num_Neighbors() << " ";
	ofs << get_WUS_concentration() << " ";
	//get CYT concentration gets CK, 
	//this is a typo but doesn't change anything
	ofs << get_CYT_concentration() << " ";
	ofs << Ti;
	ofs << endl;

}
//Returs vector: <Longest Orientation Axis (Theta), given orientation axis (Theta), AR>
vector<double> Cell::calc_Orientation_Stats() { 
	//SPEEDUP POSSIBLE: Don't calc area separately, fold into Errera usage.
	vector<double> orientation_stats;
	vector<shared_ptr<Wall_Node>> short_axis_nodes;
	vector<shared_ptr<Wall_Node>> long_axis_nodes;
	Errera_div(short_axis_nodes);
	Coord v_1 = short_axis_nodes.at(0)->get_Location();
	Coord v_2 = short_axis_nodes.at(1)->get_Location();
	Coord short_direction = (v_2-v_1)/((v_2-v_1).length());
	//Note that this is ALWAYS taken to be with theta between 0 and 180 degrees.
	Coord long_direction = short_direction.perpVector();
	find_nodes_for_div_plane(long_direction, long_axis_nodes, 11); 
	Coord v_3 = long_axis_nodes.at(0)->get_Location();
	Coord v_4 = long_axis_nodes.at(1)->get_Location();
	double long_length = (v_4 - v_3).length();
	double short_length = (v_2 - v_1).length();
	double aspect_ratio = long_length / short_length;
	if (aspect_ratio <= 0) { 
	//cout << "NONPOSITIVE ASPECT RATIO!" << endl;
	} else if (aspect_ratio < 1) { 
	//cout << "ASPECT RATIO < 1!" << endl;
		aspect_ratio = pow(aspect_ratio,-1);
	}
	double theta_act = acos(long_direction.dot(Coord(1,0))); 
	double theta_prescribed = acos( min( max( get_growth_direction().dot(Coord(1,0)),-1.0),1.0));
	orientation_stats.push_back(theta_act);
	if (this->growth_direction == Coord(0,0)) { 
		orientation_stats.push_back(-1);
	} else { 
		orientation_stats.push_back(theta_prescribed);
	}
	orientation_stats.push_back(aspect_ratio);

	return orientation_stats;
}

double Cell::calc_Long_Length() { 
	//SPEEDUP POSSIBLE: Don't calc area separately, fold into Errera usage.
	vector<shared_ptr<Wall_Node>> short_axis_nodes;
	vector<shared_ptr<Wall_Node>> long_axis_nodes;
	Errera_div(short_axis_nodes);
	Coord v_1 = short_axis_nodes.at(0)->get_Location();
	Coord v_2 = short_axis_nodes.at(1)->get_Location();
	Coord short_direction = (v_2-v_1)/((v_2-v_1).length());
	//Note that this is ALWAYS taken to be with theta between 0 and 180 degrees.
	Coord long_direction = short_direction.perpVector();
	find_nodes_for_div_plane(long_direction, long_axis_nodes, 11); 
	Coord v_3 = long_axis_nodes.at(0)->get_Location();
	Coord v_4 = long_axis_nodes.at(1)->get_Location();
	double long_length = (v_4 - v_3).length();

	return long_length;
}



//Note that negative number returned means that you're higher up that the L1 Central cell.
//This is possible if the dome is skewed or if L2 cells creep up through L1.
double Cell::calc_Depth() { 
	double top_Cell_Y = my_tissue->get_Top_Cell_Center().get_Y();
	double my_Cell_Y = get_Cell_Center().get_Y();
	return top_Cell_Y - my_Cell_Y;
}

void Cell::set_Lineage(int parent_lineage) { 
	lineage = parent_lineage;
	return;
}

void Cell::print_VTK_Points(ofstream& ofs, int& count, bool cytoplasm) {
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	shared_ptr<Wall_Node> orig = curr_wall;
	////cout << "knows left corner" << endl;
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();
		count++;
	} while (curr_wall != orig);
	if (cytoplasm) { 
		//cout << "walls worked" << endl;
		for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
			Coord loc = cyt_nodes.at(i)->get_Location();
			ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
			count++;
		};
	}

	//cout << "points worked" << endl;
	return;
}

/*void Cell::print_VTK_Points_noCyt(ofstream& ofs, int& count) {

  shared_ptr<Wall_Node> curr_wall = left_Corner;
  shared_ptr<Wall_Node> orig = curr_wall;
////cout << "knows left corner" << endl;
do {
Coord loc = curr_wall->get_Location();
ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
//cout << "maybe cant do left neighbor" << endl;
curr_wall = curr_wall->get_Left_Neighbor();
count++;
//cout << "did it  " << count << endl;
} while (curr_wall != orig);

return;
}*/

/*void Cell::print_VTK_Scalars_Wall_Pressure(ofstream& ofs, bool cytoplasm){
  this->wall_Pressure();
  Wall_Node* curr_wall = this->left_Corner;


  double pressure;
  do{
  pressure = curr_wall->get_pressure();
  curr_wall = curr_wall->get_Left_Neighbor();
  ofs << pressure << endl;
  }while(curr_wall!= left_Corner);
  if (cytoplasm) {
  for(unsigned int i=0; i<cyt_nodes.size(); i++){
  ofs << 0 << endl;
  }
  }
  return;
  }*/
void Cell::print_VTK_Scalars_Average_Pressure(ofstream& ofs, bool cytoplasm) {
	//float pressure = this->average_Pressure();
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		//concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		//	ofs << pressure << endl;

		curr_wall = curr_wall->get_Left_Neighbor();

	} while (curr_wall != left_Corner);


	//for(unsigned int i=0;i<wall_nodes.size();i++){
	//	ofs << pressure << endl;
	//}
	if (cytoplasm) { 
		for(unsigned int i=0;i<cyt_nodes.size();i++){
			//	ofs<< pressure << endl;
		}
	}
	return;
}
/*void Cell::print_VTK_Scalars_Average_Pressure_cell(ofstream& ofs, bool cytoplasm) {
  double pressure = this->average_Pressure();
  for(unsigned int i=0;i<1;i++){
  ofs << pressure << endl;
  }
  return;
  }*/

void Cell::print_VTK_Scalars_WUS(ofstream& ofs, bool cytoplasm) {

	double concentration = 0;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();

	} while (curr_wall != left_Corner);

	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();
			ofs << concentration << endl;
		}
	}
	return;
}
void Cell::print_VTK_Scalars_CK(ofstream& ofs, bool cytoplasm) {

	double concentration = 0;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	do {
		concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();

	} while (curr_wall != left_Corner);

	if (cytoplasm) { 
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
			ofs << concentration << endl;
		}
	}
	return;
}
/*void Cell::print_VTK_Scalars_WUS_cell(ofstream& ofs, bool cytoplasm) {
  for(unsigned int i = 0; i < 1; i++) {

  double concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();

  ofs << concentration << endl;
  }
  return;
  }

  void Cell::print_VTK_Scalars_CYT(ofstream& ofs, bool cytoplasm) {

//	double concentration = 0;

shared_ptr<Wall_Node> curr_wall = left_Corner;
do {
double concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
ofs << concentration << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while (curr_wall != left_Corner);


for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
double concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
ofs << concentration << endl;
}
return;
}*/
void Cell::print_VTK_Scalars_Node(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		if(currW->get_added()==1){
			color = 30.0;
		} else {
			color = 0.0;
		}
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while (currW != left_Corner);
	if (cytoplasm) { 
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			color = 0.0;
			ofs << color << endl;
		}
	}
	return;
}

void Cell::print_VTK_Tensile_Stress(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		if(cytoplasm) currW->calc_Tensile_Stress();
		color = currW->get_Tensile_Stress();
		ofs << (isnan(color)?0:color) << endl;
		//cout << "Tensile" << color << endl;
		//cout << "Angle" << currW->get_Angle() << endl;s
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) { 
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			color = CYT_COLOR;
			ofs << color << endl;
		}
	}
	return;
}

void Cell::print_VTK_Shear_Stress(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	double color;
	do {
		color = currW->calc_Shear_Stress();
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			color = CYT_COLOR;
			ofs << color << endl;
		}
	}
	return;
}

void Cell::print_VTK_Cell_Progress(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	do {
		ofs << this->Cell_Progress << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << this->Cell_Progress << endl;
		}
	}
	return;
}

void Cell::print_VTK_Neighbors(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	update_Adh_Neighbors();	
	unsigned int color;
	int neigh = num_Neighbors();
	/*switch (neigh) {
	  case 4: 
	  color = 1;
	  break;
	  case 5: 
	  color = 2;
	  break;
	  case 6: 
	  color = 3;
	  break;
	  case 7: 
	  color = 4;
	  break;
	  default: 
	  if (neigh <= 3) { 
	  color = 0;
	  } else if (neigh >= 8) {
	  color = 5;
	  } 
	  }*/
	//Didn't do the whole lokoup table thing; so color is just the native
	//number of neighbors.
	color = neigh;

	do {
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << color << endl;
		}
	}
	return;
}

void Cell::print_VTK_Curved(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	vector<pair<double,shared_ptr<Wall_Node>>> angle_wall_pairs = get_Angle_Wall_Sorted();
	unsigned int truncate_index = static_cast<int>(angle_wall_pairs.size() * HIGH_ANGLE_DISCOUNT);
	unsigned int color;
	bool curved;
	color = 0;
	do {
		curved = false;
		for (unsigned int i = truncate_index; i < angle_wall_pairs.size(); i++) { 
			if (currW == angle_wall_pairs.at(i).second) { 
				curved = true;
			}
		}
		color = (curved ? 2 : 1);
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << 0 << endl;
		}
	}
	return;
}

/*void Cell::print_VTK_Corners(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	//cout << "Corner Print - Mark 1" << endl;
	vector<shared_ptr<Wall_Node>> corners = get_Corner_Nodes();
	//cout << "Corner Print - Mark 2" << endl;
	unsigned int color;
	bool is_corner;
	color = 0;
	do {
		is_corner = false;
		for (unsigned int i = 0; i < corners.size(); i++) { 
			if (currW == corners.at(i)) 
				is_corner = true;
		}
		color = (is_corner ? 2 : 1);
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << 0 << endl;
		}
	}
	//cout << "Corner Print - Mark 3" << endl;
	return;
}*/

void Cell::print_VTK_Boundary(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	unsigned int color;
	color = 0;
	do {
		color = (currW->is_Boundary() ? 1 : 0);
		if (color > 1) { 
			cout << "Strange Color Value: " << color << endl;
			cout << "is_Boundary() Value: " << currW->is_Boundary() << endl;
			color = 0;
		}
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();

	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << 2 << endl;
		}
	}
	return;
}

void Cell::print_VTK_MD(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	unsigned int color = this->recent_div_MD;
	do {
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << color << endl;
		}
	}
	return;
}

void Cell::print_VTK_OOP(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	bool OOP_Flag = this->growing_this_cycle;
	unsigned int color = (OOP_Flag) ? 1 : 0;
	do {
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << color << endl;
		}
	}
	return;
}

void Cell::print_VTK_Lineage(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	unsigned int color = this->lineage;
	do {
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << color << endl;
		}
	}
	return;
}

vector<pair<double,shared_ptr<Wall_Node>>> Cell::get_Angle_Wall_Sorted() { 
	vector<shared_ptr<Wall_Node>> mother_walls;
	this->get_Wall_Nodes_Vec(mother_walls);
	vector<pair<double,shared_ptr<Wall_Node>>> angle_node_pairs;
	for (unsigned int i = 0; i < mother_walls.size(); i++) { 
		angle_node_pairs.push_back(
				make_pair(mother_walls.at(i)->get_Angle() , mother_walls.at(i))
				);
	}
	//Note that sort defaults to sorting by the first element of pairs.
	sort(angle_node_pairs.begin(),angle_node_pairs.end(), greater<>());

	return angle_node_pairs; 

}

//this function loops through every node by neighbors, and determines 
//whether or not the neighboring nodes are adhered to the same cells.
//Nodes that don't match their neighbors in terms of cells adhered to
//are returned as "corner" nodes, as well as a few nodes to their left and
//right.
/*
vector<shared_ptr<Wall_Node>> Cell::get_Corner_Nodes() { 
	//cout << "Corner Calc - Mark 1" << endl;
	vector<shared_ptr<Wall_Node>> pre_corners;
	vector<shared_ptr<Wall_Node>> corners;
	shared_ptr<Wall_Node> head = get_Left_Corner();
	shared_ptr<Wall_Node> curr = head;
	shared_ptr<Wall_Node> next = NULL;
	vector<shared_ptr<Cell>> curr_adh_cells;
	vector<shared_ptr<Cell>> next_adh_cells;
	shared_ptr<Cell> curr_cell;
	bool is_new_cell;
	bool is_corner;
	//cout << "Corner Calc - Mark 2" << endl;
	//The loop must be initialized with the vector of cells that next is adhering to.
	for (unsigned int i = 0; i < curr->get_adh_vec().size(); i++) { 
		is_new_cell = true;
		curr_cell = curr->get_adh_vec().at(i)->get_My_Cell();
		for (unsigned int j = 0; j < curr_adh_cells.size(); j++) { 
			if (curr_cell == curr_adh_cells.at(j)) {
				is_new_cell = false;
				break;
			}
		}
		if (is_new_cell) {
			curr_adh_cells.push_back(curr_cell);
		}
	}

	//cout << "Corner Calc - Mark 3" << endl;
	//This loop compares the list of next's cell adhesion partners with that of current's
	//and determines which nodes are "corners".
	curr = head;
	do { 
		next_adh_cells.clear();
		next = curr->get_Left_Neighbor();
		//cout << "Corner Calc - Mark 3.1" << endl;
		for (unsigned int i = 0; i < next->get_adh_vec().size(); i++) { 
			is_new_cell = true;
			curr_cell = next->get_adh_vec().at(i)->get_My_Cell();
			for (unsigned int j = 0; j < next_adh_cells.size(); j++) { 
				if (curr_cell == next_adh_cells.at(j)) {
					is_new_cell = false;
					break;
				}
			}
			if (is_new_cell) {
				next_adh_cells.push_back(curr_cell);
			}
		}
		//Compare the cells curr and next are adhering to. If the vectors aren't the same,
		//then the node is a (pre) corner.

		//cout << "Corner Calc - Mark 3.2" << endl;
		if (curr_adh_cells.size() != next_adh_cells.size()) { 
			is_corner = true;
		} else { 
			is_corner = false;
			for (unsigned int i = 0; i < curr_adh_cells.size(); i++) { 
				if (!is_corner) { 
					for (unsigned int j = 0; j < curr_adh_cells.size(); j++) { 
						if (curr_adh_cells.at(i) == next_adh_cells.at(j)) { 
							is_corner = false;
							break;
						}
						is_corner = true;
					}
				} else { 
					break;
				}
			}
		}
		//cout << "Corner Calc - Mark 3.3" << endl;
		if (is_corner) pre_corners.push_back(curr);
		curr_adh_cells.clear();
		curr_adh_cells = next_adh_cells;
		curr = next;
	} while (next != head); 
	//cout << "Corner Calc - Mark 4" << endl;

	//The actual corners are not only the nodes where cell_neighbors change, 
	//but also CORNER_RADIUS around them.
	bool already_corner;
	for (unsigned int i = 0; i < pre_corners.size(); i++) { 
		curr = pre_corners.at(i);
		//Check nodes left of current (+1 since the immediate left
		//neighbor ought to be a pre_corner by symmetry)
		for (unsigned int j = 0; j < CORNER_RADIUS + 1; j++) {
			already_corner = false;
			for (unsigned int k = 0; k < corners.size(); k++) { 
				if (corners.at(k) == curr) { 
					already_corner = true;
				}
			}
			if (!already_corner) { 
				corners.push_back(curr);
			}
			curr = curr->get_Left_Neighbor();
		}
		//Check nodes right of current
		curr = pre_corners.at(i);
		for (unsigned int j = 0; j < CORNER_RADIUS; j++) {
			already_corner = false;
			for (unsigned int k = 0; k < corners.size(); k++) { 
				if (corners.at(k) == curr) { 
					already_corner = true;
				}
			}
			if (!already_corner) { 
				corners.push_back(curr);
			}
			curr = curr->get_Right_Neighbor();
		}
	}
	//cout << "Corner Calc - Mark 5" << endl;
	return corners;
}*/

void Cell::print_VTK_Growth_Dir(ofstream& ofs, bool cytoplasm) {
	shared_ptr<Wall_Node> currW = left_Corner;
	Coord vert(0,1);
	Coord horiz(1,0);
	Coord iso(0,0);
	int color;
	//Color-number combinations are listed in tissue.cpp
	//under print_VTK_File(...), discrete_colors lookup table
	if (this->growth_direction == vert) {
		//Vertical growth gives Red 
		color = 3;
	} else if (this->growth_direction == horiz) { 
		//Horizontal growth gives Blue
		color = 4;
	} else if (this->growth_direction == iso) { 
		//Isotropic growth gives Green
		color = 2;
	} else { 
		//Error gives white
		color = 5;
	}
	do {
		ofs << color << endl;
		currW = currW->get_Left_Neighbor();
	} while(currW != left_Corner);
	if (cytoplasm) {
		for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			ofs << color << endl;
		}
	}
	return;
}
/*		   
		   void Cell::print_VTK_Scalars_Total(ofstream& ofs, bool cytoplasm) {

//	double concentration = 0;

Wall_Node* curr_wall = left_Corner;
do {
double concentration = curr_wall->get_My_Cell()->get_total_concentration();
ofs << concentration << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while (curr_wall != left_Corner);

if (cytoplasm) { 
for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
double concentration = cyt_nodes.at(i)->get_My_Cell()->get_total_concentration();
ofs << concentration << endl;
}
}
return;
}

void Cell::print_VTK_Vectors(ofstream& ofs, bool cytoplasm) {

Wall_Node* curr_wall = left_Corner;
do {
Coord force = curr_wall->get_CytForce();
ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while(curr_wall != left_Corner);

if (cytoplasm) { 
for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
Coord force = cyt_nodes.at(i)->get_Force();
}
} 
return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

Wall_Node* curr_wall = left_Corner;
do {
Coord force = curr_wall->get_CytForce();
ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

curr_wall = curr_wall->get_Left_Neighbor();

} while(curr_wall != left_Corner);

for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
Coord force = cyt_nodes.at(i)->get_Force();
ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;
}

return;
}*/




//Debugging functions

void Cell::one_To_One_Check() { 
	shared_ptr<Wall_Node> curr_wall; 
	for (unsigned int i = 0; i < wall_nodes.size(); i++) { 
		curr_wall = wall_nodes.at(i); 
		curr_wall->one_to_one_check();
		curr_wall = curr_wall->get_Right_Neighbor();
	} 
	return;
}

void Cell::NAN_CATCH(int Ti) { 

	vector<shared_ptr<Cyt_Node>> cyts;
	this->get_Cyt_Nodes_Vec(cyts);
#pragma omp parallel 
	{
#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < cyts.size(); i++) {
			if (isnan(cyts.at(i)->get_Location().get_X())) {
				//cout << "Cyt node " << i << " in cell " << rank << " X NaN Ti=" << Ti << endl;
			}
			if (isnan(cyts.at(i)->get_Location().get_Y())) {
				//cout << "Cyt node " << i << " in cell " << rank << " Y NaN Ti=" << Ti << endl;
			}
		}	
	}

	//update wall nodes
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
#pragma omp parallel 
	{	
#pragma omp for schedule(static,1)
		for (unsigned int i = 0; i < walls.size(); i++) {
			if (isnan(walls.at(i)->get_Angle())) { 
				//cout << "Wall node " << i << " in cell " << rank << " angle NaN Ti=" << Ti << endl;
			}
			if (isnan(walls.at(i)->get_Location().get_X())) {
				//cout << "Wall node " << i << " in cell " << rank << " X NaN Ti=" << Ti << endl;
			}
			if (isnan(walls.at(i)->get_Location().get_Y())) {
				//cout << "Wall node " << i << " in cell " << rank << " Y NaN Ti=" << Ti << endl;
			}
		}	
	}
	//update cell_Center
	return;
}
//////////////////////////////////
