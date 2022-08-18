//cell_div.cpp
//========================
//Forward Declarations
//


//========================
//Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <memory>
#include <bits/stdc++.h>
#include "phys.h"
//#include "rand.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
#include "externs.h"
//==========================
//find the nodes to set up the division plane
void Cell::find_nodes_for_div_plane(Coord& orientation, vector<shared_ptr<Wall_Node>>& nodes, int search_amount) {
	//we will find the longest axis of the cell and feed in
	//the vector perpendicular to that axis to this function
	//this input  vector is the orientation vector above
	//this function finds the two cell wall nodes
	//on the current cell where the line between them is
	//closest representation to the orientation vector
	//need to know all the cell wall nodes of curr cell
	int counter = 0;
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	//this vector will hold the distance of each node
	//to the the division orientaiton vector  and will be sorted on distance
	vector<pair<double,shared_ptr<Wall_Node>>> pairs;
	//need to find two points on line that goes
	//through cell center and parallel to orientation vector
	//so we get all distances and put in the above vector
	//then we will take the smallest as first node and more
	//later about how to find the second node
	Coord center = this->get_Cell_Center();
	Coord q = Coord(center.get_X(), center.get_Y());
	Coord r = Coord(center.get_X() + orientation.get_X(), center.get_Y() + orientation.get_Y());
	Coord q_to_r = Coord(r.get_X() - q.get_X(), r.get_Y() - q.get_Y());
	Coord q_to_curr_wall;
	double q_to_r_cross_q_to_curr_wall;
	double curr_distance=0;
	shared_ptr<Wall_Node> first_node;
	shared_ptr<Wall_Node> second_node;
	for(unsigned int i = 0; i< walls.size(); i++) {

		q_to_curr_wall = Coord(walls.at(i)->get_Location().get_X() - q.get_X(), walls.at(i)->get_Location().get_Y()-q.get_Y());
		q_to_r_cross_q_to_curr_wall = q_to_r.cross(q_to_curr_wall);
		curr_distance = abs(q_to_r_cross_q_to_curr_wall)/q_to_r.length();
		pairs.push_back(make_pair(curr_distance,walls.at(i)));
	} 
	sort(pairs.begin(), pairs.end()); 
	first_node = pairs[0].second;
	shared_ptr<Wall_Node> curr = NULL;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = NULL;
	Coord a_i;
	Coord a_j;
	double area_1 = 0;
	double area_2 = 0;
	double curr_area = 0;
	//double diff = 100;
	//double curr_diff;
	counter = 1;
	bool too_small = true;
	do{
		area_1 = 0;
		area_2 = 0;
		curr = first_node;
		orig = pairs[counter].second;
		do {
			next = curr->get_Left_Neighbor();
			a_i = curr->get_Location() - cell_center;
			a_j = next->get_Location() - cell_center;
			curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
			area_1 += curr_area;
			curr = next;
		} while(next != orig);

		//cout << "area 1" << endl;
		//cout << area_1 << endl;

		curr = orig;
		orig = first_node;
		do {
			next = curr->get_Left_Neighbor();
			a_i = curr->get_Location() - cell_center;
			a_j = next->get_Location() - cell_center;
			curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
			area_2 += curr_area;
			curr = next;
		} while(next != orig);

		//cout << "area two" << endl;
		//cout << area_2 << endl;

		//curr_diff = abs(area_1 - area_2);
		//cout << "curr diff " << curr_diff << endl;
		if((area_1/area_2 > .7) && (area_1/area_2 < 1.5)){
			second_node = pairs[counter].second;
			too_small = false;
		}
		counter++;
	}while(too_small);
	//cout << "diff" << diff << endl;
	//cout << "end second" << endl;
	//cout << "first" << first << endl;
	//cout << "Second"  << second << endl;
	if((first_node == NULL)||(second_node == NULL)){

		//cout << "Did not pick up div plane nodes" << endl;

		exit(1);


	}
	nodes.clear();
	nodes.push_back(first_node);
	nodes.push_back(second_node);
	//cout << "div plane nodes size " << nodes.size() << endl;
	return;
}
void Cell::Errera_div(vector<shared_ptr<Wall_Node>>& nodes){
	vector<shared_ptr<Wall_Node>> walls;
	this->get_Wall_Nodes_Vec(walls);
	double area_1 = 0;
	double area_2 = 0;
	double curr_area = 0;
	vector<shared_ptr<Wall_Node>> first_nodes;
	vector<shared_ptr<Wall_Node>> second_nodes;
	Coord a_i;
	Coord a_j;
	shared_ptr<Wall_Node> curr;
	shared_ptr<Wall_Node> curr_partner;
	shared_ptr<Wall_Node> start_possible_partners;
	shared_ptr<Wall_Node> end_possible_partners;
	shared_ptr<Wall_Node> area_curr;
	shared_ptr<Wall_Node> area_next;
	int counter1;
	int counter2;
	for(unsigned int i = 0; i< walls.size(); i++){
		curr = walls.at(i);
		counter1 = 0;
		counter2 = 0;
		do{
			start_possible_partners = curr->get_Left_Neighbor();
			counter1++;
		} while(counter1 < 35);
		do {
			end_possible_partners = curr->get_Right_Neighbor();
			counter2++;
		}while(counter2 < 35);
		do{
			curr_partner = start_possible_partners->get_Left_Neighbor();

			area_1 = 0;
			area_2 = 0;
			//from curr to curr_partner area
			area_curr = curr;
			do {
				area_next = area_curr->get_Left_Neighbor();
				a_i = area_curr->get_Location() - cell_center;
				a_j = area_next->get_Location() - cell_center;
				curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
				area_1 += curr_area;
				area_curr = area_next;
			} while(area_next != curr_partner);

			//cout << "area 1" << endl;
			//cout << area_1 << endl;
			//from curr_partner back to curr area
			area_curr = curr_partner;
			do {
				area_next = area_curr->get_Left_Neighbor();
				a_i = area_curr->get_Location() - cell_center;
				a_j = area_next->get_Location() - cell_center;
				curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
				area_2 += curr_area;
				area_curr = area_next;
			} while(area_next != curr);

			//cout << "area two" << endl;
			//cout << area_2 << endl;

			if((area_1/area_2 > .9) && (area_1/area_2 < 1.1)){
				first_nodes.push_back(curr);
				second_nodes.push_back(curr_partner);

			}
			start_possible_partners = curr_partner;
		}while(curr_partner != end_possible_partners);
	}
	//first nodes and second nodes hold pairs close to 1/2
	//calculate distance between each and choose the smallest
	double distance;
	double max_distance = 100;
	shared_ptr<Wall_Node> first;
	shared_ptr<Wall_Node> second;
	for(unsigned int i = 0; i< first_nodes.size(); i++){
		distance = (first_nodes.at(i)->get_Location()-second_nodes.at(i)->get_Location()).length();
		if(distance < max_distance){
			max_distance = distance;
			first = first_nodes.at(i);
			second = second_nodes.at(i);
		}
	}

	nodes.clear();
	nodes.push_back(first);
	nodes.push_back(second);

	return;
}

void Cell::move_start_end_points(shared_ptr<Wall_Node> first, shared_ptr<Wall_Node> second,vector<shared_ptr<Wall_Node>>& daughter_ends){
	daughter_ends.clear();
	//to delete a node we nedd to
	//making sure nothing points to them
	//and that they don't point to anything
	//left and right neighbor will be reset in div function
	//thus when reindexed they will no longer 
	//be in the ndoes vector for this cell
	//and adhesion will be re-assigned
	//for this cell and neighbor cells
	//so adh vectors for neighboring nodes will be cleared
	//and no neighboring nodes will find them in a pair
	//thus all we need to do is clear their adhesion vec
	first->clear_adhesion_vec();
	//cout << "first clear adh" << endl;
	second->clear_adhesion_vec();
	//cout << "second clear adh" << endl;

	//move once
	shared_ptr<Wall_Node> start_daughter_one = first->get_Left_Neighbor();
	//cout << "new daughter one" << endl;
	start_daughter_one->clear_adhesion_vec();
	//cout << "daughter one clear adh" << endl;
	//move twice
	//set the new starting point of daughter cell one
	start_daughter_one = start_daughter_one->get_Left_Neighbor();
	//cout << "daught one move over" << endl;
	//start_daughter_one->clear_adhesion_vec();
	//start_daughter_one = start_daughter_one->get_Left_Neighbor();
	//move once
	shared_ptr<Wall_Node>end_daughter_one = second->get_Right_Neighbor();
	//cout << "new daughter one" << endl;
	end_daughter_one->clear_adhesion_vec();
	//cout << "clear adh" << endl;
	//move twice
	//set ending point for daughter cell one
	end_daughter_one = end_daughter_one ->get_Right_Neighbor();
	//end_daughter_one->clear_adhesion_vec();
	//end_daughter_one = end_daughter_one->get_Right_Neighbor();

	//move once
	shared_ptr<Wall_Node>start_daughter_two = second->get_Left_Neighbor();
	//cout << "daughter two start" << endl;
	start_daughter_two->clear_adhesion_vec();
	//cout << "daughter two clear" << endl;
	//move twice
	//set start point for daughter two
	start_daughter_two = start_daughter_two->get_Left_Neighbor();
	//start_daughter_two->clear_adhesion_vec();
	//start_daughter_two = start_daughter_two->get_Left_Neighbor();

	//move once
	shared_ptr<Wall_Node>end_daughter_two = first->get_Right_Neighbor();
	//cout << "daughter two end" << endl;
	end_daughter_two->clear_adhesion_vec();
	//cout << "end two clear" << endl;
	//move twice
	//set end point for daughter two
	end_daughter_two = end_daughter_two ->get_Right_Neighbor();
	//cout << "Made it through deletions" << endl;
	//end_daughter_two->clear_adhesion_vec();
	//end_daughter_two = end_daughter_two->get_Right_Neighbor();

	daughter_ends.push_back(start_daughter_one);
	daughter_ends.push_back(end_daughter_one);
	daughter_ends.push_back(start_daughter_two);
	daughter_ends.push_back(end_daughter_two);

	return;
}
Coord Cell::produce_random_vec(){
	//for random orientation
	double random = ((double) rand()) / (double) RAND_MAX;
	double diff = 1 + 1;
	double  r = random * diff;
	double x = -1 + r;
	random = ((double) rand()) / (double) RAND_MAX;
	r = random + diff;
	double y = -1 + r;
	Coord orientation = Coord(x,y);
	return orientation;
}

void Cell::find_nodes_for_div_plane_mechanical(vector<shared_ptr<Wall_Node>>& nodes){
	vector<pair<double,shared_ptr<Wall_Node>>> pairs;
	vector<shared_ptr<Wall_Node>> mother_walls;
	vector<shared_ptr<Wall_Node>> curved_walls;
	this->get_Wall_Nodes_Vec(mother_walls);
	vector<pair<double,shared_ptr<Wall_Node>>> angle_node_pairs
		= get_Angle_Wall_Sorted();

	unsigned int truncate_index = static_cast<int>(mother_walls.size() * HIGH_ANGLE_DISCOUNT);
	//for (unsigned int i = 0; i < mother_walls.size(); i++) { 
	//	angle_node_pairs.push_back(
	//	make_pair(mother_walls.at(i)->get_Angle() , mother_walls.at(i))
	//	);
	//}
	////Note that sort defaults to sorting by the first element of pairs.
	//sort(angle_node_pairs.begin(),angle_node_pairs.end(), greater<>());

	//now curved walls contains wall nodes who have the top 10% in terms of angle.

	//double theta = 0;
	//double costheta = 0;
	//double curr_len = 0;
	//double growth_len = 0;
	Coord curr_vec;	
	unsigned int counter = 0;
	double curr_stress;
	//double max_stress  =0;
	//double second_max_stress = 0;
	shared_ptr<Wall_Node> first_node;
	shared_ptr<Wall_Node> second_node;

	for(unsigned int i = 0; i < mother_walls.size();i++) {	
		mother_walls.at(i)->calc_Tensile_Stress();
		curr_stress = mother_walls.at(i)->get_Tensile_Stress();
		pairs.push_back(make_pair(curr_stress,mother_walls.at(i)));
	}


	//Avoids the unlikely situation that all of the pairs nodes could
	//be considered "too curved", and yield no div plane nodes.
	while (pairs.size() < mother_walls.size() - truncate_index) { 
		truncate_index++;
		//cout << "Unlikely situation cell_div" << endl;
	}
	sort(pairs.begin(), pairs.end(), greater<>()); 
	// Prints values to terminal to check sorting.
	//for(unsigned int i=0; i< pairs.size(); i++){
	//	//cout << pairs[i].first << endl;
	//}
	//first_node = pairs[0].second;
	first_node = NULL;
	int first_iterator = 0;
	bool curved_first_node = false;
	bool curved_second_node = false;
	shared_ptr<Wall_Node> temp;
	do { 
		curved_first_node = false;
		temp = pairs.at(first_iterator).second;
		for (unsigned int i = truncate_index; i < angle_node_pairs.size(); i++) { 
			if (temp == angle_node_pairs.at(i).second) {
				curved_first_node = true;
				break;
			}
		}
		if(!curved_first_node) { 
			first_node = temp;
			break;
		}
		first_iterator++;
	} while(!first_node);
	//cout << "Mech 9" << endl;
	if (!first_node) { 
		//cout << "WARNING: first node could not be assigned without being a corner. Defaulting";
		//cout << " to highest tensile stress." << endl;
		first_node = pairs.at(0).second;
	}
	shared_ptr<Wall_Node> curr = NULL;
	shared_ptr<Wall_Node> next = NULL;
	shared_ptr<Wall_Node> orig = NULL;
	Coord a_i;
	Coord a_j;
	double area_1 = 0;
	double area_2 = 0;
	double curr_area = 0;
	//double diff = 100;
	//double curr_diff;
	counter = 1;
	bool too_small = true;
	//IN PROGRESS - CHECK
	do {
		curved_second_node = false;
		area_1 = 0;
		area_2 = 0;
		curr = first_node;
		//Orig is the candidate for second_node.
		orig = pairs[counter].second;
		for (unsigned int i = truncate_index; i < angle_node_pairs.size(); i++) { 
			if(orig == angle_node_pairs.at(i).second) {
				curved_second_node = true;
			}
		}
		do {
			next = curr->get_Left_Neighbor();
			a_i = curr->get_Location() - cell_center;
			a_j = next->get_Location() - cell_center;
			curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
			area_1 += curr_area;
			curr = next;
		} while(next != orig);

		//cout << "\narea 1 " << area_1 << "; ";

		curr = orig;
		//Orig is now first_node
		orig = first_node;
		do {
			next = curr->get_Left_Neighbor();
			a_i = curr->get_Location() - cell_center;
			a_j = next->get_Location() - cell_center;
			curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
			area_2 += curr_area;
			curr = next;
		} while(next != orig);

		//cout << "; area 2 " << area_2 << endl;

		if((area_1/area_2 > 0.9) && (area_1/area_2 < 1.1111)){
			second_node = pairs[counter].second;
			too_small = false;
		}
		if(second_node == first_node) {
			//cout << "Nodes must be distinct, second = first in mech_div" << endl;
			second_node = NULL;
			too_small = true;
		}
		counter++;
		if (counter >= pairs.size()) { 
			//cout << "Skipped all nodes.  Returning to origin and relaxing angle tolerance." << endl;
			truncate_index++;
			counter = 0;
		}
	} while(too_small || curved_second_node);
	if((first_node == NULL)||(second_node == NULL)){

		//cout << "Did not select second div plane node." << endl;

		exit(1);
	}

	nodes.clear();
	nodes.push_back(first_node);
	nodes.push_back(second_node);

	return;
}

Coord Cell::find_orientation_circ() {
	//Find neighbors in this layer
	Coord orientation;
	vector<shared_ptr<Cell>> my_adh_neighbors;
	vector<shared_ptr<Cell>> adh_neighbors_this_layer;
	this->get_ADH_Neighbors_Vec(adh_neighbors);

	for (unsigned int i = 0; i < my_adh_neighbors.size(); i++) { 
		if (my_adh_neighbors.at(i)->get_Layer() == this->get_Layer()) {
			adh_neighbors_this_layer.push_back(my_adh_neighbors.at(i));
		}
	}
	//Check to see if there are two of them.
	
	if (adh_neighbors_this_layer.size() == 2) { 
		//Circle fit division plane orientation

		Coord here = this->get_Cell_Center();
		Coord right_loc = adh_neighbors_this_layer.at(0)->get_Cell_Center();
		Coord left_loc = adh_neighbors_this_layer.at(1)->get_Cell_Center();
		double x1, y1, x2, y2, x3, y3;
		x1 = here.get_X();
		y1 = here.get_Y();
		x2 = right_loc.get_X();
		y2 = right_loc.get_Y();
		x3 = left_loc.get_X();
		y3 = left_loc.get_Y();

		double x12 = x1 - x2; 
		double x13 = x1 - x3; 

		double y12 = y1 - y2; 
		double y13 = y1 - y3; 

		double y31 = y3 - y1; 
		double y21 = y2 - y1; 

		double x31 = x3 - x1; 
		double x21 = x2 - x1; 
		// x1^2 - x3^2 
		double sx13 = pow(x1, 2) - pow(x3, 2); 

		// y1^2 - y3^2 
		double sy13 = pow(y1, 2) - pow(y3, 2); 

		double sx21 = pow(x2, 2) - pow(x1, 2); 
		double sy21 = pow(y2, 2) - pow(y1, 2); 

		double f = ((sx13) * (x12) 
				+ (sy13) * (x12) 
				+ (sx21) * (x13) 
				+ (sy21) * (x13)) 
			/ (2 * ((y31) * (x12) - (y21) * (x13))); 
		double g = ((sx13) * (y12) 
				+ (sy13) * (y12) 
				+ (sx21) * (y13) 
				+ (sy21) * (y13)) 
			/ (2 * ((x31) * (y12) - (x21) * (y13))); 

		//double c = -pow(x1, 2) - pow(y1, 2) - 2 * g * x1 - 2 * f * y1; 
		//In case I need this, the radius can be calculate with c.

		// eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0 
		// where centre is (h = -g, k = -f) and radius r 
		// as r^2 = h^2 + k^2 - c 
		double h = -g; 
		double k = -f; 

		orientation = this->get_Cell_Center() - Coord(h,k);
	} else {
		//Tissue base approximation via semicircular tissue
		vector<shared_ptr<Cell>> cells;
		this->my_tissue->get_Cells(cells);
		Coord tissue_base_location = cells.at(TISSUE_BASE)->get_Cell_Center();
		Coord pre_orientation = this->get_Cell_Center() - tissue_base_location;
		orientation = pre_orientation / pre_orientation.length();
	}


	return orientation;
}

shared_ptr<Cell> Cell::division() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as parent cell
	//	this function will create a sister cell
	//	and return it to the tissue
	//cout << "New Cell" << endl;
	shared_ptr<Wall_Node> curr_wall = left_Corner;
	/*cout << "Mother angles before: " << endl;
	  do { 
	  cout << curr_wall->get_Angle() << ", ";
	  curr_wall = curr_wall->get_Left_Neighbor();
	  } while (curr_wall != left_Corner);
	  cout << endl;
	  */
	shared_ptr<Cell> sister = make_shared<Cell> (this->my_tissue);
	sister->set_Layer(this->layer);
	sister->set_Lineage(this->lineage);
	sister->set_MD(2);
	this->set_MD(1);

	this->life_length = 0;	
	this->Cell_Progress = 0;		
	sister->reset_Cell_Progress();
	sister->reset_Life_Length();
	//this is the vector for the desired division orientation
	Coord orientation;
	//for layer and chemical division decision a vector called orientation
	//is fed into the find nodes for div plane function
	//find nodes for div plan takes as input orientation vector, all wall nodes of mother cell
	//a hard coded search amount (could this be better?)
	//return a vector with two pointers to two wall nodes nodes indicating where the cell wall 
	//should be built
	//orientation by layer
	//if((this->layer == 1) || (this->layer == 2)){
	//	orientation = Coord(0,1);
	//}
	//else{
	//	orientation = Coord(1,0);
	//}

	//orientation by chemical concentration
	vector<shared_ptr<Wall_Node>> nodes;
	//cout << "Nodes before " << nodes.size() << endl;

	bool forced_anticlinal = ((this->layer == 1)||(this->layer==2))
		&& L1_L2_FORCED_ANTICLINAL_DIV;

	vector<shared_ptr<Cell>> cells;
	my_tissue->get_Cells(cells);
	//Coord tissue_base_location = cells.at(TISSUE_BASE)->get_Cell_Center();

	if (forced_anticlinal) { 

		//orientation = Coord(0,1);
		//Coord pre_orientation = this->get_Cell_Center() - tissue_base_location;
		//orientation = pre_orientation / pre_orientation.length();
		orientation = find_orientation_circ();
		find_nodes_for_div_plane(orientation,nodes,11);


	} else { 
		switch (DIV_MECHANISM) { 
			case 1: 
				//Errera's rule
				Errera_div(nodes);
				break;
			case 2:
				//Chemically Driven Division
				if(HILL_PROB) { 
					if(my_tissue->unifRand() < hill_Prob()) {
						//hill_Prob is the probability of periclinal.
						orientation = Coord(1,0); //periclinal division
					} else {
						//Else anticlinal.
						orientation = Coord(0,1);
					}
				} else { 
					if((this->get_CYT_concentration() > this->get_WUS_concentration())){
						orientation = Coord(1,0);
					} else { 
						orientation = Coord(0,1);
					}
				}
				find_nodes_for_div_plane(orientation,nodes,11);
				break;
			case 3:
				//Mechanical Division
				find_nodes_for_div_plane_mechanical(nodes); 
				break;
			case 4:
				//Merged Mechanism
				if(my_tissue->unifRand() < hill_Prob()) {
					//hill_Prob is the probability of periclinal.
					orientation = Coord(1,0); //periclinal division
					find_nodes_for_div_plane(orientation,nodes,11);
				} else {
					//Else mechanical.
					find_nodes_for_div_plane_mechanical(nodes); 
				}
				break;
			case 5:
				//Layer-Specific Mechanism:
				//L3-4 experience CAE-M (Case 3) and L5-6 experience CED (Case 2)
				//Note that if forced L1 anticlinal is off, this defaults
				//to CAE-M for layers 1 and 2.
				if (this->layer <= 4) { 
					find_nodes_for_div_plane_mechanical(nodes);
				} else { 
					if (HILL_PROB) { 
						if (my_tissue->unifRand() < hill_Prob()) { 
							//hill_Prob is the probability of periclinal.
							orientation = Coord(1,0); //periclinal division
						} else { 
							//Else anticlinal
							orientation = Coord(0,1);
						}
					} else { 
						if ((this->get_CYT_concentration() > this->get_WUS_concentration())) { 
							orientation = Coord(1,0); 
						} else { 
							orientation = Coord(0,1);
						}
					}
					find_nodes_for_div_plane(orientation,nodes,11);
				}
				break;
			default:
				//cout << "DIV_MECH NOT ENTERED. Exiting..." << endl;
				exit(1);
		}
	}

	my_tissue->update_Divplane_Vector(nodes.at(0)->get_Location() - nodes.at(1)->get_Location(), 
			this->layer);
	//cout << "Nodes after " << nodes.size() << endl;
	//finds node on one side of cell
	//vector<shared_ptr<Wall_Node>> nodes;
	shared_ptr<Wall_Node> first = nodes.at(0);
	//cout <<"check one" << endl;
	shared_ptr<Wall_Node> second = nodes.at(1);
	//shared_ptr<Wall_Node> first = first_node;
	//shared_ptr<Wall_Node> second = second_node;
	vector<shared_ptr<Wall_Node>> daughter_ends;
	move_start_end_points(first, second, daughter_ends);
	//move start end points returns daughter_ends as
	//daughter_ends = {start_one, end_one, start_two, end_two}



	shared_ptr<Wall_Node> start_daughter_one = daughter_ends.at(0);
	shared_ptr<Wall_Node> start_daughter_two = daughter_ends.at(2);
	shared_ptr<Wall_Node> end_daughter_one = daughter_ends.at(1);
	shared_ptr<Wall_Node> end_daughter_two = daughter_ends.at(3);
	double daughter_one_length = (start_daughter_one->get_Location() - end_daughter_one->get_Location()).length();
	double daughter_two_length = (start_daughter_two->get_Location() - end_daughter_two->get_Location()).length();
	//cout << "Daughter_Lengths: One = " << daughter_one_length << " Two = " << daughter_two_length << endl;
	//cout << "d one length" << daughter_one_length << endl;
	//cout << "d two lenght" << daughter_two_length << endl;

	//this is how many new nodes will be made for each cell
	int total_num_one = static_cast<int>(daughter_one_length/(Membr_Equi_Len_Short*3));
	int total_num_two = static_cast<int>(daughter_two_length/(Membr_Equi_Len_Short*3));
	//cout << "Total num one" << total_num_one << endl;
	//cout << "Total num two" << total_num_two << endl;
	//this is how much change there is in the x direction of divison plane for each cell
	double x_one = end_daughter_one->get_Location().get_X() - start_daughter_one->get_Location().get_X();
	double x_two = end_daughter_two->get_Location().get_X() - start_daughter_two->get_Location().get_X();

	//this is how much change there is in the y direction of division plane for each cell
	double y_one = end_daughter_one->get_Location().get_Y() - start_daughter_one->get_Location().get_Y();
	double y_two = end_daughter_two->get_Location().get_Y() - start_daughter_two->get_Location().get_Y();

	//this is how much x should increment by when
	//making new cell wall nodes for 
	//each cell
	double delta_x_one = (double)x_one/(double)total_num_one;
	double delta_x_two = (double)x_two/(double)total_num_two;


	//this is how much y should increment by when 
	//making new cell wall nodes for 
	//each cell
	double delta_y_one = (double)y_one/(double)total_num_one;
	double delta_y_two = (double)y_two/(double)total_num_two;

	//some variables to be used for making new nodes
	shared_ptr<Wall_Node> currW = start_daughter_one;
	Coord location = start_daughter_one->get_Location();
	shared_ptr<Wall_Node> prevW = start_daughter_one;
	shared_ptr<Cell> this_cell = shared_from_this();
	double curr_X = start_daughter_one->get_Location().get_X()+2*delta_x_one;
	double curr_Y = start_daughter_one->get_Location().get_Y()+2*delta_y_one;
	//this counter was for debugging

	int counter = 0;

	//cout << "Make daughter one new cell wall nodes" << endl;
	double total_num_one_minus_four = total_num_one -4;
	for(unsigned int i = 0; i < total_num_one_minus_four; i++){
		curr_X = 
			curr_X + delta_x_one;
		curr_Y = curr_Y + delta_y_one;
		location = Coord(curr_X,curr_Y);
		//node made 
		//cout << "making new wall node here" << endl;
		currW = make_shared<Wall_Node>(location,this_cell);
		//set neighbors
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;
		//debugging
		counter++;
		//cout << counter << endl;	
	}
	currW->set_Right_Neighbor(end_daughter_one);
	end_daughter_one->set_Left_Neighbor(currW);
	this->left_Corner = start_daughter_one;
	//cout << "one new" << counter << endl;
	//cout << "Make daughter two new cell wall nodes" << endl;
	counter = 0;
	currW = start_daughter_two;
	location = start_daughter_two->get_Location();
	prevW = start_daughter_two;
	curr_X = location.get_X()+2*delta_x_two;
	curr_Y = location.get_Y()+2*delta_y_two;
	double total_num_two_minus_four = total_num_two - 4;
	for(unsigned int i = 0; i< total_num_two_minus_four; i++){

		curr_X = curr_X + delta_x_two;
		curr_Y =
			curr_Y + delta_y_two;
		location = Coord(curr_X,curr_Y);
		//node made
		//cout << "make new node here" << endl;
		currW = make_shared<Wall_Node>(location,sister);
		//set neighbors		
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		prevW = currW;	
		//debugging
		counter++;
		//cout<< counter << endl;

	}
	currW->set_Right_Neighbor(end_daughter_two);
	end_daughter_two->set_Left_Neighbor(currW);
	sister->set_Left_Corner(start_daughter_two);
	//cout << "two new" << counter << endl;
	//cout << "Made wall nodes" << endl;
	//here i am refilling
	//the wall node vectors for each cell 
	//and setting private member variables for 
	//each new wall node
	double new_damping = this->get_Damping();	
	sister->set_Damping(new_damping);
	shared_ptr<Wall_Node> curr = sister->get_Left_Corner();
	shared_ptr<Wall_Node> orig = curr;
	double l_thresh = 0;
	double k_lin = 0;
	double k_bend = 0;
	int number_nodes_B = 0;
	//cout << "Cell two counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		curr->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(curr);
		k_lin = compute_k_lin(curr);
		k_bend = compute_k_bend(curr);
		curr->set_membr_len(l_thresh);
		curr->set_K_LINEAR(k_lin);
		curr->set_K_BEND(k_bend);

		sister->add_Wall_Node_Vec(curr);
		curr = curr->get_Left_Neighbor();
		//cout << number_nodes_B << endl;
		if(number_nodes_B > 250) {
			//cout << "blowing up!" << endl;
			exit(1);
		}
		//cout << curr->get_Location() << endl;
		if(curr == NULL) {
			//cout << "notlinked" << endl;
			exit(1);
		}
	} while (curr != orig);
	//cout << "Done counting cell two" << endl;
	sister->set_Wall_Count(number_nodes_B);
	//cout << "Sister: " << number_nodes_B << endl;

	this->wall_nodes.clear();

	curr = this->left_Corner;
	orig = this->left_Corner;
	int number_nodes_A = 0;
	//cout << "Cell one counting" << endl;
	do {

		number_nodes_A++;
		curr->update_Cell(this_cell);
		curr->set_Damping(new_damping);
		l_thresh = compute_membr_thresh(curr);
		k_lin = compute_k_lin(curr);
		k_bend = compute_k_bend(curr);
		curr->set_membr_len(l_thresh);
		curr->set_K_LINEAR(k_lin);
		curr->set_K_BEND(k_bend);


		this->add_Wall_Node_Vec(curr);
		curr = curr->get_Left_Neighbor();
		//cout << number_nodes_A << endl;
		if(number_nodes_A > 300){
			//cout << "blowing up" << endl;
			exit(1);
		}
		//cout << curr->get_Location() << endl;
		if(curr == NULL) {

			//cout << "notlinked" << endl;
			exit(1);
		}
	} while(curr != orig);
	//cout << "Done counting cell one" << endl;	
	this->set_Wall_Count(number_nodes_A);
	//cout << "Parent nodes: " << number_nodes_A << endl;	

	//update cell level variables
	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();;
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	this->update_Cell_Center();
	sister->update_Cell_Center();
	Coord L1_AVG = this->get_Tissue()->Compute_L1_AVG();
	double acd = this->get_Tissue()->get_Avg_Cell_Diameter();
	this->calc_WUS(L1_AVG, 2 * acd);
	sister->calc_WUS(L1_AVG, 2 * acd);
	this->calc_CK(L1_AVG, 3 * acd);
	sister->calc_CK(L1_AVG, 3 * acd);

	//Set sister growth rate at random
	/*
	if((this->layer == 1)||(this->layer == 2)) {
		sister->set_growth_direction(Coord(1,0));
	} else if(rand()% 100 < 30) {
		sister->set_growth_direction(Coord(0,0));
	} else {
		if(rand()% 100 < 36){
			sister->set_growth_direction(Coord(1,0));
		} else {
			sister->set_growth_direction(Coord(0,1));
		}
	}
	*/
	//Sister updates growth rate based on chemical distribution.
	this->update_growth_direction();
	sister->update_growth_direction();

	
	if(my_tissue->unifRand() < calc_My_Oop_Prob()) { 
		sister->set_Growing_This_Cycle(false);
	} else { 
		sister->set_Growing_This_Cycle(true);
	}

	this->set_growth_rate(true);
	sister->set_growth_rate(true);


	//cout << "Reassign cyt nodes" << endl;
	vector<shared_ptr<Cyt_Node>> temp_cyts;
	this->get_Cyt_Nodes_Vec(temp_cyts);
	this->cyt_nodes.clear();
	this->num_cyt_nodes = 0;

	//this is the center of the new cell wall for each cell

	Coord vector_from_center;
	Coord new_location;
	//for debugging
	counter = 0;
	Coord center_1 = (end_daughter_one->get_Location() - start_daughter_one->get_Location())*.5;
	Coord center_2 = (end_daughter_two ->get_Location()-start_daughter_two->get_Location())*.5;
	double length_1;
	double length_2;
	for(unsigned int i = 0; i < temp_cyts.size(); i++) {
		length_1 = (temp_cyts.at(i)->get_Location() - this->get_Cell_Center()).length();
		length_2 = (temp_cyts.at(i)->get_Location() - sister->get_Cell_Center()).length();
		if(length_1 < length_2) {
			temp_cyts.at(i)->update_Cell(this_cell);
			this->update_cyt_node_vec(temp_cyts.at(i));
			//cout << temp_cyts.at(i)->get_Location() << endl;
			//this->Cell_Progress++;
			//cout << "this" << endl;
			//cout << "my cell: " << temp_cyts.at(i)->get_My_Cell() << endl;
			//if(length_1 < 4){
			//	temp_cyts.at(i)->new_location(this->get_Cell_Center() + (cyt_nodes.at(i)->get_Location() - this->get_Cell_Center())*.7);
			//}

		} else {
			temp_cyts.at(i)->update_Cell(sister);
			sister->update_cyt_node_vec(temp_cyts.at(i));
			//cout << temp_cyts.at(i)->get_Location() << endl;
			//sister->update_Cell_Progress();
			//cout << "sister" << endl;
			//cout << "My cell" << temp_cyts.at(i)->get_My_Cell() << endl;
			//if(length_2 < 4){
			//	temp_cyts.at(i)->new_location(sister->get_Cell_Center() + (cyt_nodes.at(i)->get_Location() - sister->get_Cell_Center())*.7);
			//}

		}

		//counter++;
		//cout << "counter" <<  counter << endl;
	}
	this->set_Init_Num_Nodes(get_cyt_count());
	sister->set_Init_Num_Nodes(sister->get_cyt_count());
	//move cyt nodes that are too close 
	//to the new cell wall
	//cout << "sister" << sister->get_cyt_count()<< endl;
	//cout << "this" << this->get_cyt_count()<< endl;
	//cout << "Move cyts here" << endl;
	//Coord center_1 = Coord(0,0);
	//Coord center_2 = Coord(0,0);
	this->move_cyt_nodes(center_1);
	sister->move_cyt_nodes(center_2);

	//cout << "Finished deleting old cyt nodes" << endl;

	curr_wall = left_Corner;
	//cout << "Mother angles after: " << endl;
	do { 
		//cout << curr_wall->get_Angle() << ", ";
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);
	//cout << endl;

	shared_ptr<Wall_Node> sis_left_Corner = sister->get_Left_Corner();  
	curr_wall = sis_left_Corner;
	//cout << "Sister angles after: " << endl;
	do { 
		//cout << curr_wall->get_Angle() << ", ";
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != sis_left_Corner);
	//cout << endl;

	if (my_tissue->unifRand() < OOP_PROBABILITY) {  
		set_Growing_This_Cycle(false);
	} else { 
		set_Growing_This_Cycle(true);
	}

	return sister;
}
void Cell::move_cyt_nodes(Coord center_point){
	Coord vector_from_center;
	//double length_from_center_pt;
	Coord location;
	for(unsigned int i=0; i<cyt_nodes.size(); i++){
		//length_from_center_pt = (cyt_nodes.at(i)->get_Location()-center_point).length();
		vector_from_center = cyt_nodes.at(i)->get_Location() - this->cell_center;

		//if(length_from_center_pt< 4) {

		location = cell_center + vector_from_center*.5;
		cyt_nodes.at(i)->new_location(location);
		//cout << cyt_nodes.at(i)->get_Location() << endl;
		//}
	}
	return;
}
//////=======================================================================================================
//end of cell_div.cpp










//DEPRECATED CODE

//Mechanical division V1
/*void Cell::find_nodes_for_div_plane_mechanical(vector<shared_ptr<Wall_Node>>& nodes){
  vector<pair<double,shared_ptr<Wall_Node>>> pairs;
  vector<shared_ptr<Wall_Node>> mother_walls;
  this->get_Wall_Nodes_Vec(mother_walls);

  double theta = 0;
  double costheta = 0;
  double curr_len = 0;
  double growth_len = 0;
  Coord curr_vec;	
  int counter = 0;
  double curr_stress;
//double max_stress  =0;
//double second_max_stress = 0;
shared_ptr<Wall_Node> first_node;
shared_ptr<Wall_Node> second_node;
if(this->growth_direction != Coord(0,0)){
for(unsigned int i = 0; i < mother_walls.size();i++) {	
curr_vec = mother_walls.at(i)->get_Left_Neighbor()->get_Location() - mother_walls.at(i)->get_Location();
curr_len = curr_vec.length();	
growth_len = this->growth_direction.length();
costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
theta = acos( min( max(costheta,-1.0), 1.0) );
if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
curr_stress = mother_walls.at(i)-> calc_Tensile_Stress();
pairs.push_back(make_pair(curr_stress,mother_walls.at(i)));
} 
}
}
else{
for(unsigned int i = 0; i < mother_walls.size();i++) {	
curr_stress = mother_walls.at(i)-> calc_Tensile_Stress();
pairs.push_back(make_pair(curr_stress,mother_walls.at(i)));
}
}

sort(pairs.begin(), pairs.end(), greater<>()); 
for(unsigned int i=0; i< pairs.size(); i++){
cout << pairs[i].first << endl;
}
first_node = pairs[0].second;
shared_ptr<Wall_Node> curr = NULL;
shared_ptr<Wall_Node> next = NULL;
shared_ptr<Wall_Node> orig = NULL;
Coord a_i;
Coord a_j;
double area_1 = 0;
double area_2 = 0;
double curr_area = 0;
//double diff = 100;
//double curr_diff;
counter = 1;
bool too_small = true;
do{
area_1 = 0;
area_2 = 0;
curr = first_node;
orig = pairs[counter].second;
do {
next = curr->get_Left_Neighbor();
a_i = curr->get_Location() - cell_center;
a_j = next->get_Location() - cell_center;
curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
area_1 += curr_area;
curr = next;
} while(next != orig);

cout << "area 1" << endl;
cout << area_1 << endl;

curr = orig;
orig = first_node;
do {
	next = curr->get_Left_Neighbor();
	a_i = curr->get_Location() - cell_center;
	a_j = next->get_Location() - cell_center;
	curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
	area_2 += curr_area;
	curr = next;
} while(next != orig);

cout << "area two" << endl;
cout << area_2 << endl;

//curr_diff = abs(area_1 - area_2);
//cout << "curr diff " << curr_diff << endl;
if((area_1/area_2 > .7) && (area_1/area_2 < 1.5)){
	second_node = pairs[counter].second;
	too_small = false;
}
counter++;
}while(too_small);
//cout << "diff" << diff << endl;
//cout << "end second" << endl;
//cout << "first" << first << endl;
//cout << "Second"  << second << endl;
if((first_node == NULL)||(second_node == NULL)){

	//cout << "Did not pick up div plane nodes" << endl;

	exit(1);


}
nodes.clear();
nodes.push_back(first_node);
nodes.push_back(second_node);
//cout << "div plane nodes size " << nodes.size() << endl;

return;
}*/

//Mechanical division V2
/*void Cell::find_nodes_for_div_plane_mechanical(vector<shared_ptr<Wall_Node>>& nodes){
  vector<pair<double,shared_ptr<Wall_Node>>> pairs;
  vector<shared_ptr<Wall_Node>> mother_walls;
  this->get_Wall_Nodes_Vec(mother_walls);


//cout << "Mech 1" << endl;
//contains a subset of mother_walls, hopefully without corners.
vector<shared_ptr<Wall_Node>> curved_walls;
vector<double> angles;
double truncate_angle = 0;
unsigned int truncate_index = static_cast<int>(mother_walls.size() * HIGH_ANGLE_DISCOUNT);
for (unsigned int i = 0; i < mother_walls.size(); i++) { 
//mother_walls.at(i)->update_Angle();
angles.push_back(mother_walls.at(i)->get_Angle());
cout << angles.at(i); 
}
//cout << endl;
//cout << "Mech 2" << endl;
sort(angles.begin(),angles.end());
//cout << "Mech 3" << endl;
truncate_angle = angles.at(truncate_index);
//cout << "Mech 5" << endl;
for (unsigned int i = 0; i < mother_walls.size(); i++) { 
if (mother_walls.at(i)->get_Angle() >= truncate_angle) { 
curved_walls.push_back(mother_walls.at(i));
}
if (curved_walls.size() > static_cast<float>(mother_walls.size())*(1.0-HIGH_ANGLE_DISCOUNT) + 1) { 
cout << "Too many walls with near maximal angle. Breaking search.\n";
cout << "This is expected for a circular cell.\n";
break;
}
}

//cout << "Mech 6"<< endl;
//now curved walls contains wall nodes who have the top 10% in terms of angle.

double theta = 0;
double costheta = 0;
double curr_len = 0;
double growth_len = 0;
Coord curr_vec;	
int counter = 0;
double curr_stress;
//double max_stress  =0;
//double second_max_stress = 0;
shared_ptr<Wall_Node> first_node;
shared_ptr<Wall_Node> second_node;
if(this->growth_direction != Coord(0,0)) {
for(unsigned int i = 0; i < mother_walls.size();i++) {	
curr_vec = mother_walls.at(i)->get_Left_Neighbor()->get_Location() - mother_walls.at(i)->get_Location();
curr_len = curr_vec.length();	
growth_len = this->growth_direction.length();
costheta = growth_direction.dot(curr_vec)/(curr_len*growth_len);
theta = acos( min( max(costheta,-1.0), 1.0) );
if((theta < ANGLE_FIRST_QUAD) || (theta > ANGLE_SECOND_QUAD)){
curr_stress = mother_walls.at(i)-> calc_Tensile_Stress();
pairs.push_back(make_pair(curr_stress,mother_walls.at(i)));
} 
}
} else {
for(unsigned int i = 0; i < mother_walls.size();i++) {	
curr_stress = mother_walls.at(i)-> calc_Tensile_Stress();
pairs.push_back(make_pair(curr_stress,mother_walls.at(i)));
}
}
//cout << "Mech 7" << endl;

sort(pairs.begin(), pairs.end(), greater<>()); 
// Prints values to terminal to check sorting.
//for(unsigned int i=0; i< pairs.size(); i++){
//	//cout << pairs[i].first << endl;
//}
//cout << "Mech 8" << endl;
//first_node = pairs[0].second;
first_node = NULL;
int first_iterator = 0;
bool curved_first_node = false;
shared_ptr<Wall_Node> temp;
do { 
	curved_first_node = false;
	temp = pairs.at(first_iterator).second;
	for (unsigned int i = 0; i < curved_walls.size(); i++) { 
		if (temp == curved_walls.at(i)) {
			curved_first_node = true;
			break;
		}
	}
	if(!curved_first_node) { 
		first_node = temp;
		break;
	}
	first_iterator++;
} while(!first_node);
//cout << "Mech 9" << endl;
if (!first_node) { 
	//cout << "WARNING: first node could not be assigned without being a corner. Defaulting";
	//cout << " to highest tensile stress." << endl;
	first_node = pairs.at(0).second;
}
shared_ptr<Wall_Node> curr = NULL;
shared_ptr<Wall_Node> next = NULL;
shared_ptr<Wall_Node> orig = NULL;
Coord a_i;
Coord a_j;
double area_1 = 0;
double area_2 = 0;
double curr_area = 0;
//double diff = 100;
//double curr_diff;
counter = 1;
bool too_small = true;
do{
	area_1 = 0;
	area_2 = 0;
	curr = first_node;
	orig = pairs[counter].second;
	do {
		next = curr->get_Left_Neighbor();
		a_i = curr->get_Location() - cell_center;
		a_j = next->get_Location() - cell_center;
		curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
		area_1 += curr_area;
		curr = next;
	} while(next != orig);

	//cout << "\narea 1 " << area_1 << "; ";

	curr = orig;
	orig = first_node;
	do {
		next = curr->get_Left_Neighbor();
		a_i = curr->get_Location() - cell_center;
		a_j = next->get_Location() - cell_center;
		curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
		area_2 += curr_area;
		curr = next;
	} while(next != orig);

	//cout << "; area 2 " << area_2 << endl;

	//curr_diff = abs(area_1 - area_2);
	//cout << "curr diff " << curr_diff << endl;
	if((area_1/area_2 > .7) && (area_1/area_2 < 1.5)){
		second_node = pairs[counter].second;
		too_small = false;
	}
	counter++;
}while(too_small);
//cout << "Mech 10" << endl;
//cout << "diff" << diff << endl;
//cout << "end second" << endl;
//cout << "first" << first << endl;
//cout << "Second"  << second << endl;
if((first_node == NULL)||(second_node == NULL)){

	//cout << "Did not pick up div plane nodes" << endl;

	exit(1);


}
//cout << "Mech 11" << endl;

nodes.clear();
nodes.push_back(first_node);
nodes.push_back(second_node);
//cout << "div plane nodes size " << nodes.size() << endl;

return;
}*/
