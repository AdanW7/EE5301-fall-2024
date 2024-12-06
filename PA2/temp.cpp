// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)


#include<iostream>
#include<stdio.h>
#include<map>
#include<vector>
#include <string.h>
#include "suraj_parser.h"

using namespace std;

int star_count(vector <int> &hyperEdgeStartIndexes);
int circuit_W=0;
int circuit_H = 0;
double gamma(int numEdges);
void q_matrix_init(vector <vector <double>> &Q, vector <int> &cell_pin_array_l, vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l, 
vector <int> &hyper_weights_l, int num_cells_no_pads_l, int num_cells_and_pads_l);

void vect_init (vector <double> &x , vector<double> &y, vector <SPinLocation> &pin_locations_l);

void d_vect_init (vector <double> &dx, vector <double> &dy,vector <SPinLocation> & pin_locations_l,
 vector <int> &cell_pin_array_l, vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l, 
 vector <int> &hyper_weights_l, int &num_cells_no_pads_l, int &num_cells_and_pads_l);

void fill_d_vect (vector <double> &dx, vector <double> &dy, vector <SPinLocation> &pin_locations_l, 
int edge1, int edge2, vector<int> hyper_weights_l, int num_cells_no_pads_l, int num_cells_and_pads_l, int curr_Edge, int edge_Size);

void fill_q_matrix (vector <vector <double>> &Q, int edge1, int edge2, vector<int> hyper_weights_l, 
	int num_cells_no_pads_l, int num_cells_and_pads_l, int curr_Edge, int edge_Size);




	
int main(int argv, char *argc[])
{
	char inareFileName[100];
	char innetFileName[100];
	char inPadLocationFileName[100];

	if (argv!=2) {
		cout << "Please provide a circuit file name with no extension." << endl;
		return 1;
	}
        	     
    cout << "Reading circuit file " << argc[1] << endl;

	strcpy (inareFileName, argc[1]);
	strcat(inareFileName, ".are");
	strcpy(innetFileName,argc[1]);
	strcat(innetFileName,".net");
	strcpy(inPadLocationFileName,argc[1]);
	strcat(inPadLocationFileName,".kiaPad");

	int success = parseIbmFile(inareFileName, innetFileName, inPadLocationFileName);
	if (success == -1) {
		cout << "Error reading input file(s)" << endl;
		return 0;
	}

	printf("\nNumber of vertices,hyper = %d %d\n",numCellsAndPads,numhyper);
	

	
	// call function(s) dealing with creating the Q matrix, placement, etc.
	

	//rather than have non local extern data have local variables
	// not really needed but nice for debugging 
	int num_cell_pins_l = numCellPins; // # of terminals connected to the end points of (hyper) edges
	int num_hyper_l = numhyper; // number of edges and hyperedges
	int num_cells_and_pads_l = numCellsAndPads; // total number of movable cells (generall with names starting with a), and I/O pads (generally starting with p)
	int num_cells_no_pads_l = numCells_noPads; // total number of movable cells
	int num_pins_l = num_cells_and_pads_l - num_cells_no_pads_l; // total number of I/O pads or pins
	map <const char *, int, ltstr> node_name_to_node_num_map_l = nodeNameToNodeNum_map;


	// rather than deal with annoying pointer arithmetic i am going to use vectors rather than int* for the following data
	vector <int> cell_pin_array_l;
	vector <int> hyper_edge_index_to_first_entry_in_pin_array_l;
	vector <int> hyper_weights_l;
	vector <int> vertex_size_l;
	vector <SPinLocation> pin_locations_l;


	vector <vector <double>> Q;
	vector <double> dx;
	vector <double> dy; 
	vector <double> x;
	vector <double> y;

	cell_pin_array_l.resize(num_cell_pins_l);
	hyper_edge_index_to_first_entry_in_pin_array_l.resize(num_hyper_l+1);
	hyper_weights_l.resize(num_hyper_l);
	vertex_size_l.resize(num_cells_and_pads_l);
	pin_locations_l.resize(num_pins_l);

	for (int i = 0; i <= num_hyper_l; i++) { // might need to be <num_hyper_l
		hyper_edge_index_to_first_entry_in_pin_array_l[i] = hEdge_idxToFirstEntryInPinArray[i];
	}
	int stars = star_count(hyper_edge_index_to_first_entry_in_pin_array_l);
	cout << "star count= " << stars << endl;


	Q.resize(num_cells_no_pads_l + stars);
	for (int i = 0; i < num_cells_no_pads_l + stars; i++) {
		Q[i].resize(num_cells_no_pads_l + stars);
	}

	dx.resize(num_cells_no_pads_l + stars);
	dy.resize(num_cells_no_pads_l + stars);
	x.resize(num_cells_no_pads_l + stars);
	y.resize(num_cells_no_pads_l + stars);

	//fill local data copies with extern data from parser

	for (int i = 0; i < num_cell_pins_l; i++) {
		cell_pin_array_l[i] = cellPinArray[i]; 
	}
	

	for (int i = 0; i < num_hyper_l; i++) {
		hyper_weights_l[i] = hyperwts[i];
	}

	for (int i = 0; i < num_cells_and_pads_l; i++) {
		vertex_size_l[i] = vertexSize[i];
	}

	for (int i = 0; i < num_pins_l; i++) {
		pin_locations_l[i] = pinLocations[i];
	} 
	

	q_matrix_init(Q, cell_pin_array_l, hyper_edge_index_to_first_entry_in_pin_array_l, hyper_weights_l, 
	num_cells_no_pads_l, num_cells_and_pads_l);

	d_vect_init (dx, dy, pin_locations_l, cell_pin_array_l, hyper_edge_index_to_first_entry_in_pin_array_l, 
	hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l);

	vect_init(x, y, pin_locations_l);



	cout << "Q Matrix:" << endl;
	for (const auto &row : Q) {
        for (const auto &element : row) {
            std::cout << element << " ";
        }
        std::cout << std::endl;
    }

	cout << "dx vector:" << endl;
	for (const auto &element : dx) {
        std::cout << element << " ";
    }
	std::cout << std::endl;

	cout << "dy vector" << endl;
	for (const auto &element : dy) {
        std::cout << element << " ";
    }
	std::cout << std::endl;

	cout << "x vector:" << endl;
	for (const auto &element : x) {
        std::cout << element << " ";
    }
	std::cout << std::endl;

	cout << "y vector" << endl;
	for (const auto &element : y) {
        std::cout << element << " ";
    }
	std::cout << std::endl;




    free(pinLocations);
	free(hEdge_idxToFirstEntryInPinArray);
	free(cellPinArray);
	free(hyperwts);
}


int star_count(vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l) {
	int star_count = 0;
	int last_edge_idx = 0;
	int current_edge_idx= 0;

	// to determine if we are we have a star subtract index of curr edge from last edge,
	// if the result is greater than 3 then we have a star and can inc our star count
	for (unsigned int idx = 0; idx < hyper_edge_index_to_first_entry_in_pin_array_l.size(); idx++) {
		current_edge_idx = hyper_edge_index_to_first_entry_in_pin_array_l[idx];
		if (current_edge_idx - last_edge_idx > 3) {
			star_count += 1;
		}
		last_edge_idx = current_edge_idx;
	}

	return star_count;
}



void vect_init (vector <double> &x , vector<double> &y, vector <SPinLocation> &pin_locations_l) {

	for (unsigned int pinNum = 0; pinNum < pin_locations_l.size(); pinNum++) {
		if (pin_locations_l[pinNum].x > circuit_W) {
			circuit_W = pin_locations_l[pinNum].x;
		}

		if (pin_locations_l[pinNum].y > circuit_H) {
			circuit_H = pin_locations_l[pinNum].y;
		}
	}

	//"You can assume all moveable cells are at (ChipWidth/2, ChipHeight/2)"
	for (unsigned int cellNum = 0; cellNum < x.size(); cellNum++) {
		x[cellNum] = circuit_W / 2.0;
		y[cellNum] = circuit_H / 2.0;
	}

}

double gamma(int numEdges) {
	//k-pin net
	double gamma;
	//if we have more than 3 edges we are a star 
	// otherwise we are working with a clique

	//we set γ = 1/(k-1)

	if (numEdges > 3) { // is a star == γW = k/(k-1)
		gamma = numEdges / (numEdges - 1.0);
	} else {
		// is a clique = γW = 1/(k-1) 
		gamma = 1.0 / (numEdges - 1.0);
	}

	return gamma;

}

void d_vect_init (vector <double> &dx, vector <double> &dy,vector <SPinLocation> & pin_locations_l,
 vector <int> &cell_pin_array_l, vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l, 
 vector <int> &hyper_weights_l, int &num_cells_no_pads_l, int &num_cells_and_pads_l){
	
	int curr_Star = 0;

	// loop throught every hyperedge
	for (unsigned int curr_Edge = 0; curr_Edge < hyper_edge_index_to_first_entry_in_pin_array_l.size()-1; curr_Edge++) {
		
		
		unsigned int startIndex = hyper_edge_index_to_first_entry_in_pin_array_l[curr_Edge];

		unsigned int edge_Size = hyper_edge_index_to_first_entry_in_pin_array_l[curr_Edge+1] - hyper_edge_index_to_first_entry_in_pin_array_l[curr_Edge];

		// if an edge exsists between only 2 nodes its a clique
		if (edge_Size == 2) {
			int edge_1 = cell_pin_array_l[startIndex]; //1
			int edge_2 = cell_pin_array_l[startIndex+1]; //2

			fill_d_vect(dx, dy, pin_locations_l, edge_1, edge_2, hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l, curr_Edge, edge_Size);
			
		}

		// if an edge exsists between only 3 nodes its a clique
		else if (edge_Size == 3) {
			int edge_1, edge_2;
			
			for (int side = 0; side < 3; side++) {
				if (side == 0) {
					edge_1 = cell_pin_array_l[startIndex]; //1
					edge_2 = cell_pin_array_l[startIndex+1]; //2
				} else if (side == 1) {
					edge_1 = cell_pin_array_l[startIndex+1]; //2
					edge_2 = cell_pin_array_l[startIndex+2]; //3
 				} else {
					edge_1 = cell_pin_array_l[startIndex]; //1
					edge_2 = cell_pin_array_l[startIndex+2]; //3
				}

				fill_d_vect(dx, dy, pin_locations_l, edge_1, edge_2, hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l, curr_Edge, edge_Size);
			}

		}

		// if an edge exsists between more than 3 nodes its a star
		else {
			int edge_1, edge_2;
			vector <int> star_node;
			star_node.resize(edge_Size);

			for (unsigned int node = 0; node < edge_Size; node++) {
				star_node[node] = cell_pin_array_l[startIndex + node];
			}

			edge_2 = num_cells_and_pads_l + curr_Star;
			curr_Star += 1;

			for (unsigned int node = 0; node < edge_Size; node++) {
				edge_1 = star_node[node];
				fill_d_vect(dx, dy, pin_locations_l, edge_1, edge_2, hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l, curr_Edge, edge_Size);
			}
		}
	}
}
void fill_d_vect (vector <double> &dx, vector <double> &dy, vector <SPinLocation> &pin_locations_l, 
int edge_1, int edge_2, vector<int> hyper_weights_l, int num_cells_no_pads_l, int num_cells_and_pads_l, int curr_Edge, int edge_Size){

	int numPins = num_cells_and_pads_l - num_cells_no_pads_l;

	// edge1 and edge2 are pins, but not star nodes ignore
	if (edge_1 >= num_cells_no_pads_l && edge_2 >= num_cells_no_pads_l && edge_1 < num_cells_and_pads_l && edge_2 < num_cells_and_pads_l);
	
	// edge1 is cell, edge2 is pin
	else if (edge_2 >= num_cells_no_pads_l && edge_2 < num_cells_and_pads_l) {
		// edge1 is star node, reinsert into index
		if (edge_1 >= num_cells_and_pads_l) 
			edge_1 -= numPins;
		dx[edge_1] -= hyper_weights_l[curr_Edge] * gamma(edge_Size) * pin_locations_l[edge_2 - numPins].x;
		dy[edge_1] -= hyper_weights_l[curr_Edge] * gamma(edge_Size) * pin_locations_l[edge_2 - numPins].y;
	} 
	
	// edge1 is pin, edge2 is cell
	else if (edge_1 >= num_cells_no_pads_l && edge_1 < num_cells_and_pads_l) {
		// edge 2 is star node, reinsert into index
		if (edge_2 >= num_cells_and_pads_l){ 
			edge_2 -= numPins;
		}

		dx[edge_2] -= hyper_weights_l[curr_Edge] * gamma(edge_Size) * pin_locations_l[edge_1 - numPins].x;
		dy[edge_2] -= hyper_weights_l[curr_Edge] * gamma(edge_Size) * pin_locations_l[edge_1 - numPins].y;
	} 
	else ;// edge between 2 cells, ignore

}

//Ф(x)= (1/2)*(x^t)*(Qx) + (dx^t) + constant
//Let qij be the entry in row i and column j of matrix Q. 
//From expression (1), the cost in the x-direction between two movable cells i and j is (1/2)*(Wij)*(xi^2+xj^2-2xixj). 
//The first and second terms contribute Wij to qii and qjj respectively. 
//The third term contributes - Wij to qij and qji. From expression (2), 
//the cost in the x-direction between a movable cell i and a fixed cell f is (1/2)*(Wif)*(xi^2+xf^2-2xixf). 
//The first term contributes Wif to qii. The third term contributes -Wiff to the vector de at row i and the second term 
//contributes to the constant part of equation (4). The objective function (4) is minimized by solving the system of 
// linear equations represented by: Qx+dx = 0.
void q_matrix_init(vector <vector <double>> &Q, vector <int> &cell_pin_array_l, vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l,
 vector <int> &hyper_weights_l, int num_cells_no_pads_l, int num_cells_and_pads_l){
		int curr_Star = 0;

	// iterate through each (hyper)edge
	for (unsigned int curr_Edge = 0; curr_Edge < hyper_edge_index_to_first_entry_in_pin_array_l.size()-1; curr_Edge++) {
		
		
		unsigned int startIndex = hyper_edge_index_to_first_entry_in_pin_array_l[curr_Edge];
		unsigned int edge_Size = hyper_edge_index_to_first_entry_in_pin_array_l[curr_Edge+1] - hyper_edge_index_to_first_entry_in_pin_array_l[curr_Edge];

		// edge between 2 nodes (clique)
		if (edge_Size == 2) {
			int edge_1 = cell_pin_array_l[startIndex];
			int edge_2 = cell_pin_array_l[startIndex+1];

			fill_q_matrix(Q, edge_1, edge_2, hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l, curr_Edge, edge_Size);
		}

		// edge between 3 nodes (clique)
		else if (edge_Size == 3) {
			int edge_1;
			int edge_2;
			
			for (int side = 0; side < 3; side++) {
				if (side == 0) {
					edge_1 = cell_pin_array_l[startIndex];
					edge_2 = cell_pin_array_l[startIndex+1]; 
				} else if (side == 1) {
					edge_1 = cell_pin_array_l[startIndex+1];
					edge_2 = cell_pin_array_l[startIndex+2];
 				} else {
					edge_1 = cell_pin_array_l[startIndex];
					edge_2 = cell_pin_array_l[startIndex+2];
				}

				fill_q_matrix(Q, edge_1, edge_2, hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l, curr_Edge, edge_Size);
			}

		}

		// edge between 4+ nodes (star)
		else {
			int edge_1;
			int edge_2;
			vector <int> star_node;
			star_node.resize(edge_Size);

			for (unsigned int node = 0; node < edge_Size; node++) {
				star_node[node] = cell_pin_array_l[startIndex + node];
			}

			edge_2 = num_cells_and_pads_l + curr_Star; // edge 2 is the current number star its going through
			curr_Star ++;

			for (unsigned int node = 0; node < edge_Size; node++) {
				edge_1 = star_node[node];
				fill_q_matrix(Q, edge_1, edge_2, hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l, curr_Edge, edge_Size);
			}



		}

	}

}


void fill_q_matrix (vector <vector <double>> &Q, int edge1, int edge2, vector<int> hyper_weights_l, 
 int num_cells_no_pads_l, int num_cells_and_pads_l, int curr_Edge, int edge_Size){
	int numPins = num_cells_and_pads_l - num_cells_no_pads_l;

	// edge1 and edge2 are pins, but not star nodes ignore				
	if (edge1 >= num_cells_no_pads_l && edge2 >= num_cells_no_pads_l && edge1 < num_cells_and_pads_l && edge2 < num_cells_and_pads_l) ;
	 else if (edge1 >= num_cells_no_pads_l && edge1 < num_cells_and_pads_l) {
		// edge1 is pin, edge2 is cell
		if (edge2 >= num_cells_and_pads_l) // edge 2 is star node, reinsert into index
			edge2 -= numPins;
		Q[edge2][edge2] += hyper_weights_l[curr_Edge] * gamma(edge_Size);
	} 
	
	else if (edge2 >= num_cells_no_pads_l && edge2 < num_cells_and_pads_l) {
		// edge1 is cell, edge2 is pin
		if (edge1 >= num_cells_and_pads_l) // edge1 is star node, reinsert into index
			edge1 -= numPins;
		Q[edge1][edge1] += hyper_weights_l[curr_Edge] * gamma(edge_Size);
	} 
	else {
		// edge between 2 cells
		if (edge1 >= num_cells_and_pads_l){
			edge1 -= numPins;
		}
		if (edge2 >= num_cells_and_pads_l){
			edge2 -= numPins;
		}


		Q[edge1][edge1] += hyper_weights_l[curr_Edge] * gamma(edge_Size);
		Q[edge2][edge2] += hyper_weights_l[curr_Edge] * gamma(edge_Size);
		Q[edge1][edge2] -= hyper_weights_l[curr_Edge] * gamma(edge_Size);
		Q[edge2][edge1] -= hyper_weights_l[curr_Edge] * gamma(edge_Size);
	}
}