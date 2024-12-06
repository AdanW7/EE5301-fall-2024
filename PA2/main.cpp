// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)


#include<iostream>
#include<stdio.h>
#include<map>
#include<vector>
#include <string.h>
#include "cmath"
#include "suraj_parser.h"
#include "assert.h"
#include "fstream"
#include "sstream"
using namespace std;

int star_count(vector <int> &hyperEdgeStartIndexes);
int circuit_W=0;
int circuit_H = 0;
double gamma(int numEdges);
double computeWireLength(const vector<vector<double>> &adjacencyMatrix, const vector<double> &xCoords, const vector<double> &yCoords);


void vect_init (vector <double> &x , vector<double> &y, vector <SPinLocation> &pin_locations_l);
void d_vect_init (vector <double> &dx, vector <double> &dy,vector <SPinLocation> & pin_locations_l,
	vector <int> &cell_pin_array_l, vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l, 
	vector <int> &hyper_weights_l, int &num_cells_no_pads_l, int &num_cells_and_pads_l);

void fill_d_vect (vector <double> &dx, vector <double> &dy, vector <SPinLocation> &pin_locations_l, 
	int edge1, int edge2, vector<int> hyper_weights_l, int num_cells_no_pads_l, 
	int num_cells_and_pads_l, int curr_Edge, int edge_Size);

void q_matrix_init(vector <vector <double>> &Q, vector <int> &cell_pin_array_l, 
	vector <int> &hyper_edge_index_to_first_entry_in_pin_array_l,vector <int> &hyper_weights_l, 
	int num_cells_no_pads_l, int num_cells_and_pads_l);

void fill_q_matrix (vector <vector <double>> &Q, int edge1, int edge2, vector<int> hyper_weights_l, 
	int num_cells_no_pads_l, int num_cells_and_pads_l, int curr_Edge, int edge_Size);


std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector);
double computeDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2);
std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2);
std::vector<double> subtractVectors(const std::vector<double>& vector1, const std::vector<double>& vector2);
std::vector<double> scaleVector(const std::vector<double>& vector, double scalar);
std::vector<double> solveConjugateGradient(const std::vector<std::vector<double>>& matrixA, const std::vector<double>& vectorRHS, 
	double tolerance, int maxIterations );

vector<vector<vector<int>>> initializeBinMatrix(vector<double> &xCoords, vector<double> &yCoords, int width, int height, int gridSize);
void spreadCellsHorizontally(vector<double> &spreadX, vector<double> &xCoords, vector<double> &yCoords, int width, int height, int gridSize);
void spreadCellsVertically(vector<double> &spreadY, vector<double> &xCoords, vector<double> &yCoords, int width, int height, int gridSize);

int load_CSV(vector<double> &xCoords, vector<double> &yCoords, string &fileName);
void write_CSV(const string &fileName, const vector<double> &xCoords, const vector<double> &yCoords, const vector<SPinLocation> &pinLocations, int numCells);

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
		fcloseall();
		return 0;
	}

	printf("\nNumber of vertices,hyper = %d %d\n",numCellsAndPads,numhyper);
	
	string circuit_name=(string) argc[1];

	string pre_spread_ = "pre_spread_";
	pre_spread_.append(circuit_name);
	pre_spread_.append(".csv");

	string post_spread_ = "post_spread_";
	post_spread_.append(circuit_name);
	post_spread_.append(".csv");
	
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
	
	free(pinLocations);
	free(hEdge_idxToFirstEntryInPinArray);
	free(cellPinArray);
	free(hyperwts);



	q_matrix_init(Q, cell_pin_array_l, hyper_edge_index_to_first_entry_in_pin_array_l, hyper_weights_l, 
	num_cells_no_pads_l, num_cells_and_pads_l);

	d_vect_init (dx, dy, pin_locations_l, cell_pin_array_l, hyper_edge_index_to_first_entry_in_pin_array_l, 
	hyper_weights_l, num_cells_no_pads_l, num_cells_and_pads_l);

	vect_init(x, y, pin_locations_l);



	if (!load_CSV(x,y,pre_spread_)){
		x = solveConjugateGradient(Q, dx, 1e-6,1000);
		y = solveConjugateGradient(Q, dy, 1e-6,1000);
		write_CSV(pre_spread_, x, y, pin_locations_l, num_cells_no_pads_l);
	}

	// cout << "star count= " << stars << endl;
	// cout << "Q Matrix:" << endl;
	// for (const auto &row : Q) {
    //     for (const auto &element : row) {
    //         std::cout << element << " ";
    //     }
    //     std::cout << std::endl;
    // }

	// cout << "dx vector:" << endl;
	// for (const auto &element : dx) {
    //     std::cout << element << " ";
    // }
	// std::cout << std::endl;

	// cout << "dy vector:" << endl;
	// for (const auto &element : dy) {
    //     std::cout << element << " ";
    // }
	// std::cout << std::endl;

	// cout << "x vector:" << endl;
	// for (const auto &element : x) {
    //     std::cout << element << " ";
    // }
	// std::cout << std::endl;

	// cout << "y vector:" << endl;
	// for (const auto &element : y) {
    //     std::cout << element << " ";
    // }
	// std::cout << std::endl;

	double pre_spread_wire_length= computeWireLength(Q, x, y);

	cout << "Pre_spreading W_L: " << pre_spread_wire_length << endl;

	vector <double> x_s;
	vector <double> y_s;
	x_s.resize(num_cells_no_pads_l + stars);
	y_s.resize(num_cells_no_pads_l + stars);
	spreadCellsHorizontally(x_s, x, y, circuit_W, circuit_H, 5);
	spreadCellsVertically(	y_s, x, y, circuit_W, circuit_H, 5);

	double post_spread_wire_length= computeWireLength(Q, x_s, y_s);
	cout << "Post_spreading W_L: " << post_spread_wire_length << endl;
	write_CSV(post_spread_, x_s, y_s, pin_locations_l, num_cells_no_pads_l);

    fcloseall();
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

double computeWireLength(const vector<vector<double>> &adjacencyMatrix, const vector<double> &xCoords, const vector<double> &yCoords) {
    double totalLength = 0.0;
    int numNodes = adjacencyMatrix.size();

    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            double weight = adjacencyMatrix[i][j];

            if (weight < 0) {
                double deltaX = xCoords[i] - xCoords[j];
                double deltaY = yCoords[i] - yCoords[j];
                totalLength += -weight * (pow(deltaX, 2.0) + pow(deltaY, 2.0));
            }
        }
    }

    return sqrt(totalLength);
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






// this entire code section was written with chat gpt
//  I used the following prompt with gpt 4-o
//  write c++ code to solve the upcomming matrix calculations using conjugate gradient method: 
//  Q * x + dx = 0; where Q is a n x n matrix, x and dx are of length n. 
//  the values of Q & dx are already known, find x 
// Function to perform matrix-vector multiplication O(vec.size()^2)
std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector) {
    int size = vector.size();
    std::vector<double> result(size, 0.0);
    for (int row = 0; row < size; ++row) {
        for (int col = 0; col < size; ++col) {
            result[row] += matrix[row][col] * vector[col];
        }
    }
    return result;
}

// Function to compute the dot product of two vectors (O(vec1.size()))
double computeDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    double result = 0.0;
    for (size_t index = 0; index < vector1.size(); ++index) {
        result += vector1[index] * vector2[index];
    }
    return result;
}

// Function to add two vectors (O(vec1.size()))
std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    int size = vector1.size();
    std::vector<double> result(size);
    for (int index = 0; index < size; ++index) {
        result[index] = vector1[index] + vector2[index];
    }
    return result;
}

// Function to subtract two vectors (O(vec1.size()))
std::vector<double> subtractVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    int size = vector1.size();
    std::vector<double> result(size);
    for (int index = 0; index < size; ++index) {
        result[index] = vector1[index] - vector2[index];
    }
    return result;
}

// Function to scale a vector by a scalar O(vec.size())
std::vector<double> scaleVector(const std::vector<double>& vector, double scalar) {
    int size = vector.size();
    std::vector<double> result(size);
    for (int index = 0; index < size; ++index) {
        result[index] = vector[index] * scalar;
    }
    return result;
}

// Conjugate Gradient method to solve Q * x + dx = 0
// Conjugate Gradient method to solve A * solution + rhs = 0
std::vector<double> solveConjugateGradient(const std::vector<std::vector<double>>& matrixA, const std::vector<double>& vectorRHS, double tolerance, int maxIterations ) {
    int size = vectorRHS.size();
    std::vector<double> solution(size, 0.0); // Initialize solution vector
    std::vector<double> residual = scaleVector(vectorRHS, -1.0); // Compute initial residual
    std::vector<double> direction = residual;
    double residualSquaredOld = computeDotProduct(residual, residual);

    for (int iteration = 0; iteration < maxIterations; ++iteration) {

        std::vector<double> matrixDirection = multiplyMatrixVector(matrixA, direction);
        double stepSize = residualSquaredOld / computeDotProduct(direction, matrixDirection);
        solution = addVectors(solution, scaleVector(direction, stepSize));
        residual = subtractVectors(residual, scaleVector(matrixDirection, stepSize));
        double residualSquaredNew = computeDotProduct(residual, residual);

        if (sqrt(residualSquaredNew) < tolerance) {
            break;
        }

        direction = addVectors(residual, scaleVector(direction, residualSquaredNew / residualSquaredOld));
        residualSquaredOld = residualSquaredNew;
    }

    return solution;
}




//i used chat gpt 4-0 giving int incremental information from the fast place paper to generate the following code
vector<vector<vector<int>>> initializeBinMatrix(vector<double> &xCoords, vector<double> &yCoords, int width, int height, int gridSize) {
    int cellCount = xCoords.size();

    vector<vector<vector<int>>> binGrid;
    vector<double> horizontalBounds;
    vector<double> verticalBounds;
    horizontalBounds.resize(gridSize);
    verticalBounds.resize(gridSize);

    horizontalBounds[0] = 0;
    verticalBounds[0] = 0;
    for (int i = 1; i < gridSize; i++) {
        horizontalBounds[i] = horizontalBounds[i - 1] + (width / gridSize);
        verticalBounds[i] = verticalBounds[i - 1] + (height / gridSize);
    }

    binGrid.resize(gridSize);
    for (int i = 0; i < gridSize; i++) {
        binGrid[i].resize(gridSize);
    }

    for (int cellIdx = 0; cellIdx < cellCount; cellIdx++) {
        double xValue = xCoords[cellIdx];
        double yValue = yCoords[cellIdx];
        int xBin = 0;
        int yBin = 0;

        for (int i = 0; i < gridSize; i++) {
            if (xValue >= horizontalBounds[i] && xValue < horizontalBounds[i + 1]) {
                xBin = i;
            }
        }

        if (xValue >= horizontalBounds[gridSize - 1]) {
            xBin = gridSize - 1;
        }

        for (int i = 0; i < gridSize; i++) {
            if (yValue >= verticalBounds[i] && yValue < verticalBounds[i + 1]) {
                yBin = i;
            }
        }

        if (yValue >= verticalBounds[gridSize - 1]) {
            yBin = gridSize - 1;
        }

        if (cellIdx == 8419) {
            cout << xBin << " " << yBin << endl;
            cout << xCoords[8419] << endl;
            cout << yCoords[8419] << endl;
        }

        binGrid[yBin][xBin].push_back(cellIdx);
    }

    return binGrid;
}

void spreadCellsHorizontally(vector<double> &spreadX, vector<double> &xCoords, vector<double> &yCoords, int width, int height, int gridSize) {
    vector<vector<vector<int>>> binGrid = initializeBinMatrix(xCoords, yCoords, width, height, gridSize);
    vector<double> oldBoundaries;
    vector<double> newBoundaries;
    vector<double> cellCountInBins;
    double adjustmentFactor = 1.5;
    double smoothingFactor = 0.8;
    oldBoundaries.resize(gridSize + 1);
    newBoundaries.resize(gridSize + 1);
    cellCountInBins.resize(gridSize + 1);

    newBoundaries[0] = 0;
    oldBoundaries[gridSize] = width;
    newBoundaries[gridSize] = width;

    oldBoundaries[0] = 0;
    for (int i = 1; i < gridSize; i++) {
        oldBoundaries[i] = oldBoundaries[i - 1] + (width / (double)gridSize);
    }

    for (int row = 0; row < gridSize; row++) {
        for (int col = 0; col < gridSize; col++) {
            cellCountInBins[col] = binGrid[row][col].size();
        }

        for (int col = 1; col < gridSize; col++) {
            double numeratorLeft = oldBoundaries[col - 1] * (cellCountInBins[col] + adjustmentFactor);
            double numeratorRight = oldBoundaries[col + 1] * (cellCountInBins[col - 1] + adjustmentFactor);
            double denominator = cellCountInBins[col - 1] + cellCountInBins[col] + 2.0 * adjustmentFactor;
            newBoundaries[col] = (numeratorLeft + numeratorRight) / denominator;

            assert(newBoundaries[col] <= width);
        }

        for (int col = 1; col < gridSize - 1; col++) {
            vector<int> cellList = binGrid[row][col];

            for (unsigned int cellIdx = 0; cellIdx < cellList.size(); cellIdx++) {
                int cell = cellList[cellIdx];

                double numeratorLeft = newBoundaries[col + 1] * (xCoords[cell] - oldBoundaries[col]);
                double numeratorRight = newBoundaries[col] * (oldBoundaries[col + 1] - xCoords[cell]);
                double denominator = oldBoundaries[col + 1] - oldBoundaries[col];
                spreadX[cell] = (numeratorLeft + numeratorRight) / denominator;

                double originalX = xCoords[cell];
                double adjustedX = spreadX[cell];
                spreadX[cell] = originalX + smoothingFactor * (adjustedX - originalX);

                assert(spreadX[cell] <= width);
            }
        }
    }
}

void spreadCellsVertically(vector<double> &spreadY, vector<double> &xCoords, vector<double> &yCoords, int width, int height, int gridSize) {
    vector<vector<vector<int>>> binGrid = initializeBinMatrix(xCoords, yCoords, width, height, gridSize);
    vector<double> oldBoundaries;
    vector<double> newBoundaries;
    vector<double> cellCountInBins;
    double adjustmentFactor = 1.5;
    double smoothingFactor = 0.8;
    oldBoundaries.resize(gridSize + 1);
    newBoundaries.resize(gridSize + 1);
    cellCountInBins.resize(gridSize + 1);

    newBoundaries[0] = 0;
    oldBoundaries[gridSize] = height;

    oldBoundaries[0] = 0;
    for (int i = 1; i < gridSize; i++) {
        oldBoundaries[i] = oldBoundaries[i - 1] + (height / gridSize);
    }

    for (int col = 0; col < gridSize; col++) {
        for (int row = 0; row < gridSize; row++) {
            cellCountInBins[row] = binGrid[row][col].size();
        }

        for (int row = 1; row < gridSize; row++) {
            double numeratorLeft = oldBoundaries[row - 1] * (cellCountInBins[row] + adjustmentFactor);
            double numeratorRight = oldBoundaries[row + 1] * (cellCountInBins[row - 1] + adjustmentFactor);
            double denominator = cellCountInBins[row - 1] + cellCountInBins[row] + 2.0 * adjustmentFactor;
            newBoundaries[row] = (numeratorLeft + numeratorRight) / denominator;

            assert(newBoundaries[row] <= height);
        }

        for (int row = 1; row < gridSize; row++) {
            vector<int> cellList = binGrid[row][col];

            for (unsigned int cellIdx = 0; cellIdx < cellList.size(); cellIdx++) {
                int cell = cellList[cellIdx];

                double numeratorLeft = newBoundaries[row + 1] * (yCoords[cell] - oldBoundaries[row]);
                double numeratorRight = newBoundaries[row] * (oldBoundaries[row + 1] - yCoords[cell]);
                double denominator = oldBoundaries[row + 1] - oldBoundaries[row];
                spreadY[cell] = (numeratorLeft + numeratorRight) / denominator;

                double originalY = yCoords[cell];
                double adjustedY = spreadY[cell];
                spreadY[cell] = originalY + smoothingFactor * (adjustedY - originalY);

                assert(spreadY[cell] <= height);
            }
        }
    }
}

int load_CSV(vector<double> &xCoords, vector<double> &yCoords, string &fileName) {
    ifstream inputFile (fileName);

    // Check if the file can be opened
    if (!inputFile.is_open()) {
        return 0; // Return 0 on error
    }

    string line, type;
    int positionIndex = 0;

    // Skip the header line
    getline(inputFile, line);

    // Read the file line by line
    while (getline(inputFile, line)) {
        istringstream lineStream(line);
        string token;
		if (token == "pin") {
            break; // Stop processing if type is "p"
        }
        // Read the type field
        getline(lineStream, type, ',');
        if (type == "pin") {
            break; // Stop processing if type is "p"
        }

        // Read and store the X coordinate
        getline(lineStream, token, ',');
        xCoords[positionIndex] = stod(token);

        // Read and store the Y coordinate
        getline(lineStream, token, ',');
        yCoords[positionIndex] = stod(token);

        positionIndex++;
    }

    return 1; // Return 1 on successful execution
}


void write_CSV(const string &fileName, const vector<double> &xCoords, const vector<double> &yCoords, const vector<SPinLocation> &pinLocations, int numCells) {
    ofstream outputFile(fileName);

    // Check if the file is open
    if (!outputFile.is_open()) {
        cerr <<"couldn't write results to " <<fileName << endl;
        return;
    }

    ostringstream contentStream;

    // Write the header
    contentStream << "type, x coord, y coord" << endl;

    // Write cell and star node positions
    int totalEntries = xCoords.size();
    for (int i = 0; i < totalEntries; i++) {
        if (i < numCells) {
            contentStream << "cell, " << xCoords[i] << ", " << yCoords[i] << endl;
        } else {
            contentStream << "star, " << xCoords[i] << ", " << yCoords[i] << endl;
        }
    }

    // Write pin positions
    for (const auto &pin : pinLocations) {
        contentStream << "pin, " << pin.x << ", " << pin.y << endl;
    }

    // Write to file and close
    outputFile << contentStream.str();
    outputFile.close();
}