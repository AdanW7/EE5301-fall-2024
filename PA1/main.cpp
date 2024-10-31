#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <iomanip>
#include <algorithm>
#include <limits>

#include "Circuit.hpp"
// #include "Gatedatabase.hpp"
#include "CircuitNode.hpp"

using namespace std;

void DFF_to_in_out(Circuit &ckt);
void create_fanout_list(Circuit &ckt);
double war_crime(CircuitNode &node, double t_in, int delay_or_slew);
vector<double> calculate_a_out(CircuitNode &node, Circuit &ckt);
vector<double> calculate_t_out(CircuitNode &node, Circuit &ckt);
void forward_traversal(Circuit &ckt);
void backward_traversal(Circuit &ckt);
vector<CircuitNode*> find_crit_path(Circuit &ckt,CircuitNode &node);
void write_ckt_taversal(Circuit &ckt);
void cell_delay_func(Circuit &ckt, int &out_node, int &in_node);


int main(int argc, char* argv[]) {

    // ./sta <lib_file> <ckt_file>, argc =3, argv[1]=lib_file, argv[2]=ckt_file

    //comment following lines for debug, uncomment after
    if (argc != 3) {
        cout << "Need 2 arguments: <lib_file> <ckt_file>" << endl;
        return -1;
    }


    string lib_file(argv[1]);
    string ckt_file(argv[2]);

    // string lib_file="./test/NLDM_lib_max2Inp";
                                                                            // results on my local machine
    // string ckt_file("./test/cleaned_iscas89_99_circuits/c17.isc");       //time ./sta = real    0m0.369s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/c1908_.isc");    //time ./sta = real    0m0.037s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/c2670.isc");     //time ./sta = real    0m0.049s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/c3540.isc");     //time ./sta = real    0m0.051s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/c5351.isc");     //time ./sta = real    0m0.080s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/c7552.isc");     //time ./sta = real    0m0.102s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/b15_1.isc");     //time ./sta = real    0m0.380s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/b18_1.isc");     //time ./sta = real    0m3.643s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/b19_1.isc");     //time ./sta = real    0m7.703s
    // string ckt_file("./test/cleaned_iscas89_99_circuits/b22.isc");       //time ./sta = real    0m1.038s
    

    // cout << "Library File: " << lib_file << endl;
    // cout << "Circuit File: " << ckt_file << endl;    

    Circuit ckt(ckt_file, lib_file);

    DFF_to_in_out(ckt); // change dff's to inputs and outputs
    create_fanout_list(ckt); // creates a fan out list and inits in/out degree.
    forward_traversal(ckt); // find circuit delay
    backward_traversal(ckt); // gets gate slacks , required time =1.1*circuit delay
    write_ckt_taversal(ckt);

    return 0;
}




// Handle flip-flops, labelled as DFF, by splitting them into a dummy input and output nodes.
void DFF_to_in_out(Circuit &ckt){
    for (int i =1; i< ckt.nodes_.size();i++){
        if(ckt.nodes_[i]!=NULL){
                std::string tmpstring= ckt.nodes_[i]->get_gate_type();
                if (tmpstring.compare("DFF")==0){
                    ckt.nodes_[i]->set_input_pad(true);
                    ckt.nodes_[i]->set_output_pad(false);

                    int tmp_node=0;
                    tmp_node=ckt.nodes_[i]->fanin_list_.front();
                    ckt.nodes_[tmp_node]->set_output_pad(true);

                    ckt.nodes_[i]->fanin_list_.clear();
                    // ckt.nodes_[i]->add_to_fanout_list(tmp_node); // ahhhhhhhhhhh 30 mins wasted kmn
                }
                if(ckt.nodes_[i]->is_input_pad() && !(ckt.nodes_[i]->is_output_pad()) /*&& (ckt.nodes_[i]->gate_type_.empty()) */ ){
                    ckt.nodes_[i]->gate_type_="INP";
                }
            }
        }    
}

// before we can complete a forward traversal we need to know the fan outs of every node
void create_fanout_list(Circuit &ckt){
    for (int i =1; i< ckt.nodes_.size();i++){
        if(ckt.nodes_[i]!=NULL){
            if(!(ckt.nodes_[i]->fanin_list_.empty())){
                for (int j=0; j< ckt.nodes_[i]->fanin_list_.size();j++){
                    NodeID tmpnode=0;
                    tmpnode=ckt.nodes_[i]->fanin_list_[j];
                    ckt.nodes_[i]->in_degree++;
                    ckt.nodes_[tmpnode]->out_degree++; 
                    ckt.nodes_[tmpnode]->add_to_fanout_list(i);
                }
            }
        }
    }
}

// the purpose of this god forsaken function 
// is to complete what apppendix 2 shows in the project description manual,
// to who ever reads this i crave pain and i curse your blood line
// if delay_or_slew==0 we want delays else we want slews
double war_crime(CircuitNode &node, double t_in, int delay_or_slew){
    double table[GATE_LUT_DIM][GATE_LUT_DIM];
    double Capacitive_load=node.output_load;
    double t_1, t_2, c_1, c_2,v_11,v_12,v_21,v_22; // for Vxy, x=T(1,2), and Y=C(1,2), as per the manual
	double value;
    int x;
    int y;

    if (delay_or_slew == 0) {
        for (int i = 0; i < GATE_LUT_DIM; ++i) {
            for (int j = 0; j < GATE_LUT_DIM; ++j) {
                table[i][j] = node.gate_info_->cell_delay[i][j];
            }
        }
    } 
    else{ // if (delay_or_slew == 1) {
        for (int i = 0; i < GATE_LUT_DIM; ++i) {
            for (int j = 0; j < GATE_LUT_DIM; ++j) {
                table[i][j] = node.gate_info_->output_slew[i][j];
            }
        }
    }


    if(t_in > node.gate_info_->delay_index_1[6]){
        x=6;
    }
    else if(t_in < node.gate_info_->delay_index_1[0]){
        x=1;
    }
    else{   //t1 ≤ t < t2, so as to not break array start x at index 1
            for(x = 1; x < 7; x++){
			if ((node.gate_info_->delay_index_1[x-1] <= t_in) && (t_in < node.gate_info_->delay_index_1[x] )){
				break;
			}
		}
    }


    if(Capacitive_load > node.gate_info_->delay_index_2[6]){
        y=6;
    }
    else if(Capacitive_load < node.gate_info_->delay_index_2[0]){
        y=1;
    }
    else{//c1 ≤ c < c2, so as to not break array start y at index 1
        for(y = 1; y < 7; y++){
			if ( (node.gate_info_->delay_index_2[y-1] <= Capacitive_load) && (Capacitive_load < node.gate_info_->delay_index_2[y]) ){
				break;
			}
        }
    }



    t_1 = node.gate_info_->delay_index_1[x-1];
	t_2 = node.gate_info_->delay_index_1[x];
	c_1 = node.gate_info_->delay_index_2[y-1];
	c_2 = node.gate_info_->delay_index_2[y];
	v_11 = table[x-1][y-1];
	v_12 = table[x-1][y];
	v_21 = table[x][y-1];
	v_22 = table[x][y];

    value =( ( 
    ((v_11)*(c_2 - Capacitive_load)*(t_2 - t_in)) + 
    ((v_12)*(Capacitive_load - c_1)*(t_2 - t_in)) + 
    ((v_21)*(c_2 - Capacitive_load)*(t_in - t_1)) + 
    ((v_22)*(Capacitive_load - c_1)*(t_in - t_1)) ) 
    / ( (c_2 - c_1) * (t_2 - t_1) ) );
    return value;
}

// for every node in the circuit, a_out = max(ai +di), 
//where ai is the ith arrival time in vector of input arrival times
// and di 
vector<double> calculate_a_out(CircuitNode &node, Circuit &ckt){
    std::vector<double> temp_output_arrival_time;
	int  num_inputs=0;

    if (!node.is_input_pad()){
		num_inputs = node.fanin_list_.size();
	}
    double multiplier = 1.0;
    if (num_inputs > 1.0){
		multiplier = num_inputs/2;
	}

    vector<double> tmp_delay;
    NodeID i = 0;
    if(node.is_input_pad()){
        tmp_delay.resize(1);
    }
    else{
        tmp_delay.resize(num_inputs);
    }
    
    for (auto input_node: node.fanin_list_ ){
        double tmp_cell_delay= multiplier * war_crime(node, ckt.nodes_[input_node]->t_out,0); // we want delay not slew
        double temp_arrival = node.input_arrival_time[i] + tmp_cell_delay; 
        node.cell_delay=tmp_cell_delay;
        tmp_delay[i]=tmp_cell_delay;

        

        temp_output_arrival_time.push_back(temp_arrival); // determin what the cell delay for every input Tin is and add them to 
        i++;                                             // vector of output arrival times which we will later use in forward traversal
    }


    ckt.nodes_[node.node_id_]->output_arrival_time=temp_output_arrival_time;
    double arrival=*min_element(node.output_arrival_time.begin(), node.output_arrival_time.end());
    for ( int j=0; j<i ; j++){
        if(arrival == ((node.input_arrival_time[j])+(tmp_delay[j]))){
            if (!(node.input_arrival_time[j]==0)){
                node.cell_delay=tmp_delay[j];
            }
        }
        
    }
    return temp_output_arrival_time;
}




vector<double> calculate_t_out(CircuitNode &node, Circuit &ckt){
    std::vector<double> temp_output_slew;
	int  num_inputs=0;

    if (!node.is_input_pad()){
		num_inputs = node.fanin_list_.size();
	}
    double multiplier = 1.0;
    if (num_inputs > 1.0){
		multiplier = num_inputs/2;
	}

    for (auto input_node: node.fanin_list_ ){
        double slew= multiplier * war_crime(node, ckt.nodes_[input_node]->t_out,1); 
        temp_output_slew.push_back(slew);                                   
    }

    return temp_output_slew;

}


//For finding circuit delay, first perform the forward traversal of the circuit. At each gate, using the LUT,
// find the appropriate delays (and slews) of each path from its inputs to its output. The details of obtaining
// values corresponding to particular input slew and load cap for a path is provided in Appendix 2.
// Next find the arrival time at the output of this gate by the “max” function as shown in Fig. 1. Repeat this
// until you reach the primary outputs. The maximum among the arrival times at all the primary outputs is
// the circuit delay. 
void forward_traversal(Circuit &ckt){
    ckt.circuit_delay=0.0;
    queue <CircuitNode*> node_queue;

    for (auto node: ckt.nodes_){
        if(node!=NULL){
            node_queue.push(node);
        }
	}
    while(!node_queue.empty()){
        CircuitNode* temp_node = node_queue.front();

        if (temp_node->in_degree == 0) {

            if (temp_node->is_input_pad()) {
                //The arrival times at each primary input is 0, and input slew is 2 picoseconds (ps). 
                temp_node->input_arrival_time.push_back(0.0);
                temp_node->input_slew.push_back(0.002); // since lib file is all in nano seconds match it and fix at end of funciton

                for (auto node: temp_node->fanout_list_){
                    temp_node->output_load+=ckt.nodes_[node]->gate_info_->capacitance;
				}

                temp_node->output_arrival_time.push_back(0.0);
                temp_node->output_slew.push_back(0.002);
            }

            else{ // if gate isn't an input
                
                for (auto node: temp_node->fanin_list_){
                    // input slew and input arrival time are vectors that should contain all of the t_out & a_out of the input nodes
                    // we will use this data to calculate temp_nodes output arrival time
                    temp_node->input_slew.push_back(ckt.nodes_[node]->t_out);
                    temp_node->input_arrival_time.push_back(ckt.nodes_[node]->a_out); 
				}

                if(temp_node->is_output_pad()){
                    //The load capacitance of the final stage of gates 
                    //(which are simply connected to the primaryoutputs)
                    // is equal to four times the capacitance of an inverter from the liberty file. 
                    temp_node->output_load= 4.0 * ckt.gate_db_.gate_info_lut_["INV"]->capacitance;
                }
                else{
                    for (auto node: temp_node->fanout_list_){
                        temp_node->output_load+=ckt.nodes_[node]->gate_info_->capacitance;
				    }
                }
                temp_node->output_arrival_time=calculate_a_out(*temp_node, ckt);
                temp_node->output_slew= calculate_t_out(*temp_node, ckt);  
            }
            //a_out and t_out are the maximums for respective arrivals times and slews for all of the inputs to the gate
            temp_node->a_out=*max_element(temp_node->output_arrival_time.begin(), temp_node->output_arrival_time.end());
            temp_node->t_out=*max_element(temp_node->output_slew.begin(), temp_node->output_slew.end());

            if ((temp_node->is_output_pad()) &&((temp_node->a_out* 1000)> (ckt.circuit_delay)) ){
                ckt.circuit_delay= (temp_node->a_out) *(1000);
            }
            for(auto node: temp_node->fanout_list_){// since we have complete forward traversal of this gate update its fan out's in_degree
				ckt.nodes_[node]->in_degree--;
			}
            node_queue.pop();
        }

        else{
			node_queue.pop();
			node_queue.push(temp_node);
		}

    }
}


// this could be 
void cell_delay_func(Circuit &ckt, int &out_node, int &in_node){
        if(ckt.nodes_[out_node]!=NULL){
            double multiplier=0;
            if(ckt.nodes_[out_node]->is_input_pad()){
                multiplier=1;
            }
            else{
                multiplier=(ckt.nodes_[out_node]->fanin_list_.size())/2;
            }
            double tmp_delay= multiplier * war_crime(*ckt.nodes_[out_node], ckt.nodes_[in_node]->t_out,0); 
            ckt.nodes_[out_node]->cell_delay=tmp_delay;
	}
}



// For finding the slack at each gate, you need to perform a backward traversal of the circuit. Find the
// required arrival time at the output of each gate. The difference of the required arrival time and the actual
// arrival time (obtained from the forward traversal) is the slack at this gate. 
void backward_traversal(Circuit &ckt){

    //The required arrival time is 1.1 times the total circuit delay, 
    // and is the same at each primary output of the circuit.
    double required_time;
    required_time = 1.1 *ckt.circuit_delay;
    queue <CircuitNode*> node_queue;

    for (auto node: ckt.nodes_){
        if(node!=NULL){
            node_queue.push(node);
        }
	}
    while(!node_queue.empty()){
        CircuitNode* temp_node = node_queue.front();
        node_queue.pop();
        double tmp_required=required_time;

        if(temp_node->out_degree == 0){
            if(!temp_node->is_output_pad()){
                for(auto node: temp_node->fanout_list_){
                    //where node is an ouput node of temp_node
                    //tmp = arrival time of temp node output = minimum node's input arrival time = node's required arrival time - cell delay
                    cell_delay_func(ckt, node,temp_node->node_id_);
                    double tmp=((ckt.nodes_[node]->required_arrival_time) - (ckt.nodes_[node]->cell_delay)*1000);
                    if (tmp <= tmp_required){
						tmp_required = tmp;
					}
			    }
            }
            else{
                tmp_required=required_time;
            }
            temp_node->required_arrival_time=tmp_required;
            temp_node->gate_slack = abs((temp_node->required_arrival_time)- (temp_node->a_out)*1000); // double check this line 
            for (auto node: temp_node->fanin_list_){
                ckt.nodes_[node]->out_degree--;
			}
        }

        else{
            node_queue.push(temp_node);
        }
		
    }


}



// To find the critical path of the circuit based on the slack values - start with the primary output with the
// minimum slack, and traverse backwards selecting the gates connected to this output with minimum slack,
// till you reach the primary input. If more than one primary output have the same value of the slack that is
// minimum, returns the first critical path we found even if multiple are equal.
vector<CircuitNode*> find_crit_path(Circuit &ckt){
	vector<CircuitNode*> path;
	CircuitNode* temp_node;
    double temp_min = numeric_limits<double>::infinity();

    for (auto node: ckt.nodes_){ // find lowest slack output (litterally all should be equal?)
        if(node != NULL){
            if((node->is_output_pad()) && (node->gate_slack < temp_min) ){
			temp_min = node->gate_slack;
			temp_node = node;
		    }
        }
	} 

    path.push_back(temp_node);
    while(!temp_node->is_input_pad()){
        temp_min = numeric_limits<double>::infinity();
        for(auto node: temp_node->fanin_list_){
            if(ckt.nodes_[node]->gate_slack < temp_min){
                temp_min = ckt.nodes_[node]->gate_slack;
                temp_node = ckt.nodes_[node];
            }
		}
		path.push_back(temp_node);
    }
    return path;
}


// this function puts all of the above helper functions together 
// and displays our circuit delay, every single gate and its slack 
// and finally the critical path through the circuit
void write_ckt_taversal(Circuit &ckt){
    std::string result_file = "ckt_traversal.txt";
    std::ofstream outfile(result_file);
    if (!outfile.is_open())
    {
        std::cout << "Couldn't open the result file " << result_file << std::endl;
        return ;
    }

    if (outfile.good())
    {
        // std::cout << "able to send data to output " << result_file << std::endl;

        std::streambuf* coutbuf = std::cout.rdbuf(); // Save the original cout buffer
        std::cout.rdbuf(outfile.rdbuf()); // Redirect cout to the file
        
        cout << std::fixed << setprecision (2);
    
        std::cout<<"Circuit delay: "<< ckt.circuit_delay <<" ps\n"<< std::endl;
        std::cout<< "Gate slacks:"<< std::endl;
        
        for (int i =0; i< ckt.nodes_.size();i++){
            if(ckt.nodes_[i]!=NULL){
                ckt.print_node_info(i);
            }
            
        }

        // the following code section is used to display the critical path
        vector<CircuitNode*> crit_path =find_crit_path(ckt);
        std::cout << "\nCritical path:"<< endl;
        string critical_path_string;
        if(!crit_path.empty()){
            // for (int i=0; i < crit_path.size(); i++){    // prints in descending order
            for (int i=crit_path.size() -1; i >=0 ; i--){       // prints in ascending  
                // if (i!=0){ 
                if (i!=crit_path.size() -1){ 
                    critical_path_string.append(", ");
                }

                // cout<< crit_path[i]->gate_type_<< "-n"<< crit_path[i]->node_id_ ;
                critical_path_string.append(crit_path[i]->gate_type_ );
                critical_path_string.append( "-n");
                critical_path_string.append(to_string( crit_path[i]->node_id_));
            }
            cout << critical_path_string; // make critical path one string to cout to output file
        }

        std::cout.rdbuf(coutbuf); // Restore the original buffer
        // std::cout << "sent all data to output file " << result_file << std::endl;
    }
    outfile.close();

}

