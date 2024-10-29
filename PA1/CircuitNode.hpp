#ifndef CIRCUITNODE_HPP
#define CIRCUITNODE_HPP

#include <string>
#include <list>
#include "GateDatabase.hpp"

typedef int NodeID;

class CircuitNode {
    private:
        

    public:
        // because im bad set all data to public and ignore setters and getters

        NodeID node_id_;
        bool input_pad_;
        bool output_pad_;
        std::string gate_type_;
        const GateInfo* gate_info_;
        std::vector<NodeID> fanin_list_;
        std::vector<NodeID> fanout_list_;

        double gate_slack; // use backwards traversal to aquire for each gate

        std::vector<double> input_arrival_time; 
        std::vector<double> output_arrival_time; // using the input arrival time vector we determin the ouptput arrival time for each input
        std::vector<double> input_slew; 
        std::vector<double> output_slew; // using the input slew vector we determin the ouptput slew for each input
        double output_load; 
        int in_degree;
        int out_degree;


        double t_out; // max ouput slew we get from ouput slew vector
        double a_out; // max ouput arrival we get from ouput arrival vector
        double required_arrival_time;
        double cell_delay;

        

    CircuitNode() 
    : node_id_(-1),               // Default ID indicating an unassigned node
      input_pad_(false),           // No input pad by default
      output_pad_(false),          // No output pad by default
      gate_type_(""),              // Empty gate type as default
      gate_info_(nullptr),         // No gate information initially
      fanin_list_(),               // Empty vector of fanin nodes
      fanout_list_(),              // Empty vector of fanout nodes
      gate_slack(0),               // will be set to a_out - required_arrival time
      input_arrival_time(),       
      output_arrival_time(),      
      input_slew(),             
      output_slew(),            
      output_load(0.0),            // Initial output load set to zero
      in_degree(0),                // Initial in-degree set to zero
      out_degree(0),                // Initial out-degree set to zero
      t_out(0.0),
      a_out(0.0), // this is the actual arrival time of output
      required_arrival_time(0.0) ,
      cell_delay(0.0)
{}
        
        void set_node_id(const NodeID& node_id);
        void set_input_pad(const bool& input_pad);
        void set_output_pad(const bool& output_pad);
        void set_gate_type(const std::string& gate_type);
        void set_gate_info(const GateInfo* gate_info);
        void add_to_fanin_list(const NodeID& node_id);
        void add_to_fanout_list(const NodeID& node_id);
        
        const NodeID& get_node_id() const;
        const bool& is_input_pad() const;
        const bool& is_output_pad() const;
        const std::string& get_gate_type() const;
        const GateInfo* get_gate_info() const;
        const std::vector<NodeID>& get_fanin_list() const;
};

#endif //CIRCUITNODE_HPP
