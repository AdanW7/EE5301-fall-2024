#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>

#include "GateDatabase.hpp"

using namespace std;

GateDatabase::GateDatabase(const std::string& file_name) {
    // cout << "Parsing database file: " << file_name << endl;

    // Open file. Check if opened correctly
    // Always check for this to save a lot of pain if an error occurs
    ifstream ifs(file_name.c_str());
    if (!ifs.is_open()) {
        cout << "Error opening file " << file_name << endl;
        return;
    }

    string gate_name;

    size_t values_row_idx = 0;
    size_t values_col_idx = 0;

    // These flags prevent searches and also as a state machine 
    // so we know which step we are in when parsing line by line
    struct flag_found {
        bool cell;
        bool cell_delay;
        bool output_slew;
        bool values;
        bool capacitance;
        bool delay_index_1;
        bool delay_index_2;
        bool slew_index_1;
        bool slew_index_2;
    };

    // These flags are to ensure both cell_delay and output_slew are 
    // both parsed for each entry of cell. They can be parsed in any order
    struct flag_parsed {
        bool cell_delay;
        bool output_slew;
        bool capacitance;
        bool delay_index_1;
        bool delay_index_2;
        bool slew_index_1;
        bool slew_index_2;
       
    };

    flag_found found = {0};
    flag_parsed parsed = {0};

    GateInfo* gate_info = nullptr;

    const regex cell_regex("cell\\(([a-zA-Z0-9_]*)\\).*");
    smatch cell_regex_match;

    while (ifs.good()) {
        string line;
        getline(ifs, line);

        // Remove all whitespace to make regex and parsing simpler
        line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());

        if (!found.cell) {
            // Find cell(<gate_name>)
            regex_match(line, cell_regex_match, cell_regex);
            if (cell_regex_match.size() > 0) {
                gate_name = cell_regex_match[1];
                gate_info = new GateInfo;
                found.cell = true;
            }
        } 
        
        if (found.cell && !found.capacitance && !parsed.capacitance) {
            if (line.find("capacitance") != string::npos) {
                found.capacitance = true;
                gate_info->capacitance=stod(line.substr(line.find(":")+1, line.find(";")-line.find(":")-1));
                parsed.capacitance=true;
            }
        }


        if (found.cell && !found.cell_delay && !parsed.cell_delay) {
            if (line.find("cell_delay") != string::npos) {
                found.cell_delay = true;
            }
        }
        if (found.cell && found.cell_delay &&!found.delay_index_1 && !parsed.delay_index_1) {
            if (line.find("index_1") != string::npos) {
                found.delay_index_1 = true;

                string tmp_line=line.substr(line.find("\"")+1);
                // tmp_line=line.substr(0,tmp_line.find("\""));
                string tmp_line_2 = tmp_line;
                stringstream ss (tmp_line_2);
                double tmpval;
                char delim;
                for (int i =0; i <GATE_LUT_DIM; i++){
                    ss >> tmpval;
                    gate_info->delay_index_1[i]=tmpval;
                    ss >> delim;
                }

                parsed.delay_index_1=true;
            }
        }

        if (found.cell && found.cell_delay &&!found.delay_index_2 && !parsed.delay_index_2) {
            if (line.find("index_2") != string::npos) {
                found.delay_index_2 = true;


                string tmp_line=line.substr(line.find("\"")+1);
                // tmp_line=line.substr(0,tmp_line.find("\""));

                stringstream ss(tmp_line);
                double tmpval;
                char delim;
                for (int i =0; i <GATE_LUT_DIM; i++){
                    ss>>tmpval;
                    gate_info->delay_index_2[i]=tmpval;
                    ss>>delim;
                }

                parsed.delay_index_2=true;
            }
        }


        if (found.cell && !found.output_slew && !parsed.output_slew) {
            if (line.find("output_slew") != string::npos) {
                found.output_slew = true;
            }
        }

        if (found.cell && found.output_slew &&!found.slew_index_1 && !parsed.slew_index_1) {
            if (line.find("index_1") != string::npos) {
                found.slew_index_1 = true;

                string tmp_line=line.substr(line.find("\"")+1);
                // tmp_line=line.substr(0,tmp_line.find("\""));

                stringstream ss(tmp_line);
                double tmpval;
                char delim;
                for (int i =0; i <GATE_LUT_DIM; i++){
                    ss>>tmpval;
                    gate_info->slew_index_1[i]=tmpval;
                    ss>>delim;
                }


                parsed.slew_index_1=true;
            }
        }

        if (found.cell && found.output_slew &&!found.slew_index_2 && !parsed.slew_index_2) {
            if (line.find("index_2") != string::npos) {
                found.slew_index_2 = true;
                
                string tmp_line=line.substr(line.find("\"")+1);
                // tmp_line=line.substr(0,tmp_line.find("\""));

                stringstream ss(tmp_line);
                double tmpval;
                char delim;
                for (int i =0; i <GATE_LUT_DIM; i++){
                    ss>>tmpval;
                    gate_info->slew_index_2[i]=tmpval;
                    ss>>delim;
                }

                parsed.slew_index_2=true;
            }
        }


        if ((found.cell_delay || found.output_slew) && !found.values) {
            if (line.find("values(") != string::npos) {
                found.values = true;

                // Starting new table parsing
                values_row_idx = 0;
            }
        }




        if (found.values) {
            size_t start_idx = line.find_first_of("\"");
            size_t end_idx = line.find_last_of("\"");

            if (start_idx != string::npos) {
                // Starting new column parsing
                values_col_idx = 0;

                // Get string between the ""
                string s = line.substr(start_idx + 1, end_idx - start_idx);
                
                double value;
                char delim;
                stringstream delays_str(s);

                // Keep reading our doubles, unless we hit a eof or other errors
                while (delays_str >> value) {
                    if (found.cell_delay) {
                        gate_info->cell_delay[values_row_idx][values_col_idx] = value;
                    } else if (found.output_slew) {
                        gate_info->output_slew[values_row_idx][values_col_idx] = value;
                    }
                    values_col_idx++;

                    // Read in the delimiter (,)
                    delays_str >> delim;
                }
                values_row_idx++;

                // Check if we parsed our entire table
                if (values_row_idx >= GATE_LUT_DIM) {
                    // Mark the appropriate table as parsed
                    if (found.cell_delay) {
                        parsed.cell_delay = true;
                        found.cell_delay = false;
                    } else if (found.output_slew) {
                        parsed.output_slew = true;
                        found.output_slew = false;
                    }

                    // Reset values, so we parse the next table
                    found.values = false;

                    // Insert into our map if we have parsed both
                    if (parsed.cell_delay && parsed.output_slew) {
                        insert(gate_name, gate_info);
                        gate_info = nullptr;

                        // Done with parsing a single cell. Reset all flags
                        found = {0};
                        parsed = {0};
                    }
                }
            }
        }
    }
    ifs.close();
}

GateDatabase::~GateDatabase() {
    map<string, GateInfo*>::iterator itr;
    for (itr = gate_info_lut_.begin(); itr != gate_info_lut_.end(); ++itr) {
        // Delete our allocated GateInfo in the map
        delete itr->second;
        itr->second = nullptr;
    }
}

void GateDatabase::insert(const std::string& gate_name, GateInfo* gate_info) {
    gate_info_lut_.insert({gate_name, gate_info});
}

const GateInfo* GateDatabase::get_gate_info(const std::string& gate_name) {
    map<string, GateInfo*>::iterator itr;
    itr = gate_info_lut_.find(gate_name);
    
    // If we found it, return reference to the GateInfo as a const pointer
    if (itr != gate_info_lut_.end())
        return itr->second;
    
    return nullptr;
}

void GateDatabase::test() {
    for (const auto& key_value: gate_info_lut_) {
        cout << key_value.first << '\t' << key_value.second->cell_delay[GATE_LUT_DIM-1][GATE_LUT_DIM-1] 
             << '\t' << key_value.second->output_slew[GATE_LUT_DIM-1][GATE_LUT_DIM-1] << '\n';
    }
}
