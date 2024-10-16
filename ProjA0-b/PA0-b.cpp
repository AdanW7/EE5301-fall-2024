#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <array>

struct ParenCommaEq_is_space : std::ctype<char>
{
    ParenCommaEq_is_space() : std::ctype<char>(get_table()) {}
    static mask const *get_table()
    {
        static mask rc[table_size];
        rc['('] = std::ctype_base::space;
        rc[')'] = std::ctype_base::space;
        rc[','] = std::ctype_base::space;
        rc['='] = std::ctype_base::space;
        rc[' '] = std::ctype_base::space;
        rc['\t'] = std::ctype_base::space;
        return &rc[0];
    }
};

//shitty helper funciton that I will use to determine if we should skip forward in circuit library
bool contains_digit(const std::string& str) {
    // Use std::any_of to check if any character in the string is a digit
    return std::any_of(str.begin(), str.end(), ::isdigit);
}

struct Gate
{
    std::string type;
    std::vector<int> inputs;
    double slew_value; // Third row, second column [2][1] from the output_slew table

};

class CircuitParser
{

private:
    std::map<std::string, std::pair<std::array<std::array<double, 7>, 7>, std::array<std::array<double, 7>, 7> > >  gate_library;
    std::vector<Gate> circuit;

public:
    CircuitParser(const std::string &library_file, const std::string &circuit_file)

    {
        circuit.resize(1000); // Resize the vector to hold 1000 elements
        // as code is written we must always parse library before circuit,
        // this is because im bad and some of parsecircuit relies on parselib lmaooooooo
        parseLibraryFile(library_file);
        parseCircuitFile(circuit_file);
    }
    void parseLibraryFile(const std::string &library_file)
    {

        std::cout << "Parsing library file: " << library_file << std::endl;
        std::ifstream lib_file(library_file);
        std::string line;

        // Open file. Check if opened correctly
        // Always check for this to save a lot of pain if an error occurs
        if (!lib_file.is_open())
        {
            std::cout << "Error opening file " << library_file << std::endl;
            return;
        }

        // double slew_data[7][7];
        // double delay_data[7][7];
        std::array<std::array<double, 7>, 7> delay_data;
        std::array<std::array<double, 7>, 7> slew_data;
        std::string gate_name;
        while (lib_file.good())
        {

            std::string first_word;
            lib_file >> first_word;

            if (first_word.compare("cell") == 0)
            {
                std::string second_word;
                lib_file >> second_word;

                // The gate name has the format of "(<gate_name>)".
                // So find the index till the closing bracket.
                // Make sure to ignore the first index which has the opening bracket (hence the 1 and -1)
                size_t delim_pos = second_word.find(")");
                // std::string gate_name;
                gate_name = second_word.substr(1, delim_pos - 1);
                gate_library[gate_name];
            }
            else if (first_word.compare("cell_delay(Timing_7_7)") == 0)
            {
                // Read 3 lines that contain the rest of above match, index 1 and index 2
                std::string tmp;
                getline(lib_file, tmp);
                getline(lib_file, tmp);
                getline(lib_file, tmp);

                // From here on the next 7 lines will contain our delays
                for (size_t i = 0; i < 7; i++)
                {
                    getline(lib_file, tmp);

                    // The delays will be between " ". Find the opening ".
                    size_t start_delim_idx = tmp.find("\"");

                    // Find the closing ".
                    // The second argument is where we want to start our search
                    // Ignore the first match so we don't get the same index again
                    size_t end_delim_idx = tmp.find("\"", start_delim_idx + 1);

                    // The second arg in substr in no. of characters, not the ending index
                    std::string data_str = tmp.substr(start_delim_idx + 1, end_delim_idx - start_delim_idx - 1);

                    // Convert this remaining string to a stream so we can parse our data in doubles
                    std::istringstream data_stream(data_str);
                    for (size_t j = 0; j < 7; j++)
                    {
                        double delay;
                        char delim;
                        data_stream >> delay >> delim;
                        delay_data[i][j] = delay;
                        // gate_library[gate_name].first[i][j] = delay_data[i][j];
                    }
                }
            }
            else if (first_word.compare("output_slew(Timing_7_7)") == 0)
            {
                // Read 3 lines that contain the rest of above match, index 1 and index 2
                std::string tmp;
                getline(lib_file, tmp);
                getline(lib_file, tmp);
                getline(lib_file, tmp);

                // From here on the next 7 lines will contain our delays
                for (size_t i = 0; i < 7; i++)
                {
                    getline(lib_file, tmp);

                    // The delays will be between " ". Find the opening ".
                    size_t start_delim_idx = tmp.find("\"");

                    // Find the closing ".
                    // The second argument is where we want to start our search
                    // Ignore the first match so we don't get the same index again
                    size_t end_delim_idx = tmp.find("\"", start_delim_idx + 1);

                    // The second arg in substr in no. of characters, not the ending index
                    std::string data_str = tmp.substr(start_delim_idx + 1, end_delim_idx - start_delim_idx - 1);

                    // Convert this remaining string to a stream so we can parse our data in doubles
                    std::istringstream data_stream(data_str);
                    for (size_t j = 0; j < 7; j++)
                    {
                        double slew;
                        char delim;
                        data_stream >> slew >> delim;
                        slew_data[i][j] = slew;
                    }
                }
                // std::cout << gate_name << std::endl ; // debug line to state all the different gate types in our lib
                gate_library[gate_name] = std::make_pair(delay_data, slew_data);             
            }
        }
        std::cout << "we have parsed all of the data in "<< library_file << std::endl ; // debug line to state all the different gate types in our li
        lib_file.close();
    }

    void parseCircuitFile(const std::string &circuit_file)
    {
        std::ifstream cir_file(circuit_file);
        // cir_file.imbue(std::locale(std::cin.getloc(), new ParenCommaEq_is_space));
        std::string line;

        if (!cir_file.is_open())
        {
            std::cout << "Error opening file " << circuit_file << std::endl;
            return;
        }
        else{
            std::cout << "parsing circuit file " << circuit_file << std::endl;
        }
        while (cir_file.good())
        {
            while (std::getline(cir_file, line))
            {
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                std::istringstream iss(line);
                iss.imbue(std::locale(std::cin.getloc(), new ParenCommaEq_is_space));
                std::string first_word;

                iss>> first_word;

                // Skip comment lines
                if (first_word[0] == '#')
                    continue;

                // Handle INPUT and OUTPUT
                else if (first_word == "OUTPUT")
                {        
                    continue;
                }
                else if(first_word == "INPUT" ){
                    int node;
                    std::string tmp = line;
                    


                    size_t start_delim_idx = tmp.find("( ");
                    size_t end_delim_idx = tmp.find(" )", start_delim_idx);

                    std::string data_str = tmp.substr(start_delim_idx + 2, end_delim_idx - start_delim_idx - 2);

                    std::istringstream data_stream(data_str);
                    int tmp_num = std::stoi(data_str);
                    while (tmp_num >= circuit.size())
                        {
                            circuit.resize(circuit.size() + 1000);
                        }
                    circuit[tmp_num].inputs.push_back(tmp_num);
                    
                    continue;
                }
                else{
                    if (contains_digit(first_word)) {
                            // Parse gate and inputs
                        int output_gate_num=std::stoi(first_word);
                        std::string gate_type = "";
                        int temp_input=0;
                        
                        iss.imbue(std::locale(std::cin.getloc(), new ParenCommaEq_is_space));
                        iss >> gate_type;
                        while (output_gate_num >= circuit.size())
                        {
                            circuit.resize(circuit.size() + 1000);
                        }
                        circuit[output_gate_num].type.append(gate_type);

                        std::string token;

                        while (iss >> token) {
                            std::stringstream ss(token);
                            int num;
                            if (ss >> num) {
                                while (num >= circuit.size())
                                {
                                    circuit.resize(circuit.size() + 1000);
                                }
                                circuit[output_gate_num].inputs.push_back(num);
                            }
                        }
                        if (gate_library.count(gate_type) <= 0) {
                            std::cout << "this gate type " << gate_type << " doesn't exist in the library" << std::endl;
                            for (int i =0;i<7;i++){
                                for (int j=0; j<7; j++){
                                    gate_library[gate_type].first[i][j]=0;
                                    gate_library[gate_type].second[i][j]=0;
                                }
                            }
                            circuit[output_gate_num].slew_value=0;
                        }
                        else{
                            double temp_slew;
                            for (const auto& pair : gate_library) {
                                if (pair.first == gate_type) {
                                    // temp_slew=gate_library[gate_type].second.__elems_[2][1];
                                    temp_slew=gate_library[gate_type].second[2][1];
                                    circuit[output_gate_num].slew_value=temp_slew;
                                }
                            }
                        }
                    }
                }
            }
        }
        std::cout << "we have parsed all of " << circuit_file <<"\n" << std::endl; 
        cir_file.close();
    }

    void printGateInfo(int node) const
    {
        if(node >= circuit.size()) {
                std::cerr << "Node " << node << " does not exist in the provided circuit file.\n";
                return;
        }
        const Gate &gate = circuit[node];
        if (circuit[node].type.empty())
        {
            if(gate.inputs.front() == gate.inputs.back()){
            std::cout << node<< " INPUT" << std::endl;
            }
            
        }
        else{
            std::cout << node << " " << gate.type << " "
                << gate.inputs.front() << " " << gate.inputs.back() << " "
                << gate.slew_value << std::endl;
        }
        
    }
};

int main(int argc, char **argv)
{

    // this bit of code is from the example code in FileIO.cpp on canvas
    if (argc < 3)
    {
        std::cout << "I need atleast 3 paramaters, one library file, one circuit file, and atleast one node number" << std::endl;
        return -1;
    }

    std::string library_file = argv[1]; // uncomment after debug
    std::string circuit_file = argv[2]; // uncomment after debug
    
    
    // std::string library_file = "./NLDM_lib_max2Inp"; // for debug
    // std::string circuit_file = "./cleaned_iscas89_99_circuits/b22.isc"; // for debug

    CircuitParser parser(library_file, circuit_file);

    for (int i = 3; i < argc; ++i)
    {
        int node = std::stoi(argv[i]);
        parser.printGateInfo(node);
    }

    // std::string test_array[1] = {"8731"};
    //     for (int i = 0; i < 1; ++i)
    // {
    //     int node = stoi(test_array[i]);
    //     parser.printGateInfo(node);
    // }

    std::cout << std::endl;
    return 0;
}