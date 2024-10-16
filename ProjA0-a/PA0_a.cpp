// Author: Adan Wodzinski wodzi003@umn.edu
// Created: 9/4/24

#include "./PA0_a.h"


void parse_file(std::string input_file, std::string result_file, Gate *my_gate)
{

    char line_char_buf[1024]; // every char in a line will enter this buffer one line a time
    std::cout << "Attempting to parse input file and write to result file " << input_file << std::endl;

    // to parse data from lib file open it through fstream
    std::ifstream infile;
    infile.open(input_file);
    std::string templine; // stores entire lines as we progress through the document likely can remove later.
    char *temp_c_str = new char[1023];

    // state if we failed to open library
    if (!infile.is_open())
    {
        std::cout << "Couldn't open the library file " << input_file << std::endl;
        return;
    }

    // state if we failed to create/open our results file
    std::ofstream outfile(result_file);
    if (!outfile.is_open())
    {
        std::cout << "Couldn't open the result file " << result_file << std::endl;
        return;
    }

    if (infile.good() && outfile.good())
    {
        std::cout << "Able to parse input file" << input_file << " and write to result file " << result_file << std::endl;
    }

    while (infile.good() && outfile.good())
    {
        outfile.seekp(0);
        outfile << "<" << my_gate->GetGateCount() << " cell types>" << std::endl;
        outfile.seekp(0, std::ios::end);

        infile.getline(line_char_buf, 1023);    // read one line
        std::string line_string(line_char_buf); // convert to C++ string

        if (line_string.empty())
        { // is empty line?
            continue;
        }

        else
        {
            std::istringstream ss(line_string);
            std::string firstWord;
            ss >> firstWord;

            if (firstWord.find("cell") != std::string::npos)
            { // found the word cell

                // std::cout << "Found cell "<< std::endl;

                std::string cellName;
                ss >> cellName;

                if (firstWord.find("delay") != std::string::npos)
                {

                    // std::cout << "1"<< std::endl;
                    std::getline(infile, line_string);

                    // while (line_string.find("}") == std::string::npos && line_string.find("slew") == std::string::npos)
                    while (line_string.find("}") == std::string::npos)
                    {

                        std::strcpy(line_char_buf, line_string.c_str());
                        char *data_c_str = strtok(line_char_buf, " () , '\' \"\" ");

                        if (!line_string.empty())
                        {

                            if (line_string.find("values") != std::string::npos)
                            { // we are at the first row of delay values

                                // input 7x7 array into mygate delay
                                data_c_str = strtok(NULL, " () , '\' \"\"  \t");
                                for (int i = 0; i < 7; i++)
                                {

                                    while (data_c_str != NULL)
                                    {

                                        for (int j = 0; j < 7; j++)
                                        {
                                            // data_c_str = strtok(NULL, " () , '\' \"\"  ");
                                            std::string data_str = data_c_str;

                                            my_gate->SetDelay((my_gate->GetGateCount() - 1), i, j, std::stod(data_str));
                                            data_c_str = strtok(NULL, " () , '\' \"\"  \t");
                                        }
                                        data_c_str = NULL; // after we get all 7 elements of data use this to break while loop
                                    }

                                    std::getline(infile, line_string);
                                    std::strcpy(line_char_buf, line_string.c_str());
                                    temp_c_str = strtok(line_char_buf, " () , '\' \"\"  \t");

                                    if (temp_c_str == nullptr)
                                    {
                                        templine = "";
                                    }
                                    else
                                    {
                                        templine = temp_c_str;
                                        std::strcpy(line_char_buf, line_string.c_str());
                                        temp_c_str = strtok(line_char_buf, " () , '\' \"\" \t ");
                                        data_c_str = temp_c_str;
                                    }

                                    while (templine.empty())
                                    {
                                        std::getline(infile, line_string);
                                        std::strcpy(line_char_buf, line_string.c_str());
                                        temp_c_str = strtok(line_char_buf, " () , '\' \"\"  \t");

                                        if (temp_c_str == nullptr)
                                        {
                                            templine = "";
                                        }
                                        else
                                        {
                                            templine = line_string;
                                            data_c_str = temp_c_str;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                std::getline(infile, line_string); // havent reached values go to next line
                            }
                        }
                        std::getline(infile, line_string); // if we have reached cell_delay but there is an
                                                           // empty line before we get to data we want
                    }

                    // output 7x7 double array to result file
                    for (int i = 0; i < 7; i++)
                    {
                        for (int j = 0; j < 7; j++)
                        {
                            outfile << my_gate->GetDelay(my_gate->GetGateCount() - 1, i, j);
                            if (j < 6)
                            { // used to seperate data chuncks
                                outfile << ", ";
                            }
                        }
                        outfile << std::endl;
                    }
                    continue;
                }

                else
                {

                    char *temp_c_name = new char[line_string.size() + 1]; // +1 for the null terminator
                    std::strcpy(temp_c_name, line_string.c_str());
                    char *data_c_name = temp_c_name;
                    data_c_name = strtok(temp_c_name, " () ,  {");
                    while (temp_c_name != NULL)
                    {
                        data_c_name = strtok(NULL, "  ( ) ,  {");
                        temp_c_name = NULL;
                    }

                    std::string data_name = data_c_name;
                    outfile << "\n<" << data_name << ">" << std::endl;
                    my_gate->SetName(data_name, my_gate->GetGateCount());

                    my_gate->SetGateCount(my_gate->GetGateCount() + 1);
                    outfile.seekp(0);
                    outfile << "<" << my_gate->GetGateCount() << " cell types>" << std::endl;
                    outfile.seekp(0, std::ios::end);
                }
            }
        }
    }

    // outfile << "\n<EOF_test>" ;

    // make sure we close fstreams
    std::cout << my_gate->GetGateCount() << " gates were detected and corresponding delay data was saved and output in result file " << result_file << std::endl;
    infile.close();
    outfile.close();
}

int main(int argc, char **argv)
{

    // this bit of code is from the example code in FileIO.cpp on canvas
    if (argc < 2)
    {
        std::cout << "I need one parameter, which is the file name." << std::endl;
        return -1;
    }

    std::string input_file;
    std::string result_file = "wodzi003.txt";
    input_file = argv[1]; // uncomment after done debuging

    // input_file = "sample_NLDM.lib"; // this is only used for debug remove this later

    Gate gate;
    parse_file(input_file, result_file, &gate);
    return 0;
}