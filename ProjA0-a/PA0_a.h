#ifndef _PA0_A_H
#define _PA0_A_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

class Gate
{

private:
    int gate_count;
    std::string gate_name[100];
    double delay[100][7][7]; // we can have a maximum of 100 gates each containing 7 x 7 arrays

public:
    // constructor
    Gate()
    {
        gate_count = 0;
        for (int i = 0; i < 100; i++)
        {
            for (int j = 0; j < 7; j++)
            {
                for (int k = 0; k < 7; k++)
                {
                    this->delay[i][j][k] = 0; // init all gate delay timmings to 0
                }
            }
            gate_name[i] = "undeclared";
        }
    }

    // basics setter and getter methods to access private data
    void SetName(std::string Gatename, int count)
    {
        this->gate_name[count] = Gatename;
    }
    void SetGateCount(int new_count)
    {
        this->gate_count = new_count;
    }

    void SetDelay(int i, int j, int k, double gate_delay)
    {
        this->delay[i][j][k] = gate_delay;
    }

    int GetGateCount(void)
    {
        return gate_count;
    }
    double GetDelay(int i, int j, int k)
    {
        return delay[i][j][k];
    }
    std::string GetName(int gate_num)
    {
        return gate_name[gate_num];
    }
};

#endif