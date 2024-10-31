#ifndef GATEDATABASE_HPP
#define GATEDATABASE_HPP

#include <string>
#include <vector>
#include <map>

#define GATE_LUT_DIM 7

struct GateInfo {
    double cell_delay[GATE_LUT_DIM][GATE_LUT_DIM];
    double output_slew[GATE_LUT_DIM][GATE_LUT_DIM];
    double capacitance;
    double delay_index_1 [GATE_LUT_DIM];
    double delay_index_2 [GATE_LUT_DIM];
    // delay_index_1 should always = slew_index_1 and the same for index 2 will alway just use delay_index_x for simplicity
    double slew_index_1 [GATE_LUT_DIM];
    double slew_index_2 [GATE_LUT_DIM];
};

class GateDatabase {
    private:
        // Stores the pointer to the GateInfo. Indexed using the name of gate
        // If you want to avoid pointers entirely, you may have to make GateInfo
        // with vectors instead of simple arrays
        // std::map<std::string, GateInfo*> gate_info_lut_;
    public:
        std::map<std::string, GateInfo*> gate_info_lut_;
        GateDatabase(const std::string& file_name);
        ~GateDatabase();

        void insert(const std::string& gateName, GateInfo* gateInfo);
        const GateInfo* get_gate_info(const std::string& gateName);

        void test();
};

#endif //GATEDATABASE_HPP
