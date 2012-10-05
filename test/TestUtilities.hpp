
#ifndef TEST_UTILITIES_HPP
#define TEST_UTILITIES_HPP

#include "CubitBox.hpp"
#include <string>

// function to get the path to a data file in the data directory
std::string data_file(char* filename);

// compare if 2 CubitBoxes are identical
bool cubit_box_identical(const CubitBox& box1, const CubitBox& box2, double tol,
    bool print_data = false);



#endif

