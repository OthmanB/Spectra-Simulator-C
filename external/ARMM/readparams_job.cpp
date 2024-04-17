/**
 * @file readparams_job.cpp
 * @brief Functions for reading parameters for make_star
 *
 * This file contains a function that can read a parameter file for make_star
 */
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "string_handler.h"

/**
 * @brief Reads a parameter file and returns a map of key-value pairs.
 *
 * This function reads a parameter file and returns a map of key-value pairs. The parameter file should be in the format "key = value".
 *
 * @param filename The name of the parameter file.
 * @return std::unordered_map<std::string, std::string> A map of key-value pairs.
 */
std::unordered_map<std::string, std::string> readParameterFile(const std::string& filename) {
    std::unordered_map<std::string, std::string> parameters;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
		std::cerr << "Check that the file exist and that the path is correct" << std::endl;
   		std::cerr << "Cannot proceed" << std::endl;
   		std::cerr << "The program will exit now" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string line;
    while (std::getline(file, line)) {
        line =strtrim(line);
        if (line.empty() || line[0] == '#') {
            continue; // Ignore empty lines and comment lines
        }
        std::istringstream iss(line);
        std::string key;
        std::getline(iss, key, '=');
        key = key.substr(0, key.find_last_not_of(" \t") + 1);
        std::string value;
        std::getline(iss, value);
        value = value.substr(value.find_first_not_of(" \t"));
        parameters[key] = value;
        //std::cout << "key : START|" << key << "|END    value : START|" << value << "|END" << std::endl;
    }
    file.close();
    return parameters;
}
