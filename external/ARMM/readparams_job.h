/**
 * @file readparams_job.h
 * @brief Functions for reading parameters for make_star
 *
 * This file contains a function that can read a parameter file for make_star
 */
#pragma once 
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

/**
 * @brief Reads a parameter file and returns a map of key-value pairs.
 *
 * This function reads a parameter file and returns a map of key-value pairs. The parameter file should be in the format "key = value".
 *
 * @param filename The name of the parameter file.
 * @return std::unordered_map<std::string, std::string> A map of key-value pairs.
 */
std::unordered_map<std::string, std::string> readParameterFile(const std::string& filename);
