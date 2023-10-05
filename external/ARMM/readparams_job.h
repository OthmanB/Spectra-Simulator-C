#pragma once 
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

std::unordered_map<std::string, std::string> readParameterFile(const std::string& filename);
