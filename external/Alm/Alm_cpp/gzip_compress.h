#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

bool compress_gzip(const std::string& source, const std::string& destination);
bool writeToFile(std::vector<double> x, std::vector<double> y, std::vector< std::vector<double> > z, std::string filename);
bool eraseFile(const std::string& filename);