#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio> 
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>

bool compress_gzip(const std::string& source, const std::string& destination) {
    // Set up input stream
    std::ifstream file_in(source, std::ios_base::in | std::ios_base::binary);
    if (!file_in) {
        std::cerr << "Error: Could not open input file " << source << std::endl;
        return 1;
    }

    // Set up output stream
    std::ofstream file_out(destination, std::ios_base::out | std::ios_base::binary);
    if (!file_out) {
        std::cerr << "Error: Could not open output file " << destination << std::endl;
        return 1;
    }

    // Set up gzip compression filter
    boost::iostreams::filtering_ostream output;
    output.push(boost::iostreams::gzip_compressor());
    output.push(file_out);

    // Copy input to output using gzip compression
    try {
        boost::iostreams::copy(file_in, output);
    } catch (const boost::iostreams::gzip_error& e) {
        std::cerr << "Error: Gzip compression error: " << e.what() << std::endl;
        return 1;
    }

    // Close streams
    file_in.close();
    output.reset();
    file_out.close();
    return 0;
}


bool writeToFile(std::vector<double> x, std::vector<double> y, std::vector< std::vector<double> > z, std::string filename) {
    std::ofstream file;
    file.open(filename);
    if (file.is_open()) {
        file << "x=";
        for (double i : x) {
            file << i << " ";
        }
        file << "\n";
        
        file << "y=";
        for (double i : y) {
            file << i << " ";
        }
        file << "\n";
        file << "z=\n";
        for (std::vector<double> vec : z) {
            for (double i : vec) {
                file << i << " ";
            }
            file << "\n";
        }
        file.close();
    }
    else {
        std::cout << "Unable to open file";
        return 1;
    }
    return 0;
}

bool eraseFile(const std::string& filename) {
    if (std::remove(filename.c_str()) != 0) {
        // Failed to delete file
        std::cerr << "Error deleting file " << filename << std::endl;
        return 1;
    } else {
        //std::cout << "File " << filename << " deleted successfully." << std::endl;
        return 0;
    }
}
