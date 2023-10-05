#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/program_options.hpp>

#include "string_handler.h"
#include "readparams_job.h"

namespace po = boost::program_options;

void modifyParameter(std::unordered_map<std::string, std::string>& parameters,
                    const std::string& fileout,
                    const std::string& keyToModify,
                    const std::string& newValue) {

    // Check if the provided key exists in the parameters
    if (parameters.find(keyToModify) == parameters.end()) {
        std::cerr << "Key: " << keyToModify << " does not exist in the parameter file" << std::endl;
        exit(EXIT_FAILURE);
    }
    // Update the value associated with the key
    parameters[keyToModify] = newValue;
    // Write the modified parameters to a new file
    std::ofstream outFile(fileout);
    if (!outFile.is_open()) {
        std::cerr << "Failed to create new file for writing: " << fileout << std::endl;
        exit(EXIT_FAILURE);
    }
    for (const auto& parameter : parameters) {
        outFile << parameter.first << "=" << parameter.second << std::endl;
    }
    outFile.close();
}

int main(int argc, char* argv[]){
    double keyvalmin=0., keyvalmax=1., keystep=0.1;
    std::string key = "q_star";
    std::string file_in = "../tests/template.cfg";
    std::string cfg_file_tmp="../tests/tmp/template_mod.cfg";
    std::string file_out_root="out_";
    std::string path_solver="";
    int i;
    double keyval;
    std::string file_out;
    std::string keyval_str;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,H", "produce help message")
        ("param,P", po::value<std::string>(&key)->default_value("q_star"), "parameter name to be scanned")
        ("pmin", po::value<double>(&keyvalmin)->default_value(0.0), "minimum value for the scanned parameter")
        ("pmax", po::value<double>(&keyvalmax)->default_value(1.0), "maximum value for the scanned parameter")
        ("step,S", po::value<double>(&keystep)->default_value(0.1), "step size for the scanned parameter")
        ("input-template-file,I", po::value<std::string>(&file_in)->default_value("../config/scanner.cfg"), "input template file used for fixed parameters")
        ("tmp-cfg-file", po::value<std::string>(&cfg_file_tmp)->default_value("../tests/scanner/tmp/template_mod.cfg"), "file of the temporary config file that is created iteratively")
        ("file-results-root", po::value<std::string>(&file_out_root)->default_value("../tests/scanner/out/out_"), "output file root name")
        ("path-solver", po::value<std::string>(&path_solver)->default_value(""), "path to the ARMMsolver program");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << "This program allows you to run iteratively the ARMMSolver in order to scan a range of a given parameter" << std::endl;
        std::cout << "By default, it scans q_star and allows you to get the mixed modes frequencies as well as the zeta function " << std::endl;
        std::cout << "For further details on options, refers to the options explanations below " << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    // Read the template file and obtain the parameters
    std::unordered_map<std::string, std::string> parameters = readParameterFile(file_in);
    
    // The main loop
    keyval=keyvalmin;
    i=0;
    while(keyval <= keyvalmax){
        file_out = file_out_root + int_to_str(i) + ".res";
        std::cout << "[" << i << "] " << file_out << std::endl;
        // Make modification and write them in the file
        keyval_str=dbl_to_str(keyval);
        modifyParameter(parameters, cfg_file_tmp, key, keyval_str);
        // Make a call to the ARMMSolver
        std::string command = path_solver + "./ARMMSolver -I " + cfg_file_tmp + " -O " + file_out + " --verbose=false";
        system(command.c_str());
        keyval=keyval + keystep;
        i=i+1;
    }
}