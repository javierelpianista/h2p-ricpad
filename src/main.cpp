#include <iostream>
#include <array>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

using std::cout;
using std::endl;

namespace po = boost::program_options;

po::options_description desc("Allowed options");
po::variables_map vm;

const std::array<const std::string, 2> required_options { "Dmin", "variable" }; 

int main(int argc, char* argv[]) {
    desc.add_options() 
        ("help", "produce help message")
        ("Dmin", po::value<int>(), "Minimum D value")
        ("Dmax", po::value<int>()->default_value(-1), "Maximum D value")
        ("variable", po::value<std::vector<std::string>>()->composing(), 
            "Symbols representing the independent variables")
        ;

    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("Dmin")) {
        cout << "Dmin was set to " 
        << vm["Dmin"].as<int>() << ".\n";
    } else {
        cout << "Dmin was not set.\n";
    }

    if (vm.count("Dmax")) {
        cout << "Dmax was set to " 
        << vm["Dmax"].as<int>() << ".\n";
    } else {
        cout << "Dmax was not set.\n";
    }

    for ( auto var : vm["variable"].as< std::vector<std::string>>() ) {
        cout << var << " ";
    }
    cout << endl;

    return 0;
}
