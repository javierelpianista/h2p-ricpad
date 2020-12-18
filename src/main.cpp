#include <iostream>
#include <array>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/multiprecision/mpc.hpp>
#include <ginac/ginac.h>

#include <h2p.hpp>
#include <req.hpp>
#include <fixed.hpp>

using std::cout;
using std::endl;

namespace po = boost::program_options;
namespace mp = boost::multiprecision;
namespace gi = GiNaC;

using mp::mpc_complex;
using mp::mpfr_float;

po::options_description 
    mandatory("Required options"), 
    optional("Non-mandatory options"), 
    hidden, opts;
po::positional_options_description p_desc;
po::variables_map vm;

const auto help_message = "Usage: h2p-ricpad MODE OPTIONS\n\n"
        "MODE can be either 'minimum' for computing U, A, and "
        " R simultaneously for the equilibrium internuclear distance, "
        " or 'fixed', for computing U and A for a given R.\n\n";

const std::array<const std::string, 2> required_options { "Dmin", "variable" }; 

void print_help_message(std::string const &a = "") {
    cout << help_message;
    cout << mandatory << endl;
    cout << optional << endl;
}

int main(int argc, char* argv[]) {
    mandatory.add_options() 
        ("Dmin", po::value<int>(), "Minimum D value")
        ("d", po::value<int>()->default_value(0), "Value of d")
        ("m", po::value<int>()->default_value(0), 
         "Value of the quantum number associated with the angular part of the"
         " spheroidal equations")
        ("s", po::value<int>()->default_value(0),
         "Value stating if the second spheroidal equation's eigenfuncions are"
         " even (s=0) or odd (s=1).")
        ("U0", po::value<std::string>(), 
         "Initial value of the electronic+nuclear energy")
        ("A0", po::value<std::string>(),
         "Initial value of the coupling constant")
        ("R0", po::value<std::string>(),
         "Initial value of the internuclear distance")
        ;
    optional.add_options()
        ("help", po::value<std::string>()
         ->zero_tokens()
         ->implicit_value("")
         ->notifier(&print_help_message)
         , "Print this message.")
        ("Dmax", po::value<int>()->default_value(-1), "Maximum D value")
        ("ndigits", po::value<int>()->default_value(1000),
         "Number of digits for the numerical calculations")
        ("tol", po::value<std::string>()->default_value("1E-100"),
         "Tolerance for the Newton-Raphson method")
        ("h", po::value<std::string>()->default_value("1E-200"),
         "Step size for the Newton-Raphson method")
        ("h2", po::value<std::string>()->default_value("1E-400"),
         "Step size for the Newton-Raphson for the computation of numerical derivations outside the Newton-Raphson iterations")
        ("use-E", po::bool_switch()->default_value(false), 
         "Set this option if the provided value of U0 is the electronic "
         "energy, instead of the electronic + nuclear one")
        ;

    hidden.add_options()
        ("mode", "" ) 
        ;

    p_desc.add("mode", -1);

    opts.add(mandatory).add(optional).add(hidden);

    po::store(
        po::command_line_parser(argc, argv)
        .options(opts).positional(p_desc).run(), vm
        );

    po::notify(vm);

    // First check if the user asked for help or if they didn't set the 
    // mode correctly.
    const std::vector<std::string> available_modes = {
        "minimum", "fixed"
    };

    if ( vm.count("mode") ) {
        std::string mode = vm["mode"].as<std::string>();
        if ( std::none_of(
                available_modes.begin(),
                available_modes.end(),
                [&mode](const std::string &a) -> bool {return a == mode;}
                )
           ) {
            cout << "Mode " << mode << " not available." << endl;
            print_help_message("");
        }
    }

    // Mandatory values
    const std::vector<std::string> mandatory_options = {
        "Dmin", "U0", "A0", "R0"
    };

    for ( auto what: mandatory_options ) {
        if ( vm.count(what) == 0 ) {
            cout << "Mandatory option " << what << " not provided." << endl;
            print_help_message();
            return 1;
        }
    }

    // Dmin
    int Dmin = vm["Dmin"].as<int>();

    // Dmax
    int Dmax;
    if ( vm.count("Dmax") ) {
        Dmax = vm["Dmax"].as<int>();
    } else { 
        Dmax = -1;
    }

    // Quantum numbers
    int m = vm["m"].as<int>();
    int s = vm["s"].as<int>();

    mpfr_float tol, h, h2;

    // Digits for calculations
    { 
        int ndigits = vm["ndigits"].as<int>();
        mpfr_float::default_precision(ndigits);
        mpc_complex::default_precision(ndigits);
        gi::Digits = ndigits;
        tol = mpfr_float(vm["tol"].as<std::string>());
        h   = mpfr_float(vm["h"].as<std::string>());
        h2  = mpfr_float(vm["h2"].as<std::string>());
    }

    // Initial values
    mpfr_float 
        U0(vm["U0"].as<std::string>()),
        A0(vm["A0"].as<std::string>()),
        R0(vm["R0"].as<std::string>());

    // Check if the U0 is the electronic+nuclear energy or just the electronic
    // one
    bool use_E = vm["use-E"].as<bool>();

    cout.precision(40);

    if ( vm["mode"].as<std::string>() == "minimum" ) {
        R_eq<mpfr_float>(Dmin, Dmax, U0, A0, R0, 0, false, tol, h, h2);
        return 0;
    } else if ( vm["mode"].as<std::string>() == "fixed" ) {
        UA<mpfr_float>(Dmin, Dmax, m, s, U0, A0, R0, 0, tol, h, h2, use_E);
        return 0;
    } else {
        return 1;
    }
}
