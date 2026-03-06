// Fit KT amplitudes for π1(1600) decay with only P-wave to data from [1]
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------
// REFERENCES: 
// [1] - https://inspirehep.net/literature/1898933
// ------------------------------------------------------------------------------

#include <algorithm>
#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "colors.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"
#include "TRandom.h"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

inline std::string bootstrap_dir()
{
    // Find the correct data file using the top level repo directory
    std::string top_dir;
    char const * env = std::getenv("BOOTSTRAP");
    if ( env == NULL || std::string(env) == "" )
    {
        return iterateKT::error("bootstrap_dir(): Cannot find environment variable BOOTSTRAP!", "");
    }
    return std::string(env);  
};

int main(int argc, char** argv)
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    if (argc == 1) fatal("main()", "Need a FIT_ID argument!");
    int fit_id = atoi(argv[1]);

    std::ofstream out;
    out.open("log.txt");
    int  precision = 12, spacing = precision + 10;
    
    // Preamble info
    out << std::left << "# I ran and successfully exited and wrote to file" << std::endl;
    out << std::left << "# I was able to find the ITERATEKT directory (" << main_dir() << ")" <<std::endl;
    out << std::left << "# I was able to find the BOOTSTRAP directory (" << bootstrap_dir() << ")" <<std::endl;
    out.close();
};