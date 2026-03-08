// Some utilities to find directories
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

#ifndef PI1_UTILITIES
#define PI1_UTILITIES

#include <string>

namespace iterateKT
{
    std::string data_dir()
    {
        // Find the correct data file using the top level repo directory
        std::string top_dir;
        char const * env = std::getenv("BOOTSTRAP");
        if ( env == NULL || std::string(env) == "" ) exit(1);
        return std::string(env) + "/data/";  
    };
};

#endif
