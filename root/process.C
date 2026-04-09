// ROOT script to parse fit_result.txt files into a big root tree for 
// easier histogram visualization
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Autónoma Nacional de México (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include <iostream> 
#include <cmath>
#include <string>
#include <sstream>
#include "TTree.h"

#define N_BIN_MIN 11
#define N_BIN_MAX 49

void process()
{
    bool minimal  = false;
    std::cout << "Minimal ? "; 
    std::cin  >> minimal; 
    // Operating parameters
    std::string path     = std::string(std::getenv("BOOTSTRAP"));
    int nBS              = 1E4;  // Total number of bootstraps

    // Output name
    std::string out_file = (minimal) ? "minimal.root" : "nonminimal.root";

    // Set up ROOT file
    TH1::AddDirectory( kFALSE );
    gSystem->Load( "libTree" );

    std::unique_ptr<TFile> file( TFile::Open( out_file.c_str(), "recreate" ));
    std::unique_ptr<TTree> tree = std::make_unique<TTree>( "bootstrap", "bootstrap");

    // Temporary save locations for parameters
    std::array<double,N_BIN_MAX-N_BIN_MIN+1> Nc, Ncp, mNd, aNd;
    double c, cp, d, fcn;

    tree->Branch("fcn", &fcn, "fcn/D");
    for (int index = 0; index <= N_BIN_MAX-N_BIN_MIN; index++)
    {
        std::string si = std::to_string(index+N_BIN_MIN);
        tree->Branch( ("modNc_"  +si).c_str(),  &Nc[index], ("|Nc_"+si+"|/D").c_str());
        if (!minimal)
        {
            tree->Branch( ("reNcp_" +si).c_str(), &Ncp[index], ("reNcp_"+si+"/D").c_str());
        };
        tree->Branch( ("modNd_"+si).c_str(), &mNd[index], ("modNd_"+si+"/D").c_str());
        tree->Branch( ("argNd_"+si).c_str(), &aNd[index], ("argNd_"+si+"/D").c_str());
    };
    tree->Branch("c", &c, "c/D");
    if (!minimal) tree->Branch("cp", &cp, "cp/D");
    tree->Branch("d", &d, "d/D");
    
    // Import data
    for (int fit = 0; fit <= nBS; fit++)
    {
     
	if (fit % 1000 == 0) std::cout << "Processed " << fit << " bootstraps!" << std::endl;
	std::string in_file  = path+"/fit_"+std::to_string(fit)+"/fit_results.dat";
        std::ifstream infile(in_file);
        if (!infile) continue;

        std::string line;
        int n = 0;
        while (std::getline(infile, line))
        {   
            if (line.empty()) continue; // skips empty lines
            
            std::istringstream is(line);  
            std::string trash;
            if (line.front() == '#')
            {
                is >> trash >> trash; 
                if (trash == "average") is >> trash >> trash >> trash >> fcn;
                continue;
            };

            if (n <= N_BIN_MAX-N_BIN_MIN)
            {
                double rNd, iNd;
                is >> trash >> trash >> Nc[n] >> Ncp[n] >> trash >> rNd >> iNd;
                std::complex<double> Nd(rNd, iNd);
                mNd[n] = abs(Nd); 
                double argNd = arg(Nd);
                aNd[n] = (argNd > 0) ? argNd - 2*M_PI : argNd;
                n++;
                continue;
            };
            is >> c >> cp >> d;  
            break;
        };
        tree->Fill();
    };
    tree->Write();
    tree->AutoSave();
};
