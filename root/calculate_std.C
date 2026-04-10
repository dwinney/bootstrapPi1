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
#include <limits>
#include "TTree.h"

#define N_BIN_MIN 11
#define N_BIN_MAX 49

void calculate_std()
{
    string branch_prefix  = "modNc";
    string file_prefix    = "nonminimal";

    cout << "Variable ? ";
    cin  >> branch_prefix;

    cout << "Minimal or nonminimal ? ";
    cin  >> file_prefix;
    cout << endl;
    cin.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' ); 

    // ----------------------------------------------------------------------------------
    // Read in the best fit pars to compare
    
    string in_file  = file_prefix+".root" ;
    string bff_file = "../data/"+file_prefix+".dat";
    string out_file = file_prefix+ "_"+branch_prefix+".dat";
    
    unique_ptr<TFile> rootFile( TFile::Open(in_file.c_str()) );
    auto tree = rootFile->Get<TTree>("bootstrap");

    // ----------------------------------------------------------------------------------
    // Read in the best fit pars to compare
    
   vector<double>  m3pis, bffs;
    
   std::ifstream infile(bff_file);
   if (!infile) {cout << "Cannot open file " + bff_file + "!"; exit(0);};
   std::string line;
    
   int nimported = 0; // Mark how many lines we've imported
   while (std::getline(infile, line))
   {   
       if (line.empty())        continue; // skips empty lines
       if (line.front() == '#') continue; // Skip comment lines 
       std::istringstream is(line);  
       
       if (nimported <= N_BIN_MAX-N_BIN_MIN)
       {
           double trash, m3pi, alpha, beta, redelta, imdelta;
           // Dont care about first two columns
           is >> trash >> m3pi;
           m3pis.push_back(m3pi);

           // we do about these though
           is >> alpha >> beta >> trash >> redelta >> imdelta;
           
           if (branch_prefix == "modNc")  bffs.push_back(alpha);
           if (branch_prefix == "modNcp") bffs.push_back(beta);
           
           complex<double> delta(redelta, imdelta);
           if (branch_prefix == "modNd") bffs.push_back(abs(delta));
           double argd = (arg(delta) > 0) ? arg(delta) - 2*M_PI : arg(delta);
           if (branch_prefix == "argNd") bffs.push_back(argd);
           nimported++;
           continue;
        };
    };

    // ----------------------------------------------------------------------------------
    double par;
    double min, max;
    vector<double> stds;
    
    for (int n = N_BIN_MIN; n <= N_BIN_MAX; n++)
    {
        
        tree->SetBranchStatus("*", false);
        std::string branch = branch_prefix+"_"+to_string(n);
        tree->SetBranchStatus(branch.c_str(), true);
        tree->SetBranchAddress(branch.c_str(), &par);
        
        tree->Draw(branch.c_str());
        gPad->Update();
        min = tree->GetMinimum(branch.c_str());
        max = tree->GetMaximum(branch.c_str());
        
        cout << "\n" << branch << ": \n";
        cout << left << setw(15) << "Best fit found: " << setw(15) << bffs[n-N_BIN_MIN] << endl;
        cout << left << setw(15) << "Cut ? ";
        
        bool cut = false;
        string input = "";
        std::getline( cin, input );
        if ( !input.empty() )
        {
            istringstream stream (input);
            stream >> cut;
        };

        if (cut)
        {
            cout << left << setw(15) << "Minumum ? ";
            cin >> min;
            cout << left << setw(15) << "Maximum ? ";
            cin >> max;
            cin.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' ); 
        };

        TH1D *h1 = new TH1D("h1", "after selections", 1000, min, max);
        for (int i = 0; i < tree->GetEntries(); ++i)
        {
            tree->GetEntry(i);
            if (cut)
            {
                if (min <= par && max >= par) h1->Fill(par);
            }
            else h1->Fill(par);
        };
        double std =  h1->GetStdDev();
        cout << left << setw(15) << "mean: " << setw(20) << h1->GetMean() << "\n";
        cout << left << setw(15) << "std: "  << setw(20) << std << "\n";
        cout << endl;
        stds.push_back(std);

        delete h1;
    };

    std::ofstream out;
    out.open(out_file);
    int precision = 9;
    int spacing = precision + 10;
    out << left << std::setw(spacing) << "# m3pi" << std::setw(spacing) << "BFF" << std::setw(spacing) << "std" << endl;
    for (uint i = 0; i < m3pis.size(); i++)
    {
        out << std::left << std::setprecision(precision);
        out << std::setw(spacing) << m3pis[i] << std::setw(spacing) << bffs[i] << std::setw(spacing) << stds[i] << endl;
    };
    out.close();

    gApplication->Terminate(0);
};
