// Methods to interface with COMPASS data sets in JSON format
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2025)
// Affiliation:  Universitat Bonn
//               Helmholtz Institute (HISKP)
// Email:        daniel.winney@gmail.com
// ------------------------------------------------------------------------------

#ifndef COMPASS_DATA_HPP
#define COMPASS_DATA_HPP

#include "nlohmann/json.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <tuple>
#include "constants.hpp"
#include "kinematics.hpp"
#include "utilities.hpp"
#include "data_set.hpp"
#include "TRandom.h"

#include "pi1_utilities.hpp"

using json = nlohmann::json;

namespace iterateKT { namespace COMPASS
{
    // Static identifiers for data_set types
    static const int kReal = 0, kImag = 1, kAbs = 2, kReal1D = 3, kImag1D = 4, kDalitz = 5;
    
    // Parse a JSON file importing everything in a data_set object
    // Columns correspond to: s, t, Abs(M), Err(M)
    inline data_set  parse_JSON(std::string input)
    {
        // Final outputs
        data_set out;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = data_dir() + "raw_files/" + input;
        std::ifstream raw_file(path_to_file);
        if (!raw_file) fatal("Could not open file: " + path_to_file);
        json data = json::parse(raw_file);
    
        // Calculate central m3pi in bin
        auto m3pi_upper = data["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

        auto bins      = data["bin_centers"];
        auto abs_M     = data["abs_M"];
        auto std_abs_M = data["std_abs_M"];
        int N          = bins.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(m3pi, M_PION);
        std::vector<double> sig1, sig2, absM, errM;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double s1 = bins[i];
                double s2 = bins[j];
                s1 *= s1; s2 *= s2; // mass squared
    
                if (!kin->in_decay_region(s1, s2)) continue;
                if (are_equal(s1, s2))             continue;
                
                // Add data
                sig1.push_back(s1); sig2.push_back(s2);
                absM.push_back(     abs_M[i][j] ); 
                errM.push_back( std_abs_M[i][j] );

                // Because the plot is symmetric, add two copies with s1<->s2
                sig1.push_back(s2); sig2.push_back(s1);
                absM.push_back(     abs_M[i][j] ); 
                errM.push_back( std_abs_M[i][j] );
            };
        };
        int N_actual = sig1.size();

        // ---------------------------------------------------------------------------
        //  Organize everything
        out._N    = N_actual;         
        out._id   = id;               
        out._type = kDalitz;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"] = m3pi; 
        out._extras["t"]    = t;    
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out;
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set parse_JSON(uint m3pi_bin, uint t_bin)
    {
        std::string sm3pi = to_string(m3pi_bin), st = to_string(t_bin);
        std::string filename = "tBin_"+st+"/dalitz_m3piBin_"+sm3pi+"_tBin_"+st+".json";
        auto out = parse_JSON(filename);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };

    // ------------------------------------------------------------------------------------------------------------
    // These methods produce pseudodata from the original data file to run a bootstrap

    // Take in a data file from input and resample the data assuming gaussian errors
    inline data_set generate_pseudodata(std::string input, TRandom * rand)
    {
        // Final outputs
        data_set out;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = data_dir() + "raw_files/" + input;
        std::ifstream raw_file(path_to_file);
        if (!raw_file) fatal("Could not open file: " + path_to_file);
        json data = json::parse(raw_file);
    
        // Calculate central m3pi in bin
        auto m3pi_upper = data["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

        auto bins      = data["bin_centers"];
        auto abs_M     = data["abs_M"];
        auto std_abs_M = data["std_abs_M"];
        int N          = bins.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(m3pi, M_PION);
        std::vector<double> sig1, sig2, absM, errM;

        // Random number generator with random seed

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double s1 = bins[i];
                double s2 = bins[j];
                s1 *= s1; s2 *= s2; // mass squared
    
                if (!kin->in_decay_region(s1, s2)) continue;
                if (are_equal(s1, s2))             continue;
                
                
                // Here instead of saving abs_M, we resample it
                double mean      = abs_M[i][j];
                double std_dev   = std_abs_M[i][j];
                if ( is_zero(std_dev) ) continue;

                double resampled = rand->Gaus(mean, std_dev);

                sig1.push_back(s1); sig2.push_back(s2);
                absM.push_back(resampled); 
                errM.push_back(std_dev);

                // Symmetric errors
                sig1.push_back(s2); sig2.push_back(s1);
                absM.push_back(resampled); 
                errM.push_back(std_dev);
            };
        };
        int N_actual = sig1.size();

        // ---------------------------------------------------------------------------
        //  Organize everything
        out._N    = N_actual;         
        out._id   = id;               
        out._type = kDalitz;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"]  = m3pi; 
        out._extras["t"]     = t;    
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out; 
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set generate_pseudodata(uint m3pi_bin, uint t_bin, TRandom * rand)
    {
        std::string sm3pi = to_string(m3pi_bin), st = to_string(t_bin);
        std::string filename = "tBin_"+st+"/dalitz_m3piBin_"+sm3pi+"_tBin_"+st+".json";
        auto out = generate_pseudodata(filename, rand);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };

        // Take in a data file from input and resample the data assuming gaussian errors
    inline data_set generate_pseudodata(std::string input, TRandom * rand, amplitude amp)
    {
        // Final outputs
        data_set out;

        // ---------------------------------------------------------------------------
        // Read in json and organize everything 

        std::string path_to_file = data_dir() + "raw_files/" + input;
        std::ifstream raw_file(path_to_file);
        if (!raw_file) fatal("Could not open file: " + path_to_file);
        json data = json::parse(raw_file);
    
        // Calculate central m3pi in bin
        auto m3pi_upper = data["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

        auto bins      = data["bin_centers"];
        auto abs_M     = data["abs_M"];
        auto std_abs_M = data["std_abs_M"];
        int N          = bins.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = amp->get_kinematics();
        std::vector<double> sig1, sig2, absM, errM;

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                double s1 = bins[i];
                double s2 = bins[j];
                
                if (are_equal(s1, s2))             continue;
                if (!kin->in_decay_region(s1, s2)) continue;
                
                s1 *= s1; s2 *= s2; // mass squared
                
                // Here instead of saving abs_M, we resample it
                double mean      = abs_M[i][j];
                double std_dev   = std_abs_M[i][j];
                if (iszero(std_dev)) continue;
                
                double model     = abs(amp->evaluate(s1, s2, kin->Sigma()-s1-s2));
                double chi2      = std::norm( (mean-model)/std_dev );
                double resampled = (chi2 > 1) ? rand->Gaus(mean, sqrt(chi2)*std_dev) : rand->Gaus(mean, std_dev);

                sig1.push_back(s1); sig2.push_back(s2);
                absM.push_back(resampled); 
                errM.push_back(std_dev);

                // Symmetric errors
                sig1.push_back(s2); sig2.push_back(s1);
                absM.push_back(resampled); 
                errM.push_back(std_dev);
            };
        };
        int N_actual = sig1.size();

        // ---------------------------------------------------------------------------
        //  Organize everything
        out._N    = N_actual;         
        out._id   = id;               
        out._type = kDalitz;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"]  = m3pi; 
        out._extras["t"]     = t;    
        out._x = sig1;  
        out._y = sig2;             
        out._z = absM; out._dz = errM;               

        return out; 
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set generate_pseudodata(uint m3pi_bin, uint t_bin, TRandom * rand, amplitude amp)
    {
        std::string sm3pi = to_string(m3pi_bin), st = to_string(t_bin);
        std::string filename = "tBin_"+st+"/dalitz_m3piBin_"+sm3pi+"_tBin_"+st+".json";

        // Set the amplitude to the correct bin
        amp->set_option(option::set_tbin,         t_bin);
        amp->set_option(option::set_mbin_COMPASS, m3pi_bin);

        auto out = generate_pseudodata(filename, rand, amp);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };
}; };

#endif 