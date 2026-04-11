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
    inline data_set  parse_JSON(uint m3pibin, std::string input)
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
        std::string bin = "m3pi_bin_number_" + to_string(m3pibin);
        auto m3pi_upper = data["bins"][bin]["bin_ranges"]["m3pi_upper_limit"].template get<double>();
        auto m3pi_lower = data["bins"][bin]["bin_ranges"]["m3pi_lower_limit"].template get<double>();
        double m3pi = (m3pi_upper + m3pi_lower)/2;
        
        // Calculate central t in bin
        auto t_upper = data["bins"][bin]["bin_ranges"]["t_upper_limit"].template get<double>();
        auto t_lower = data["bins"][bin]["bin_ranges"]["t_lower_limit"].template get<double>();
        double t = -(t_upper + t_lower)/2;
    
        std::string id = "m3π = " + to_string(m3pi,3) + ", t' = " + to_string(-t,3);

        auto bins      = data["bins"][bin]["bin_centers"];
        auto widths    = data["bins"][bin]["bin_widths"];
        auto abs_M     = data["bins"][bin]["abs_M"];
        auto std_abs_M = data["bins"][bin]["std_abs_M"];
        int N          = bins.size();

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(m3pi, M_PION);
        std::vector<double> sig1, sig2, absM, errM, bin_area, incorrect_bin_area;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double s1 = bins[i];
                double s2 = bins[j];
                s1 *= s1; s2 *= s2; // mass squared
    
                if (!kin->in_decay_region(s1, s2)) continue;
                if (are_equal(s1, s2))             continue;

                double z = abs_M[i][j];
                if (is_zero(z))                    continue;

                // Add data
                sig1.push_back(s1); sig2.push_back(s2);
                absM.push_back(     abs_M[i][j] ); 
                errM.push_back( std_abs_M[i][j] );

                // Because the plot is symmetric, add two copies with s1<->s2
                sig1.push_back(s2); sig2.push_back(s1);
                absM.push_back(     abs_M[i][j] ); 
                errM.push_back( std_abs_M[i][j] );
                double ds1 = widths[i], ds2 = widths[j];
                bin_area.push_back(ds1*ds2*4*sqrt(s1*s2));
                incorrect_bin_area.push_back(ds1*ds2);
            };
        };
        int N_actual = sig1.size();
        
        // Normalization of the data files is wrong here we fix it
        double wrong_norm = 0, corrected_norm = 0;
        for (int i = 0; i < absM.size(); i++)
        {
            wrong_norm += std::norm(absM[i])*incorrect_bin_area[i];
            corrected_norm += std::norm(absM[i])*bin_area[i];
        };
        double norm = sqrt(wrong_norm / corrected_norm);
        print(norm);
        // ---------------------------------------------------------------------------
        //  Organize everything
        out._N    = N_actual;         
        out._id   = id;               
        out._type = kDalitz;     
        out._extras["Nbins"] = N; 
        out._extras["m3pi"] = m3pi; 
        out._extras["t"]    = t;    
        out._extras["normalization"] = norm;
        out._x = sig1;  
        out._y = sig2;             
        out._z = norm*absM; out._dz = norm*errM;               

        return out;
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set parse_JSON(uint m3pi_bin, uint t_bin)
    {
        std::string st = to_string(t_bin);
        std::string filename = "tBin_"+st+".json";
        auto out = parse_JSON(m3pi_bin, filename);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };

    // ------------------------------------------------------------------------------------------------------------
    // These methods produce pseudodata from the original data file to run a bootstrap

    // Take in a data file from input and resample the data assuming gaussian errors
    inline data_set generate_pseudodata(uint m3pi_bin, std::string input, TRandom * rand)
    {
        // Final outputs
        data_set data = parse_JSON(m3pi_bin, input);

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = new_kinematics(data._extras["m3pi"], M_PION);
        std::vector<double> resampled_absM;

        // Random number generator with random seed

        for (int i = 0; i < data._z.size(); i++)
        {
            double mean      = data._z[i];
            double std_dev   = data._dz[i];

            double resampled = rand->Gaus(mean, std_dev);
            resampled_absM.push_back(resampled); 
        };

        data._z = resampled_absM;
        return data; 
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set generate_pseudodata(uint m3pi_bin, uint t_bin, TRandom * rand)
    {
        std::string st = to_string(t_bin);
        std::string filename = "tBin_"+st+".json";
        auto out = generate_pseudodata(m3pi_bin, filename, rand);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };

        // Take in a data file from input and resample the data assuming gaussian errors
    inline data_set generate_pseudodata(uint m3pi_bin, std::string input, TRandom * rand, amplitude amp)
    {
        // Final outputs
        data_set data = parse_JSON(m3pi_bin, input);

        // ---------------------------------------------------------------------------
        // Need to filter out any data outside of the physical kinematic region
    
        kinematics kin = amp->get_kinematics();
        std::vector<double> resampled_absM;

        for (int i = 0; i < data._z.size(); i++)
        {        
            // Here instead of saving abs_M, we resample it
            double mean      = data._z[i];
            double std_dev   = data._dz[i];
            
            double model     = abs(amp->evaluate_in_dalitz(data._x[i], data._y[i]));
            double chi2      = std::norm( (mean-model)/std_dev );
            double resampled = (chi2 > 1) ? rand->Gaus(mean, sqrt(chi2)*std_dev) : rand->Gaus(mean, std_dev);

            resampled_absM.push_back(resampled); 
        };
        data._z = resampled_absM;
        return data; 
    };

    // Do the above but input bin numbers IDs which are subsequently saved in the data_set
    inline data_set generate_pseudodata(uint m3pi_bin, uint t_bin, TRandom * rand, amplitude amp)
    {
        std::string st = to_string(t_bin);
        std::string filename = "tBin_"+st+".json";

        // Set the amplitude to the correct bin
        amp->set_option(option::set_tbin,         t_bin);
        amp->set_option(option::set_mbin_COMPASS, m3pi_bin);

        auto out = generate_pseudodata(m3pi_bin, filename, rand, amp);
        out._extras["t_bin"]    = t_bin;
        out._extras["m3pi_bin"] = m3pi_bin;
        return out;
    };
}; };

#endif 
