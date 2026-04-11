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

#include "pi1_utilities.hpp"
#include "pi1_isobar.hpp"
#include "pi1_amplitude.hpp"
#include "pi1_fitter.hpp"
#include "pi1_isobar.hpp"

int main()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // If to use three basis functions or only two
    bool minimal = false;

    // Which range of m3pi bins to consider
    int min = 11, max = 49; 

    // Fitter's stopping tolerance
    double tolerance = 0.0001;

    // Path to precalculated isoabrs
    std::string iso_path    = data_dir()+"basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CCD";
    // Are we taking initial values from file? if so which?
    std::string in_pars_file = (minimal) ? data_dir()+"minimal.dat" : data_dir()+"nonminimal.dat";
    // Where do we export the fit parameter values
    std::string out_pars_file  = current_dir()+"fit_results.dat";
    // Put a file description at the beginning
    std::string description = "contact + contact + deck, no form factor";
    
    // -----------------------------------------------------------------------
    // Import initial values
    
    auto x = COMPASS::import_parameters({min, max}, in_pars_file, minimal);
    auto m3pis            = std::get<0>(x);
    auto initial_fit_pars = std::get<1>(x);

    auto y   = COMPASS::import_parameters({min, max}, in_pars_file, false);
    auto BFF = std::get<1>(y);

    // Import our data sets  
    std::vector<std::string> labels;        // parameter labels
    for (int i = min; i <= max; i++)
    {
        labels.push_back("alpha_"+to_string(i));
        labels.push_back("beta_" +to_string(i));
        labels.push_back("delta_"+to_string(i));
    };
    labels.push_back("b_alpha");
    labels.push_back("b_beta");
    labels.push_back("b_delta");
    
    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    // Set up our amplitude (uniterated)
    auto args     = std::make_tuple(m3pis, COMPASS::t_bins);
    amplitude amp = new_amplitude<pi1_binned>(nullptr, args);
    amp->set_name("π₁ → 3π");
    
    // and import the pre-calculated isobars
    amp->import_solution(iso_path+file_prefix);
    
    // precalculate isobars in the decay region
    amp->precompute_dalitz(300);
    
    // Load best-fit parameters into amplitude
    COMPASS::fit_2D::process_parameters(BFF, amp);
    
    // -----------------------------------------------------------------------
    // Data set up

    
    TRandom * rand = new TRandom(0);
    rand->WriteRandom("rand_generator");

    // Import our data sets  
    std::vector<data_set> data; 
    for (int i = min; i <= max; i++)
    {
       for (int j = 0; j < 4; j++) data.emplace_back(COMPASS::generate_pseudodata(i, j, rand, amp));
    };

    // -----------------------------------------------------------------------
    // Set up fitter

    fitter<amplitude,COMPASS::fit_2D> fitter(amp, "Combined");
    fitter.set_tolerance(tolerance/2*1E3);
    fitter.set_print_level(2);
    fitter.set_strategy(0);

    // Add all bins
    fitter.add_data(data);

    // Add three t-slopes in addition to three subtraction coeffs
    fitter.add_extra_parameters(3);

    fitter.set_parameter_labels(labels);
    // Fix alphas to all be real (and positive)
    for (int i = min; i <= max; i++)
    {
        fitter.fix_argument("alpha_"+to_string(i), 0.); 
        fitter.fix_argument("beta_"+to_string(i) , 0.);
    };
    // t-slopes as well
    fitter.make_real("b_alpha"); 
    fitter.make_real("b_beta"); 
    fitter.make_real("b_delta"); 

    if (minimal)
    {
        for (int i = min; i <= max; i++) fitter.fix_parameter("beta_"+to_string(i), 0.);
        fitter.fix_parameter("b_beta", 0.);
    };

    fitter.do_fit(initial_fit_pars);
    
    // -----------------------------------------------------------------------
    // Print fit results to out_file
    COMPASS::export_parameters({min, max}, fitter.pars(), 
                               "average χ² per Dalitz: "+iterateKT::to_string(fitter.fcn()), 
                               out_pars_file);

    return 1;
};
