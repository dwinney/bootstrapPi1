// Take results from a multi-dimensional fit and plot the integrated widths
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Nacional Autonoma de Mexico (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include "kinematics.hpp"
#include "amplitude.hpp"
#include "utilities.hpp"
#include "constants.hpp"
#include "plotter.hpp"
#include "fitter.hpp"

#include "amplitudes/pi1.hpp"
#include "isobars/pi1.hpp"
#include "COMPASS_pi1/fitter.hpp"
#include "COMPASS_pi1/data.hpp"

void plot_width_m3pi()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which range of m3pi bins to consider
    int min = 13, max = 49; 

    // Path to precalculated isoabrs
    std::string iso_path    = main_dir()+"/analysis/COMPASS_pi1/basis_functions/";
    // and the prefix given to each file
    std::string file_prefix = "CD";
    // path to par file
    std::string path =  gSystem->Getenv("BOOTSTRAP_SRC");
    std::string par_file    = path+"/data/pars/best_fit.dat";

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::array<std::vector<data_set>,4> data;
    for (int j = 0; j < 4; j++)
    {
        for (int i = min; i <= max; i++) data[j].emplace_back(COMPASS::parse_JSON(i, j));
    };

    // -----------------------------------------------------------------------
    // Import fit values

    std::vector<complex> pars    = COMPASS::import_parameters({min,max}, par_file);
    std::vector<complex> pars_c, pars_d;
    for (int i = 0; i < (pars.size()-2)/2; i++)
    {
        // Contact only
        pars_c.push_back( pars[2*i]   );
        pars_c.push_back( 0. );
        // Deck only
        pars_d.push_back( 0. );
        pars_d.push_back( pars[2*i+1]  );
    };
    for (int n = 0; n < 2; n++)
    {
        int m = pars.size()-2;
        pars_c.push_back(pars[m+n]);
        pars_d.push_back(pars[m+n]);
    };

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::vector<double> m3pis = COMPASS::m3pi_bins(min, max);

    // Set up our amplitude (uniterated)
    amplitude amp  = new_amplitude<pi1_binned>(nullptr, std::make_tuple(m3pis, COMPASS::t_bins));
    amp->import_solution(iso_path+file_prefix);
    amp->precompute_dalitz(300);
    
    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution
    
    std::array<std::vector<double>,4> mws, ews, mws_c, mws_d;
    for (int j = 0; j < 4; j++)
    {
        for (auto bin : data[j])
        {
            double ew = 0, mw = 0, nw = 0;
            for (int i = 0; i < bin._z.size(); i++) ew += norm(bin._z[i])*bin._dx[i];
            ews[j].push_back(ew/1E4);
            
            COMPASS::fit_2D::process_parameters(pars, amp);
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            mws[j].push_back(amp->width()/1E4);
            COMPASS::fit_2D::process_parameters(pars_c, amp);
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            mws_c[j].push_back(amp->width()/1E4);
            COMPASS::fit_2D::process_parameters(pars_d, amp);
            amp->set_option(option::set_tbin,         bin._extras["t_bin"]);
            amp->set_option(option::set_mbin_COMPASS, bin._extras["m3pi_bin"]);
            mws_d[j].push_back(amp->width()/1E4);
        };
    };
        
    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    auto blue  = square(jpacColor::Blue, "Contact #times 0.1");
    blue._draw_opt = "PL";
    auto red   = square(jpacColor::Red, "Deck #times 0.1");
    red._draw_opt = "PL";
    auto green = square(jpacColor::Green, "Full");
    green._draw_opt = "PL";

    // Plot widths as a function of m3pi
    plot p1 = plotter.new_plot();
    p1.set_legend(0.2, 0.65);
    p1.set_ranges({m3pis[0], m3pis.back()}, {0, 8.7});
    p1.add_header("#minus #it{t} = 0.12 GeV^{2}");
    p1.add_data(m3pis, ews[0],    dot(jpacColor::DarkGrey, "Data"));
    p1.add_data(m3pis, mws[0],    green); 
    p1.add_data(m3pis, mws_d[0]/10,   red); 
    p1.add_data(m3pis, mws_c[0]/10,   blue); 
    p1.set_labels("#it{m}_{3#pi}  [GeV]", "#Gamma(#it{t}, #it{m}_{3#pi}^{2})    [a.u.]");
    p1.save("intensity_m3pi_0.pdf");

    plot p2 = plotter.new_plot();
    p2.set_legend(0.2, 0.65);
    p2.set_ranges({m3pis[0], m3pis.back()}, {0, 8.7});
    p2.add_header("#minus #it{t} = 0.17 GeV^{2}");
    p2.add_data(m3pis, ews[1],    dot(jpacColor::DarkGrey, "Data"));
    p2.add_data(m3pis, mws[1],    green); 
    p2.add_data(m3pis, mws_d[1]/10,   red); 
    p2.add_data(m3pis, mws_c[1]/10,   blue); 
    p2.set_labels("#it{m}_{3#pi}  [GeV]", "#Gamma(#it{t}, #it{m}_{3#pi}^{2})    [a.u.]");
    p2.save("intensity_m3pi_1.pdf");

    plot p3 = plotter.new_plot();
    p3.set_legend(0.2, 0.65);
    p3.set_ranges({m3pis[0], m3pis.back()}, {0, 8.7});
    p3.add_header("#minus #it{t} = 0.26 GeV^{2}");
    p3.add_data(m3pis, ews[2],    dot(jpacColor::DarkGrey, "Data"));
    p3.add_data(m3pis, mws[2],    green); 
    p3.add_data(m3pis, mws_d[2]/10,   red); 
    p3.add_data(m3pis, mws_c[2]/10,   blue); 
    p3.set_labels("#it{m}_{3#pi}  [GeV]", "#Gamma(#it{t}, #it{m}_{3#pi}^{2})    [a.u.]");
    p3.save("intensity_m3pi_2.pdf");

    plot p4 = plotter.new_plot();
    p4.set_legend(0.2, 0.65);
    p4.set_ranges({m3pis[0], m3pis.back()}, {0, 8.7});
    p4.add_header("#minus #it{t} = 0.66 GeV^{2}");
    p4.add_data(m3pis, ews[3],    dot(jpacColor::DarkGrey, "Data"));
    p4.add_data(m3pis, mws[3],    green); 
    p4.add_data(m3pis, mws_d[3]/10,   red); 
    p4.add_data(m3pis, mws_c[3]/10,   blue); 
    p4.set_labels("#it{m}_{3#pi}  [GeV]", "#Gamma(#it{t}, #it{m}_{3#pi}^{2})    [a.u.]");
    p4.save("intensity_m3pi_3.pdf");
};