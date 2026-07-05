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

void plot_width_t()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // Which tbins to plot
    std::array<int,3> m3pibins = {14, 24, 34};

    // -----------------------------------------------------------------------
    // Data set up

    // Import our data sets  
    std::vector<double> ets = COMPASS::mt_bins(), twidths;
    std::array<std::vector<data_set>,3> data;
    for (int j = 0; j < 3; j++)
    {
        for (int i = 0; i < 4; i++)
        {
            data[j].emplace_back(COMPASS::parse_JSON(m3pibins[j], i));
            if (j==0) twidths.push_back(data[j].back()._extras["tbin_width"]);
        };
    };

    // -----------------------------------------------------------------------
    // Set up amplitude and iterative solution

    std::array<std::vector<double>,3> ews;
    for (int j = 0; j < 3; j++)
    {
        for (auto bin : data[j])
        {
            double ew = 0;
            for (int i = 0; i < bin._z.size(); i++) ew += norm(bin._z[i])*bin._dx[i];
            ews[j].push_back(ew);
        };
    };

    std::string path =  gSystem->Getenv("BOOTSTRAP_SRC");
    auto x = import_data<4>(path+"/data/intensity_14.dat");
    auto y = import_data<4>(path+"/data/intensity_24.dat");
    auto z = import_data<4>(path+"/data/intensity_34.dat");

    // -----------------------------------------------------------------------
    // Plot results

    plotter plotter;

    // Plot widths as a function of m3pi
    plot p1 = plotter.new_plot();
    p1.set_legend(0.25, 0.1);
    p1.set_logscale(false, true);
    p1.add_header("#it{m}_{3#pi} = " + to_string(COMPASS::m_bins[m3pibins[0]-11]) + " GeV");
    p1.set_labels("#minus #it{t}  [GeV^{2}]", "Integrated Intensity  [a.u.]");
    p1.add_curve(x[0], x[2], "Contact");
    p1.add_curve(x[0], x[3], "Deck");
    p1.add_curve(x[0], x[1], "Full");
    p1.add_data ({ets, twidths},  ews[0]/1E4,  dot(jpacColor::DarkGrey, "Data"));


    plot p2 = plotter.new_plot();
    p2.set_legend(0.25, 0.075);
    p2.set_logscale(false, true);
    p2.add_header("#it{m}_{3#pi} = " + to_string(COMPASS::m_bins[m3pibins[1]-11]) + " GeV");
    p2.set_labels("#minus #it{t}  [GeV^{2}]", "Integrated Intensity  [a.u.]");
    p2.add_curve(y[0], y[2]);
    p2.add_curve(y[0], y[3]);
    p2.add_curve(y[0], y[1]);
    p2.add_data ({ets, twidths},  ews[1]/1E4,  dot(jpacColor::DarkGrey));

    plot p3 = plotter.new_plot();
    p3.set_legend(0.25, 0.20);
    p3.set_logscale(false, true);
    p3.add_header("#it{m}_{3#pi} = " + to_string(COMPASS::m_bins[m3pibins[2]-11]) + " GeV");
    p3.set_labels("#minus #it{t}  [GeV^{2}]", "Integrated Intensity  [a.u.]");
    p3.add_curve(z[0], z[2]);
    p3.add_curve(z[0], z[3]);
    p3.add_curve(z[0], z[1]);
    p3.add_data ({ets, twidths},  ews[2]/1E4,  dot(jpacColor::DarkGrey));

    plotter.stack({p1,p2,p3}, "widths.pdf");   
};