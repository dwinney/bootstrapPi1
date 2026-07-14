// Take results from a multi-dimensional fit and output distributions of parameters
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2026)
// Affiliation:  Instituto de Ciencias Nucleares (ICN)
//               Universidad Nacional Autonoma de Mexico (UNAM)
// Email:        daniel.winney@nucleares.unam.mx
// ------------------------------------------------------------------------------

#include "constants.hpp"
#include "plotter.hpp"

void plot_pars()
{
    using namespace iterateKT;
    using iterateKT::complex;
    using iterateKT::to_string;

    // -----------------------------------------------------------------------
    // Operating options

    // File containing parameters
    std::string path =  gSystem->Getenv("BOOTSTRAP_SRC");

    // -----------------------------------------------------------------------
    // Import fit values

    auto Nc   = import_data<3>(path+"/data/pars/modNc.dat");
    auto Nd   = import_data<3>(path+"/data/pars/modNd.dat");
    auto phi  = import_data<3>(path+"/data/pars/argNd.dat");

    // -----------------------------------------------------------------------
    // Plot results
    
    plotter plotter;

    // Plot distributions of parameters
    plot p1 = plotter.new_plot();
    p1.set_legend(0.675, 0.65);
    p1.set_labels("#it{m}_{3#pi}   [GeV]", "|#it{N}| / 10^{3}");
    p1.add_data( Nc[0],     {Nc[1]/1E3, Nc[2]/1E3},   dot(jpacColor::Blue, "Contact"));
    p1.add_data( Nd[0],     {Nd[1]/1E3, Nd[2]/1E3},   dot(jpacColor::Red,  "Deck"));
    p1.save("Ns.pdf");

    plot p2 = plotter.new_plot();
    p2.add_horizontal(-PI, {2011, kDashed});
    p2.set_labels("#it{m}_{3#pi}   [GeV]", "#it{#phi}#kern[-0.1]{_{#it{d}}} - #it{#phi}#kern[-0.1]{_{#it{c}}}");
    p2.add_data(phi[0], {phi[1], phi[2]},  dot(jpacColor::Green));
    p2.save("phi.pdf");
};