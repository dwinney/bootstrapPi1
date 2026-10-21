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

    auto xpars = import_data<6>(path+"/data/pars/extpars.dat");
    std::vector<double> xNc, xNcp, xNd, xphi;
    for (int i = 0; i < xpars[0].size(); i++)
    {
        xNc.push_back( xpars[2][i]);
        xNcp.push_back(xpars[3][i]);
        complex Nd = xpars[4][i] + I*xpars[5][i];
        xNd.push_back( std::abs(Nd) );
        double xp =  std::arg(Nd);
        if (xp > 0) xphi.push_back(xp-2*PI);
        else xphi.push_back(xp);
    };

    // -----------------------------------------------------------------------
    // Plot results
    
    plotter plotter;

    // Plot distributions of parameters
    plot p1 = plotter.new_plot();
    p1.set_legend(0.65, 0.65);
    p1.set_labels("#it{m}_{3#pi}   [GeV]", "|#it{N}| / 10^{3}");
    p1.add_data( Nc[0],     {Nc[1]/1E3, Nc[2]/1E3},  solid(jpacColor::Blue));
    p1.add_data( Nc[0],     xNc/1E3,    dot(jpacColor::Blue, "Contact"));
    p1.add_data( Nd[0],     {Nd[1]/1E3, Nd[2]/1E3},   solid(jpacColor::Red));
    p1.add_data( Nc[0],     xNd/1E3,   dot(jpacColor::Red, "Deck"));
    p1.save("Ns.pdf");

    // Plot distributions of parameters
    plot p3 = plotter.new_plot();
    p3.set_legend(0.6, 0.65);
    p3.set_labels("#it{m}_{3#pi}   [GeV]", "|#it{N}| / 10^{3}");
    p3.add_data( Nc[0],      xNcp/1E3,   dot(jpacColor::Orange, "Second contact term"));
    p3.save("Nsp.pdf");

    plot p2 = plotter.new_plot();
    p2.add_horizontal(-PI, {2011, kDashed});
    p2.set_labels("#it{m}_{3#pi}   [GeV]", "#it{#phi}#kern[-0.1]{_{#it{d}}} - #it{#phi}#kern[-0.1]{_{#it{c}}}");
    p2.add_data(phi[0], {phi[1], phi[2]},  solid(jpacColor::Green));
    p2.add_data(phi[0], xphi,   dot(jpacColor::Green));
    p2.save("phi.pdf");
};