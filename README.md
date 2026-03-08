# bootstrapPi1
Statically-compiled implementation of the $\pi_1 \to 3\pi$ decay fitting script amenable for uncertainty analysis through bootstrapping. We require the [iterateKT](https://github.com/dwinney/iterateKT) library and thus refer to the dependencies and installation instructions there. 

This library compiles iterateKT as a static library and a single executable to run a single global fit of the unitarized contact-plus-Deck model of [1] to pseudo-data generated from the COMPASS freed-isobar data set [2]. Additional [scripts](./bash) are included to streamline submitting arbitrary instances of this executable to a high-performance computing cluster using SLUM (specifically, the TOCHTLI cluster at ICN-UNAM). 

## References
- [1] [Production and Rescattering in $\pi_1 \to 3\pi$](https://youtu.be/-f71ClZ00ts?si=P3TTd5_tJO23zXzj)
- [2] [Exotic meson $\pi_1(1600)$ with $J^{PC}=1^{-+}$ and its decay into $\rho(770)\pi$](https://inspirehep.net/literature/1898933)
