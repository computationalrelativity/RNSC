# RNS

Library to create initial data for rotating neutron stars using Stergioulas RNS code and import them to Cartesian grids for 3+1 evolutions.

The original code is based on the method developed by [Komatsu, Eriguchi & Hachisu (1989)](https://academic.oup.com/mnras/article/237/2/355/976460) and modifications introduced by [Cook, Shapiro & Teukolsky (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...422..227C/abstract) and detailed in [Stergioulas and Friedman (1995)](https://arxiv.org/abs/astro-ph/9411032) and [Nozawa, Stergioulas, Gourgoulhon & Eriguchi (1998)](https://arxiv.org/abs/gr-qc/9804048).

The code current is an early version by SB of the `Cactus/WhiskyRNSID` module that finished in the [Cactus/Einstein Toolkit](https://einsteintoolkit.org/thornguide/EinsteinInitialData/Hydro_RNSID/documentation.html).

WARNING

 * Error handling is very poor and problematic with multi-threading or MPI.
 * The use of EOS tables should be tested in some detail.

HISTORY

 * SB 07/2022 started

