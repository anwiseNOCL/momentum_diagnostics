# momentum_diagnostics

Context
* At some point the momentum diagnostics that could be output from NEMO simulations stopped being updated to align with changes in the underlying calculations of the momentum equations. David Storkey created a branch that corrected this at NEMO version 4.0.4 https://forge.ipsl.jussieu.fr/nemo/wiki/2020WP/ENHANCE-01_davestorkey_fix_3D_momentum_trends
* These changes have been ported to version 4.2.0 https://forge.nemo-ocean.eu/nemo/nemo/-/tree/RP_momentum_trend_diag_fixes/src/OCE?ref_type=heads
* This repo is a port to NEMO version 4.2.1 with one further modification (see below)
* When the model is run with an inverse barometer effect due to spatially varying atmospheric pressure then this can be seperated out in the momentum diagnostics. In this version of the code this term is included as its own term in the 2d budget.
* Note that all changes relative to version 4.2.0 are commented with a prefix "AW"
