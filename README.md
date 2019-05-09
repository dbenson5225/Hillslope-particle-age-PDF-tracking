# Hillslope-particle-age-PDF-tracking
Contains files to create 2-D heterogeneous K fields, solve velocity for idealized hillslope, and track particles with age PDFs through aquifer.
Run 3 matlab files in sequence: Bridge2d_exp.m, vcalc_sparse.m, and het_v_aniso.m.
The first will create one or more heterogeneous K fields.  You may condition this field, and control the heterogeneity statistics.  Right now it creates exponential correlation function (power law has been disabled). The K is output as kout_exp_#, where # is the number of the realization.
The second solves a BVP for constant infiltration along the top of the previously created K field (after reading and loading it).  It uses a constant head 0f zero to represent the stream on the right.  Other boundaries are no-flow.  The K, v fields and other parameters that you control (Infiltration, porosity) are written to disk in "Velocities_exp.mat"
The third file reads these velocities and tracks particles through the aquifer.  It keeps track of an age PDF on each particle and mixes the ages between particles.
Search for the paper or preprint entitled "Aging and mixing as pseudo-chemical-reactions between, and on, particles: Perspectives on particle interaction and multi-modal ages in hillslopes and streams"  for more details.
