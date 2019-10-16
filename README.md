# N2-FIXATION-NATURE
doit.m       is the main driver

neglogpost.m computes the negative of the logarithm of the posterior
             pdf for the parameters conditioned on the DIN, DIP and
	     DON data 
	     
buildPFD.m   computes the sinking particle flux divergence operator

uptake.m     computes the rate of uptake of DIN given the DIP uptake
             as well as the 1st and 2nd derivatives w.r.t. the model
	     parameters
	     
denit_don.m  computes the denitrification rate as well as the 1st and
             2nd derivative w.r.t. the model parameters
	     
fixit.m      computes the N uptake limiter that implies N fixation,the
             derivative of the limiter with respect to the N model
	     state as well as the 1st and 2nd derivatives with respect
	     to the model parameters
	     
eqPcycle.m   computes the equilibrium state of the P-cycle model as
             well as the 1st and 2nd derivaties of the equilibrium
	     solution  w.r.t. the model parameters
	     
eqNcycle_v2.m computes the equilibrium of the N-cycle model state as
              well as the 1st and 2nd derivative of the solution w.r.t
	      the model parameters

doit.m 
|
+--> neglogpost.m
       |
       +--> eqPcycle.m
       |       |
       |       +-->[buildPFD.m]
       +--> eqNcycle_v2.m
             |
	     +-->[buildPFD.m uptake.m denit_don.m fixit.m]


utility scripts:

nsnew.m   C.T. Kelly's Newton Armijo solver

mfactor.m Timothy A. Davis' LINFACTOR VERSION 1.1.0, Nov 1, 2007
          Copyright 2007, Timothy A. Davis, University of Florida
	  
d0.m      makes a sparse diagonal matrix given a 3d field
	
data needed to drive the model:

The tracer transport operators and binned WOA13 and/or GLODAP data sets used for the optimization are available from (fprimeau@uci.edu).
