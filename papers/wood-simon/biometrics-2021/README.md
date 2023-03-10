# Data and code for "Inferring UK COVID-19 fatal infection trajectories from daily mortality data" S.N. Wood (2020/21)

* The code all runs in R. Setting ps <- TRUE causes figures to be written to eps files. 

* Use setwd to set the working directory to the directory where the code is located. EFS-full.r needs to
  be able to locate the output from o2d.r and reported-daily.r

* o2d.r produce the draws from the combined infection to death distribution model (and paper figure 2).

* reported-daily.r contains the daily death data used, and will be sourced by EFS-full.r

* EFS-full.r contains the main analysis code. It expects o2d.r to have been run to produce a file of
  infection to death distribution model distribution parameters. To get the Flaxman model dilation check
  requires modifying the Flaxman run as described in the code comments. 

The methods should take minutes to run, not hours, although the analysis of the full run of data to early
2021 takes close to half an hour.

Other files:

* death-age.r and death-age-2020.txt are for the supporting material age effects investigation.

* PCR-test-sim.r is for the supporting material PCR incidence estimation feasibility check.

* wackamole.r is the simulation based demonstration that the lack of independence of disease
  duration and time of death is not a detail that can be ignored.
