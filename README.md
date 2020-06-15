# gm2-EDMsim

Python simulation code for g-2 EDM studies, without the faff of running the official code. 

# Main options: 
1. Option: main running mode. Two options:
      - MCgen: generates events according to an exponential time distribution, bins counts for both g-2 wiggle with an energy cut and average angle modulo the g-2 period for EDM analysis.
      - test: used for a single time snapshot, good for debugging
2. n_events: the total number of events to generate. 
3. t_start and t_end: define the time range over which to generate events. 
4. n_bins: the number of time bins to use. 

The n_events and n_bins have the usual tradeoff for good statistics per bin. 

# WiggleFitter.py

Run this on the output text files generated by gm2Sim to fit the wiggle plot or vertical angle oscillation. Calculates parameters, their errors, and the chi-squared of the fit. If the input files have no injected EDM signal, will also calculate the limit set on the EDM by the generate data, if there is an EDM the results are nonsensical at the moment (planning to implement a check later). 

# hadd.py

Copying the ROOT 'hadd' method for combining historgrams and datasets to allow gm2sim's output txt files to be combined sensibly. Only works on gm2sim output patterns atm.

# gm2.mplstyle 

This is the official matplotlib style file for the g-2 experiment, needed to make the plots look nice. 

# CLs.py and FC.py

Two methods for setting more rigorious EDM limits than the vanilla confidence interval.


