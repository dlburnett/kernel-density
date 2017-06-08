# Kernel-density analysis
Kernel density analysis for each farm; based on simulated data from the particle tracking model. 

File names and what they are:
1. Kernel density code from Raph.R
      This script does a few things....
      i.) This will take a csv for the time frame you are analyzing (i.e. one week of pulses released, and followed for the entire 11 days        of tracking), then filter them down to whichever farm you are doing the kernel density analysis for (i.e. farm 1-20). Note that           another script (I THINK IT IS COMPILE_CSV.R BUT DOUBLE CHECK IT)  is the one that will turn the netcdf files into csvs, and then           compile them into one csv. 
      ii.) Filter particles so we are only doing the density for particles that are still alive and infectious
      iii.) Load shape files of the coast line to be used a boundary for the kernel density analysis 
      iv.) Compute a bandwidth, which is based on the "guideline of 8". This guideline makes a rectangle that encompasses all points from       the pulses, and then divides the longest side by 8. This is a good rule of thumb to establish an initial bandwidth. Note that code         can eaisly be changed to accomodate different bandwidths.  
      v.) Generate's the kernel density surface one farm at a time and using Diggle's edge correction
      vi.) Extract the kernel density value FOR each farm FROM each farm to generate a connectivity matrix 
