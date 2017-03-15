# ORBILT - Oak Ridge BILayer analysis Toolkit

### Mean Squared displacement functions:
  * pMSD.py - a set of functions to compute the MSD of MD Analysis selections and lists of selections
  * pMSD_cyton.pyx - a Cython version of the MSD_list function from pMSD.py
  * mda_msd.py - combines most of the features of the individual functions in pMSD.py into a single function

### Profile Desnity functions
  * mda_density_profile.py - A set of fucntions to compute common profile density measurements
      on an MDAnalysis trajectory. Includes a fucntion to compute mass profile density, as well
      two estimators of the electron profile density. 
 

### Center of Mass (COM) representation functions and classes:
  * MemSys.py - processes the trajectory and selection of the membrane lipids and
	      rebuilds the trajectory with just the centers of mass of the lipids
	      Includes various anlyses based on the COM mass motion of the lipids.	
  *	Version 2 - The frames of COM representations are stored on disk using the shelve module.
	            This allows large trajectories to be processed without eating up all the RAM 
                    (although this requires up to a few GB of free disk space).
### Lipid Gridding and Analysis
  * mem_sys_leaflet_gridding.py - reads in a COM frame (from MemSys) and converts it to 
        assigns the lipids to a grid. The grid can then be analyzed in order to estimate
        properties such as area per lipid. The gridding and anlaysis procedures are based on 
        the decriptions given in Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858,
        which is itself a modified version of the GridMAT-MD method by 
        Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry.
				
# Dependencies:
   * NumPy,
   * SciPy,
   * Matplotlib,
   * Seaborn,
   * MDAnalysis
