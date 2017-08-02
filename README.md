# PyBILT         ![docstring-coverage badge](https://img.shields.io/badge/docstring%20coverage-31.4%25-yellow.svg)
## *Py*thon based lipid *BIL*ayer analysis *T*oolkit

![alt text](./_images/7percentCL_sideview_b.jpg "Lipid Bilayer")
#### PyBILT is a Python toolkit developed to analyze molecular simulation trajectories of lipid bilayers systems. The toolkit includes a variety of analyses from various lipid bilayer molecular simulation publications.

The analyses include:
   * Mean Squared Displacement (MSD)
   * Diffusion coefficent estimators (from MSD curves) - includes Einstein relation, linear fit, and anomalous diffusion fit.
   * Area per lipid estimators
   * Bilayer thickness
   * Displacement Vector (flow) maps and correlations
   * Deuterium order parameter
   * Orientation parameters
   * Mass and Electron Density Estimators
   * and more!

### Install

#### Warning: PyBILT is still under heavy development and may rapidly change.

Download PyBILT from the github repo (https://github.com/blakeaw/PyBILT.git)
and then add the path to the PyBILT directory to your PYTHONPATH. From a
terminal you can type
```
export PYTHONPATH="path_to/PyBILT:$PYTHONPATH
```
to add it to the current shell environment. For persistence add the line to your
.bashrc file.

PyBILT has the following major dependencies:
   * MDAnalysis 0.16.0
   * NumPy  1.11.3
   * SciPy 0.18.1,
   * Matplotlib 2.0.0
   * Seaborn 0.7.1


In addition, it is highly recommended that you install
[Anaconda Python](https://www.continuum.io/)
version 4.3.1 Python 2.7 before installing PyBILT. PyBILT
has yet to be tested outside of an Anaconda environment.

#### Setup using Anaconda's conda tool
The file environment.yml has been provided to allow for easy setup of a new
environment with all the appropriate dependencies using the conda tool. Run
```
conda env create -f environment.yml
```
which will create a new conda environment named *pybilt* with the appropriate
dependencies. Then activate the environment
```
source activate pybilt
```
before running PyBILT modules.

#### Quick overview of PyBILT
**PyBILT** is composed of 5 major modules:
  * bilayer_analyzer -- The bilayer_analyzer is the primary analysis module.
                        This object has the most comprehensive set of built-in
                        features (types of analyses and use of different bilayer
                        representations) and makes use of much of the
                        functionality from proceeding modules.  
  * mda_tools -- This module has various sets of functions for directly
                 analyzing and operating on MDAnalysis trajectories and objects.
                 e.g. functions to compute density profiles.
  * lipid_grid -- The lipid grid module can be used construct "lipid grid" grid
                  representations of lipid bilayers, which can be used to
                  accurately estimate quantities such as area per lipid.
  * com_trajectory -- This module can be used to construct a center of mass
                      trajectory (COMTraj) out of an MDAnalysis trajectory,
                      which is useful for computing quantities like mean squared
                      displacement.
  * plot_generation -- This module has several pre-written plotting functions
                       (using matplotlib and seaborn) for some of the properties
                       that can be computed from functions in the other modules.
                       e.g. mean squared displacement and area per lipid.
