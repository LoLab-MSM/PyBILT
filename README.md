![alt text](./_images/PyBILT_logo.png "PyBILT Logo")
# *Py*thon based lipid *BIL*ayer molecular simulation analysis *T*oolkit
------
![Python version badge](https://img.shields.io/badge/python-2.7%2C3.6%2C3.7-blue.svg)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](LICENSE)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/2d6da71328a24ef6930ad8f554074292)](https://www.codacy.com/manual/blakeaw1102/PyBILT?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=LoLab-VU/PyBILT&amp;utm_campaign=Badge_Grade)
![version](https://img.shields.io/badge/version-0.2.0-orange.svg)
[![release](https://img.shields.io/github/release-pre/LoLab-VU/PyBILT.svg)](https://github.com/LoLab-VU/PyBILT/releases/tag/v0.2.0)
[![Anaconda-Server Badge](https://anaconda.org/blakeaw/pybilt/badges/version.svg)](https://anaconda.org/blakeaw/pybilt)
[![Documentation Status](https://readthedocs.org/projects/pybilt/badge/?version=latest)](http://pybilt.readthedocs.io/en/latest/?badge=latest)

------
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

------
![alt text](./_images/7percentCL_sideview_b.jpg "Lipid Bilayer")

------

# Install

### Core dependencies
PyBILT has the following core dependencies:
   * [MDAnalysis](https://www.mdanalysis.org/)
   * [NumPy](http://www.numpy.org/)
   * [SciPy](https://www.scipy.org/)
   * [Matplotlib](https://matplotlib.org/)
   * [Seaborn](https://seaborn.pydata.org/)
   * [six](https://pypi.org/project/six/)
   * [future](http://python-future.org/)

### Python version support
The pybilt package has been tested using [Anaconda Python](https://www.anaconda.com/) 2.7, 3.6, and 3.7.

#### Sunsetting of Python 2
Please be aware that Python 2 is scheduled to be sunset on January 1 2020. You can read about it here: [https://www.python.org/doc/sunset-python-2/](https://www.python.org/doc/sunset-python-2/)
Parallel to the sunsetting of Python 2 many open source packages are also dropping support for Python 2 ([https://python3statement.org/](https://python3statement.org/)), including some of PyBILT's core dependencies. As such, after January 1, 2020, PyBILT will also likely sunset its support for Python 2.7.  

### pip install
You can install the latest version of the `pybilt` package using `pip` sourced from the GitHub repo:
```
pip install -e git+https://github.com/LoLab-VU/PyBILT@v0.2.0#egg=pybilt
```
However, this will not automatically install the core dependencies. You will have to do that separately:
```
pip install MDAnalysis numpy scipy matplotlib seaborn six future
```

### conda install
First make sure you have the `conda-forge` channel in your channel list; that is the channel from which MDAnalysis is installed. You can use the following command to add it to the bottom of your channel list:
```
conda config --append channels conda-forge
```

Then you can install the `pybilt` package from the `blakeaw` Anaconda Cloud channel,
```
conda install -c blakeaw pybilt
```
The core dependencies will be automatically installed.

### Recommended additional software

The following software is not required for the basic operation of **PyBILT**, but provides extra capabilities and features when installed.

#### pytest
The pybilt test suite is designed to be run with [pytest](https://docs.pytest.org/en/latest/), so if you want to run the tests then you will need to install pytest.

#### Jupyter
 PyBILT comes with a set of [Jupyter IPython notebooks](./jupyter_notebooks) which supplement the doc pages. If you want to run these notebooks locally then you will need to intall [Jupyter](https://jupyter.org/) (or at least the IPython kernel).

Note that the notebooks have not been updated for Python 3 yet.

#### sphinx, sphinx_rtd_theme, and recommonmark
If you want to build local versions of doc pages install the following additional packages:
  * [sphinx](http://www.sphinx-doc.org/en/master/)
  * [sphinx_rtd_theme](https://sphinx-rtd-theme.readthedocs.io/en/stable/)
  * [recommonmark](https://recommonmark.readthedocs.io/en/latest/)


------

# Documentation and Usage

## Quick overview of PyBILT
**PyBILT** is composed of 2 primary analysis packages:
  * bilayer_analyzer -- The [bilayer_analyzer](http://pybilt.readthedocs.io/en/latest/pybilt.bilayer_analyzer.html#module-pybilt.bilayer_analyzer.bilayer_analyzer) is an analysis package that
                        is designed to analyze (quasi) planar lipid bilayer
                        systems. It is accessed through the BilayerAnalyzer
                        object, which can be imported via: ```from
                        pybilt.bilayer_analyzer import BilayerAnalyzer```. The
                        BilayerAnalyzer features automatic dynamic unwrapping of
                        coordinates and leaflet detection. The bilayer_analyzer
                        works on a multiple-representation model, whereby the
                        various analyses are conducted using different
                        representations of the bilayer lipids. Bilayer lipids
                        can be represented using the following four
                        representations:
    * All atom
    * Centers-of-mass -- Each lipid (or selection of atoms from each lipid) is reduced to a
center-of-mass.
    * Grid (or lipid grid) -- The lipids are mapped to two-dimensional grids (one for each leaflet) in the
style of the [GridMAT-MD method](http://www.bevanlab.biochem.vt.edu/GridMAT-MD/)
    * Vectors - Each lipid is converted to a vector representation using select reference atoms (or sets of reference atoms) that are used to compute the head and tail of the vector; e.g., a lipid tail atom to lipid head atom, or P-N vectors.

The bilayer_analyzer features various types of analyses and the use of different
representations is handled internally based the requirements and design of each
analysis type. See the [documentation](https://pybilt.readthedocs.io/en/latest/ba_analyses.html) for list of analyses that can be added to intances of the BilayerAnalyzer.   

  * [mda_tools](http://pybilt.readthedocs.io/en/latest/pybilt.mda_tools.html) -- This package includes various modules and functions for directly
                 analyzing and operating on MDAnalysis trajectories and objects.
                 e.g. functions to compute density profiles.

 Additional packages include:
   * [lipid_grid](http://pybilt.readthedocs.io/en/latest/pybilt.lipid_grid.html) -- The lipid grid module can be used construct "lipid grid" grid
                  representations of lipid bilayers, which can be used to
                  accurately estimate quantities such as area per lipid.

  * [com_trajectory](http://pybilt.readthedocs.io/en/latest/pybilt.com_trajectory.html) -- This module can be used to construct a center of mass
                      trajectory (COMTraj) out of an MDAnalysis trajectory,
                      which is useful for computing quantities like mean squared
                      displacement. The COMTraj is designed to work with bilayers.

  * [plot_generation](http://pybilt.readthedocs.io/en/latest/pybilt.plot_generation.html) -- This module has several pre-written plotting functions
                       (using matplotlib and seaborn) for some of the properties
                       that can be computed from functions in the other modules.
                       e.g. mean squared displacement and area per lipid.
## Docs
Visit the PyBILT docs on [Read the Docs](http://pybilt.readthedocs.io/en/latest/index.html).
Docs can also be viewed offline/locally by opening the [PyBILT/docs/build/html/index.html](docs/build/html/index.html) file from the
repo in a web browser; however, this build of the docs is not updated often.

## Jupyter IPython notebooks
In addition
to the Docs, there are currently a few Jupyter IPython
[notebooks](jupyter_notebooks) that provide some examples and show some basic
usage (these have not been updated/tested for/with python 3 yet); updates and more of these are in the pipeline.

------

# Contact

To report problems or bugs please open a
[GitHub Issue](https://github.com/LoLab-VU/PyBILT/issues). Additionally, any
comments, suggestions, or feature requests for PyBILT can also be submitted as
a
[GitHub Issue](https://github.com/LoLab-VU/PyBILT/issues).

------

# Contributing

If you would like to contribute directly to PyBILT's development please
 1. Fork the repo (https://github.com/LoLab-VU/PyBILT/fork)
 2. Create a new branch for your feature (git checkout -b feature/foo_bar)
 3. Create test code for your feature
 4. Once your feature passes its own test, run all the tests using [pytest](https://docs.pytest.org/en/latest/) (python -m pytest)
 5. Once your feature passes all the tests, commit your changes (git commit -am 'Add the foo_bar feature.')
 6. Push to the branch (git push origin feature/foo_bar)
 7. Create a new Pull Request

------

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

------

# Acknowledgments

* A special thanks to James Pino (https://github.com/JamesPino) for his inciteful comments and suggestions that have helped improve the quality of this code, and thanks to him for pointing out some very useful coding tools.   
* Thanks to my advisors, Carlos F. Lopez and Arvind Ramanathan, for catalyzing this project and for providing me with the space and means to pursue it.  

------

# Citing

If you use the **PyBILT** software as a part of your research, please cite the GitHub repo.

Also, please cite the following references as appropriate for scientific/research software used with/via **PyBILT**:

#### MDAnalysis
See: https://www.mdanalysis.org/pages/citations/

#### Packages from the SciPy ecosystem

These include NumPy, SciPy, and Matplotlib for which references can be obtained from:
https://www.scipy.org/citing.html

#### seaborn
Reference can be exported from the [seaborn Zeondo DOI entry](https://doi.org/10.5281/zenodo.592845)

#### Jupyter
See: https://github.com/jupyter/jupyter/issues/190

#### IPython

See: https://ipython.org/citing.html

------
