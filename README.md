![alt text](./_images/PyBILT_logo.png "PyBILT Logo")
# *Py*thon based lipid *BIL*ayer molecular simulation analysis *T*oolkit
------
![Python version badge](https://img.shields.io/badge/python-2.7%2C3.6-blue.svg)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](LICENSE)
[![Code Health](https://landscape.io/github/LoLab-VU/PyBILT/master/landscape.svg?style=flat)](https://landscape.io/github/Lolab-VU/PyBILT/master)
[![Documentation Status](https://readthedocs.org/projects/pybilt/badge/?version=latest)](http://pybilt.readthedocs.io/en/latest/?badge=latest)
[![docstring-coverage badge](https://img.shields.io/badge/docstring--coverage-49.5%25-orange.svg)](.docstring-coverage_report.txt)
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

| **! Warning** |
| :--- |
|  PyBILT is still under heavy development and may rapidly change. |

#### PyBILT run dependencies
PyBILT has the following major dependencies:
   * MDAnalysis - https://www.mdanalysis.org/
   * NumPy - http://www.numpy.org/
   * SciPy - https://www.scipy.org/
   * Matplotlib - https://matplotlib.org/
   * Seaborn - https://seaborn.pydata.org/
   * six - https://pypi.org/project/six/
   * future - http://python-future.org/

To run the pybilt test suite:
   * pytest - https://docs.pytest.org/en/latest/

To use the PyBILT Jupyter Notebooks:
   * Jupyter - http://jupyter.org/

To build the docs locally requires the following additional packages:
   * Sphinx - http://www.sphinx-doc.org/en/master/
   * recommonmark - https://recommonmark.readthedocs.io/en/latest/
   * sphinx_rtd_theme - https://sphinx-rtd-theme.readthedocs.io/en/latest/



The pybilt package has been tested using [Anaconda Python](https://www.anaconda.com/) 2.7 and Python 3.6;

The following section describes the process for setting up the dependencies and
installing the 'pybilt' package using a conda environment and the setup.py
script.

## Setup and install using Anaconda's conda tool

### Method 1 (manual package installation)
First, clone or download the GitHub repo
```
git clone https://github.com/blakeaw/PyBILT.git
```
Then create a new conda environment for pybilt and activate it:
```
conda create --name pybilt
conda activate pybilt
```
The install the preferred python version:
 * for Python 2.7
```
conda install python=2.7
```
 * for Python 3.6
 ```
 conda install python=3.6
 ```

Then install all the pybilt run dependencies:
```
conda install numpy scipy matplotlib seaborn six future
conda install -c conda-forge MDAnalysis
```
Then install pybilt:
```
python PyBILT/setup.py install
```

If you want to run the pybilt tests you can install pytest:
```
conda install pytest
```

If you want to run the pybilt Jupyter notebooks (PyBILT/jupyter_notebooks), then install Jupyter:
```
conda install jupyter
```
Note that the notebooks have not been updated for Python 3 yet.

If you want to build local versions of doc pages install the following packages:
```
conda install sphinx
pip install sphinx_rtd_theme
pip install recommonmark
```

### Method 2 (From environment yaml)

The files environment_py27.yml and environment_py36.yml have been provided to allow for easy setup of a new conda
environment with all of the most recently tested versions of dependencies.
Run:
 * For Python 2.7:
```
conda env create -f environment_py27.yml
```
 * For Python 3.6
 ```
 conda env create -f environment_py36.yml
 ```

which will create a new conda environment named *pybilt* with the appropriate
dependencies. Then activate the environment
```
conda activate pybilt
```
Next, run the setup.py script with install,
```
python PyBILT/setup.py install
```
to install the 'pybilt' package into the *pybilt* environment.  

If you want to run the pybilt Jupyter notebooks (PyBILT/jupyter_notebooks), then install Jupyter:
```
conda install jupyter
```
Note that the notebooks have not been updated for Python 3 yet.

------

# Quick overview of PyBILT
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

------

# Documentation

Visit the PyBILT docs on [Read the Docs](http://pybilt.readthedocs.io/en/latest/index.html).
Docs can also be viewed offline/locally by opening the [PyBILT/docs/build/html/index.html](docs/build/html/index.html) file from the
repo in a web browser; however, this build of the docs is not updated often.
In addition
to the doc pages, there are currently a few Jupyter IPython
[notebooks](jupyter_notebooks) that provide some examples and show some basic
usage (these have not been updated/tested for/with python 3 yet); updates and more of these are in the pipeline.

------

# Contact

To report problems or bugs please open a
[GitHub Issue](https://github.com/blakeaw/PyBILT/issues). Additionally, any
comments, suggestions, or feature requests for PyBILT can also be submitted as
a
[GitHub Issue](https://github.com/blakeaw/PyBILT/issues).

For any other inquiries, including questions about PyBILT use or
implementation, you can contact Blake directly via e-mail at either
blake.a.wilson@vanderbilt.edu or blakeaw1102@gmail.com; please include "PyBILT
inquiry" in the e-mail subject line.

------

# Contributing

If you would like to contribute directly to PyBILT's development please
 1. Fork the repo (https://github.com/blakeaw/PyBILT/fork)
 2. Create a new branch for your feature (git checkout -b feature/foo_bar)
 3. Create test code for your feature
 4. Once your feature passes its own test, run all the tests using [pytest](https://docs.pytest.org/en/latest/) (python -m pytest)
 5. Once your feature passes all the tests, commit your changes (git commit -am 'Add the foo_bar feature.')
 6. Push to the branch (git push origin feature/foo_bar)
 7. Create a new Pull Request

------

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

------

## Acknowledgments

* A special thanks to James Pino (https://github.com/JamesPino) for his inciteful comments and suggestions that have helped improve the quality of this code, and thanks to him for pointing out some very useful coding tools.   
* Thanks to my advisors, Carlos F. Lopez and Arvind Ramanathan, for catalyzing this project and for providing me with the space and means to pursue it.  

------

## Built With

* [ANACONDA](https://www.continuum.io/) - ANACONDA Python distribution and CONDA package and environment manager
* [PyCharm](https://www.jetbrains.com/pycharm/) - Text Editor/IDE
* [ATOM](https://atom.io/) - Text Editor/IDE
* [Sublime Text](https://www.sublimetext.com/) - Text Editor used in earlier work
* [Landscape](https://landscape.io/) - Code quality analysis and tracking
* [Git](https://git-scm.com/) - Version control
* [GitHub](https://github.com/) - Development Platform and repository storage
* [Sphinx](http://www.sphinx-doc.org/en/stable/index.html) - Python documentation generator
* [recommonmark](https://github.com/rtfd/recommonmark) - A docutils-compatibility bridge to CommonMark
* [Read the Docs](https://readthedocs.org/) - Documentation hosting
* [docstring-coverage](https://bitbucket.org/DataGreed/docstring-coverage/wiki/Home) -  A simple audit tool for examining python source files for missing docstrings.
* [Python-Modernize](https://python-modernize.readthedocs.io/en/latest/) - Automatic modernization of Python 2 code for dual Python 2 and 3 support.

------

## Core Developers

* **Blake A Wilson** - Currently a Postdoctoral Fellow at Vanderbilt University
  * Vandy e-mail: blake.a.wilson@vanderbilt.edu
  * Gmail: blakeaw1102@gmail.com
  * [Blake's VU Website]( https://my.vanderbilt.edu/blakeaw/)
  * Also find me on [LinkedIn](https://www.linkedin.com/in/blakewilson3/) and [Research Gate](https://www.researchgate.net/profile/Blake_Wilson3)

------
