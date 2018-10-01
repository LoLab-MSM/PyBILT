
## BilayerAnalyzer analysis: ```msd_multi``` - Multiple time origin mean squared displacement.

## Description
Estimate the mean squared displacement using multiple time origins.

The MSDMultiProtocol is used to compute the mean squared displacement
(MSD) of the centers of mass of the specified lipids using multiple time
origins and thus multiple time blocks.
The MSD is given by

    MSD = <<(r(tau) - r_0)^2>\_i>\_tau

where the inner angle
brackets denote averaging over all lipids i and the
outer brackets denote averaging over all time origins tau. The diffusion
coefficient is estimated from the MSD using a simplified version
Einstein's relation

    D_i ~ MSD_i/(4.0*tau)

#### Initiated by instance of:

    <class 'pybilt.bilayer_analyzer.analysis_protocols.MSDMultiProtocol'>



## Syntax

```
msd_multi analysis-ID keyword value
```
* msd_multi = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in
        this analysis. Default: 'all', averages over all lipid types.
    * n_tau (int): Specify the time block size in number of frames.
        Default: 50 ( 1000 picoseconds for timestep of 2 fs and frame
        output ever 100 timesteps).
    * n_sigma (int): Specify the time between origins in number of frames.
        Default: 50.


### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
Add by string - use default settings:

```python
analyzer.add_analysis('msd_multi msd_multi_1')
```

Add by string - adjust leaflet setting:

```python
analyzer.add_analysis('msd_multi msd_multi_1 leaflet upper')
```

Add by list:
```python
analyzer.add_analysis(list(['msd_multi', 'msd_multi_1', dict{'leaflet':'lower', 'resname':'POPC'}]))
```
Add by dict:
```python
ba.add_analysis(dict({'analysis_key': 'msd_multi', 'analysis_id': 'msd_multi_1',
                 'analysis_settings':dict({'n_tau':100, 'n_sigma':200})}))
```
### Output Info:
Retrieve output after running analyses:
```python
msd_multi_output = analyzer.get_analysis_data('ID')
```
The output has is a numpy.ndarray with shape (4,):
* output[0] = float: The average mean squared displacement over time blocks of size tau.
* output[1] = float: The standard error of output[0].
* output[2] = float: The estimate of the diffusion coefficient.
* output[3] = float: The standard error of output[2].


| ! Note |
| --- |
| According to results in Ref 3 the time blocks used to estimate the MSD should not be overlapping. Therefore, it is recommended to use n_sigma >= n_tau. |

## Related analyses
* [msd]()

## References
1. Christofer Hofsab, Erik Lindahl, and Olle Edholm, "Molecular
    Dynamics Simulations of Phospholipid Bilayers with Cholesterol",
    Biophys J. 2003 Apr; 84(4): 2192-2206.
    doi:  10.1016/S0006-3495(03)75025-5

2. Orsi, Mario, Julien Michel, and Jonathan W. Essex.
    "Coarse-grain modelling of DMPC and DOPC lipid bilayers."
    Journal of Physics: Condensed Matter 22.15 (2010): 155106.
    http://iopscience.iop.org/article/10.1088/0953-8984/22/15/155106/meta

3. Gaurav Pranami and Monica H. Lamm, Estimating Error in Diffusion
    Coefficients Derived from Molecular Dynamics Simulations,
    Journal of Chemical Theory and Computation 2015 11 (10),
    4586-4592, DOI: 10.1021/acs.jctc.5b00574,
    http://pubs.acs.org/doi/full/10.1021/acs.jctc.5b00574
