## BilayerAnalyzer analysis: ```msd``` - Single time origin mean squared displacement.
 
## Description
 
Estimate the mean squared displacement from a single time origin.

This analysis computes the mean squared displacement (MSD) versus time
of the centers of mass of the specified lipids for a single time origin
(i.e. the first frame in the trajectory analysis):
MSD(t) = <(r(t) - r_0)^2>\_i
where the angle brackets denote averaging over all lipids of type i;
the lipid type used for the analysis is controlled by the ```leaflet```
setting.

This protocol is identified by the analysis key: 'msd'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.MSDProtocol'>

## Syntax

```
msd analysis-ID keyword value
```
* msd = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', averages over all lipid types.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('msd msd_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('msd msd_1 leaflet both')
```
 
Add by list:
```python
analyzer.add_analysis(list(['msd', 'msd_1', dict({'leaflet':'both'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'msd', 'analysis_id': 'msd_1','analysis_settings':dict({'leaflet':'both'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('msd_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('msd_1')
```
 
The output is type ```<class 'numpy.ndarray'>```
 
## Related analyses
* [msd_multi](msd_multi.html)

## References

1. Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical
Properties of a Hydrated Lipid Bilayer from a Multinanosecond
Molecular Dynamics Simulation, Biophysical Journal, Volume 81,
Issue 5, 2001, Pages 2484-2494, ISSN 0006-3495,
http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
(http://www.sciencedirect.com/science/article/pii/S0006349501758948)

2. Yoshimichi Andoh, Susumu Okazaki, Ryuichi Ueoka, Molecular
dynamics study of lipid bilayers modeling the plasma membranes
of normal murine thymocytes and leukemic GRSL cells, Biochimica
et Biophysica Acta (BBA) - Biomembranes, Volume 1828, Issue 4,
April 2013, Pages 1259-1270, ISSN 0005-2736,
https://doi.org/10.1016/j.bbamem.2013.01.005.
(http://www.sciencedirect.com/science/article/pii/S0005273613000096)

3. Section 8.7,
http://manual.gromacs.org/documentation/5.1.4/manual-5.1.4.pdf
