## BilayerAnalyzer analysis: ```apl_box``` - Area per lipid using box dimensions.
 
## Description
 
Estimate the area per lipid using the lateral area.

This analysis is used to estimate the area per lipid (APL)
using the lateral box dimensions,
A_l = 2<A_xy>/N_l ,
where A_xy is area of the lateral box dimensions (used to approximate the surface area of the bilayer)
and N_l is the number of lipids in the bilayer.
As a molar quantity this approach is only accurate for
homogenous lipid bilayers. If the bilayer is inhomogenous then this
estimate represents a composite average of the area per lipid.

This protocol is identified by the analysis key: 'apl_box'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.APLBoxProtocol'>

## Syntax

```
apl_box analysis-ID
```
* apl_box = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('apl_box apl_box_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('apl_box apl_box_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['apl_box', 'apl_box_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'apl_box', 'analysis_id': 'apl_box_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('apl_box_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('apl_box_1')
```
 
The output is type ```<class 'numpy.ndarray'>```
 
## Related analyses
* [apl_grid](apl_grid.html)

## References

1. Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical Properties of a Hydrated Lipid Bilayer
from a Multinanosecond Molecular Dynamics Simulation, Biophysical Journal, Volume 81, Issue 5, 2001,
Pages 2484-2494, ISSN 0006-3495, http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
(http://www.sciencedirect.com/science/article/pii/S0006349501758948)
