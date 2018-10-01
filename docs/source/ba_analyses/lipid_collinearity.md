## BilayerAnalyzer analysis: ```lipid_collinearity``` - Lipid-lipid collinearity.
 
## Description
 
Estimate the lipid-lipid collinearity angles.

This analysis computes the mean lipid-lipid collinearity angle (or order
parameter) using the vector represetation of the specified lipids.

This protocol is identified by the analysis key: 'lipid_collinearity'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.LipidCollinearityProtocol'>

## Syntax

```
lipid_collinearity analysis-ID keyword value
```
* lipid_collinearity = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'upper'
    * resname_1 (str): Specify the resname of the reference lipid type to include in this analysis. Default: 'first', the first lipid in the list pulled from the com_frame representation.
    * resname_2 (str): Specify the resname of the target lipid type to include in this analysis. Default: 'first', the first lipid in the list pulled from the com_frame representation.
    * style (str: 'angle', 'order'): Specify whether to compute the tilt angle ('angle') or the tilt angle order parameter ('order'). Default: 'angle'

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('lipid_collinearity lipid_collinearity_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('lipid_collinearity lipid_collinearity_1 leaflet upper')
```
 
Add by list:
```python
analyzer.add_analysis(list(['lipid_collinearity', 'lipid_collinearity_1', dict({'leaflet':'upper'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'lipid_collinearity', 'analysis_id': 'lipid_collinearity_1','analysis_settings':dict({'leaflet':'upper'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('lipid_collinearity_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('lipid_collinearity_1')
```
 
The output is type ```<class 'numpy.ndarray'>```
 
## Related analyses
* [lipid_tilt](lipid_tilt.html)

## References

1. Anton O. Chugunov,  Pavel E. Volynsky, Nikolay A. Krylov,
Ivan A. Boldyrev, and Roman G. Efremov,  Liquid but Durable:
Molecular Dynamics Simulations Explain the Unique Properties of
Archaeal-Like Membranes, Scientific Reports, 4:7462, 2014,
doi:10.1038/srep07462
(https://www.nature.com/articles/srep07462)

