## BilayerAnalyzer analysis: ```lipid_tilt``` - Lipid tilt using the lipid vectors.
 
## Description
 
Estimate the lipids tilt angles using the defined lipid vector.

This analyis estimates the mean lipid tilt using the vector represetation
of the specified lipids in reference to a particular axis, typically the
bilayer normal.

This protocol is identified by the analysis key: 'lipid_tilt'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.LipidTiltProtocol'>

## Syntax

```
lipid_tilt analysis-ID keyword value
```
* lipid_tilt = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'upper'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', averages over all lipid types.
    * style (str: 'angle', 'order'): Specify whether to compute the tilt angle ('angle') or the tilt angle order parameter ('order'). Default: 'angle'
    * ref_axis (str: 'x', 'y', or 'z'): Specify the reference axis that should be used to estimate the tilt. This is typically the axis along the bilayer normal. Default: 'z'

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('lipid_tilt lipid_tilt_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('lipid_tilt lipid_tilt_1 leaflet upper')
```
 
Add by list:
```python
analyzer.add_analysis(list(['lipid_tilt', 'lipid_tilt_1', dict({'leaflet':'upper'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'lipid_tilt', 'analysis_id': 'lipid_tilt_1','analysis_settings':dict({'leaflet':'upper'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('lipid_tilt_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('lipid_tilt_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
## Related analyses
* [loa](loa.html)
* [lipid_collinearity](lipid_collinearity.html)

## References

1. Anton O. Chugunov,  Pavel E. Volynsky, Nikolay A. Krylov,
Ivan A. Boldyrev, and Roman G. Efremov,  Liquid but Durable:
Molecular Dynamics Simulations Explain the Unique Properties of
Archaeal-Like Membranes, Scientific Reports, 4:7462, 2014,
doi:10.1038/srep07462
(https://www.nature.com/articles/srep07462)

