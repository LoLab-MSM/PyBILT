## BilayerAnalyzer analysis: ```vcm``` - Volume compressibility modulus.
 
## Description
 
Estimate the isothermal volume compressibility modulus.

This protocol is used to estimate the volume compressibility modulus,
K_V = (<V>kT) / var(V),
where V is the volume.

This protocol is identified by the analysis key: 'vcm'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.VolumeCompressibilityModulusProtocol'>

## Syntax

```
vcm analysis-ID keyword value
```
* vcm = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * temperature (float): The absolute temperature that the simulation was run at (i.e. in Kelvin). Default: 298.15 K

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('vcm vcm_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('vcm vcm_1 temperature 298.15')
```
 
Add by list:
```python
analyzer.add_analysis(list(['vcm', 'vcm_1', dict({'temperature':298.15})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'vcm', 'analysis_id': 'vcm_1','analysis_settings':dict({'temperature':298.15})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('vcm_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('vcm_1')
```
 
The output is type ```<class 'numpy.ndarray'>```
 
## Related analyses
* [acm](acm.html)
* [ac](ac.html)

## References

1. Christofer Hofsab, Erik Lindahl, and Olle Edholm, "Molecular
Dynamics Simulations of Phospholipid Bilayers with Cholesterol",
Biophys J. 2003 Apr; 84(4): 2192-2206.
doi:  10.1016/S0006-3495(03)75025-5
