## BilayerAnalyzer analysis: ```ac``` - Isothermal area compressibility.
 
## Description
 
Estimate the isothermal area compressibility.

This protocol is used to estimate the area compressibility modulus,
K_A^-1 = [(<A>kT) / var(A)]^-1 ,
where A is the area in the lateal dimension of the bilayer.

This protocol is identified by the analysis key: 'ac'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.AreaCompressibilityProtocol'>

## Syntax

```
ac analysis-ID keyword value
```
* ac = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('ac ac_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('ac ac_1 temperature 298.15')
```
 
Add by list:
```python
analyzer.add_analysis(list(['ac', 'ac_1', dict({'temperature':298.15})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'ac', 'analysis_id': 'ac_1','analysis_settings':dict({'temperature':298.15})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('ac_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('ac_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> Area compressibility is the inverse of area compressibility modulus (see Ref.1).  </p> 
</div> 
 
## Related analyses
* [acm](acm.html)
* [vcm](vcm.html)

## References

1. Yoshimichi Andoha, Susumu Okazakia, Ryuichi Ueokab, "Molecular
dynamics study of lipid bilayers modeling the plasma membranes of
normal murine thymocytes and leukemic GRSL cells", Biochimica et
Biophysica Acta (BBA) - Biomembranes, Volume 1828, Issue 4, April
2013, Pages 1259-1270. https://doi.org/10.1016/j.bbamem.2013.01.005
