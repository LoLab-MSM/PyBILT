## BilayerAnalyzer analysis: ```lipid_length``` - Lipid length using the lipid vectors.
 
## Description
 
Estimate the lipids length using the defined lipid vector.

This analysis is used to compute the mean lipid length using
the vector represetation of the specified lipids.

This protocol is identified by the analysis key: 'lipid_length'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.LipidLengthProtocol'>

## Syntax

```
lipid_length analysis-ID keyword value
```
* lipid_length = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('lipid_length lipid_length_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('lipid_length lipid_length_1 leaflet both')
```
 
Add by list:
```python
analyzer.add_analysis(list(['lipid_length', 'lipid_length_1', dict({'leaflet':'both'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'lipid_length', 'analysis_id': 'lipid_length_1','analysis_settings':dict({'leaflet':'both'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('lipid_length_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('lipid_length_1')
```
 
The output is type ```<class 'numpy.ndarray'>```
 
## Related analyses
* None

## References

1. Anton O. Chugunov,  Pavel E. Volynsky, Nikolay A. Krylov,
Ivan A. Boldyrev, and Roman G. Efremov,  Liquid but Durable:
Molecular Dynamics Simulations Explain the Unique Properties of
Archaeal-Like Membranes, Scientific Reports, 4:7462, 2014,
doi:10.1038/srep07462
(https://www.nature.com/articles/srep07462)

