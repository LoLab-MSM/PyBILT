## BilayerAnalyzer analysis: ```ald``` - Average lateral displacement.
 
## Description
 
Estimate the average lateral displacement of lipids.

This protocol is identified by the analysis key: 'ald'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.ALDProtocol'>

## Syntax

```
ald analysis-ID keyword value
```
* ald = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', includes all lipid types.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('ald ald_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('ald ald_1 leaflet both')
```
 
Add by list:
```python
analyzer.add_analysis(list(['ald', 'ald_1', dict({'leaflet':'both'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'ald', 'analysis_id': 'ald_1','analysis_settings':dict({'leaflet':'both'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('ald_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('ald_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
## Related analyses
* [msd](msd.html)

## References

1. Kenichiro Koshiyama, Tetsuya Kodama, Takeru Yano, Shigeo
Fujikawa, "Molecular dynamics simulation of structural changes
of lipid bilayers induced by shock waves: Effects of incident
angles", Biochimica et Biophysica Acta (BBA) - Biomembranes,
Volume 1778, Issue 6, June 2008, Pages 1423-1428
