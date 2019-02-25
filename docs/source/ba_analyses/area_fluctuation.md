## BilayerAnalyzer analysis: ```area_fluctuation``` - Bilayer lateral box area fluctuation.
 
## Description
 
Estimate the area fluctuation in the box along the bilayer laterals.

This protocol is identified by the analysis key: 'area_fluctuation'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.AreaFluctuationProtocol'>

## Syntax

```
area_fluctuation analysis-ID
```
* area_fluctuation = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('area_fluctuation area_fluctuation_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('area_fluctuation area_fluctuation_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['area_fluctuation', 'area_fluctuation_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'area_fluctuation', 'analysis_id': 'area_fluctuation_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('area_fluctuation_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('area_fluctuation_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
## Related analyses
* [acm](acm.html)
* [ac](ac.html)

## References
None