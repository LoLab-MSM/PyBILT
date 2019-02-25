## BilayerAnalyzer analysis: ```apl_grid``` - Area per lipid using 2D lipid grids.
 
## Description
 
Estimate the indvidual area per lipid for each lipid type using a gridding procedure.

This protocol is identified by the analysis key: 'apl_grid'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.APLGridProtocol'>

## Syntax

```
apl_grid analysis-ID
```
* apl_grid = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('apl_grid apl_grid_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('apl_grid apl_grid_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['apl_grid', 'apl_grid_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'apl_grid', 'analysis_id': 'apl_grid_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('apl_grid_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('apl_grid_1')
```
 
The output is type ```<type 'dict'>```
 
## Related analyses
* [apl_box](apl_box.html)

## References

1. Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry
2. Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858
