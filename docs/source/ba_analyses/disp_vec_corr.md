## BilayerAnalyzer analysis: ```disp_vec_corr``` - Displacement vector correlation matrix.
 
## Description
 
Comute the pair-wise cross correlation matrix for the displacement vectors for each lipid in the specified leaflet(s) of bilayer.

This analysis computes the displacement vectors as in 'disp_vec' analysis
(DispVecProtocol), but then continues to compute the pair-wise cross
correlation matrix between each vector. i.e. the cos(theta) for the angle
theta between the vectors.

This protocol is identified by the analysis key: 'disp_vec_corr'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.DispVecCorrelationProtocol'>

## Syntax

```
disp_vec_corr analysis-ID keyword value
```
* disp_vec_corr = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', includes all lipid types.
    * wrapped (bool): Specify whether to use the wrapped ('True') or un-wrapped ('False') coordintes for the base of the vectors. Default: False
    * interval (int): Sets the frame interval over which to compute the displacement vectors. f

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('disp_vec_corr disp_vec_corr_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('disp_vec_corr disp_vec_corr_1 leaflet both')
```
 
Add by list:
```python
analyzer.add_analysis(list(['disp_vec_corr', 'disp_vec_corr_1', dict({'leaflet':'both'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'disp_vec_corr', 'analysis_id': 'disp_vec_corr_1','analysis_settings':dict({'leaflet':'both'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('disp_vec_corr_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('disp_vec_corr_1')
```
 
The output is type ```<class 'list'>```
 
## Related analyses
* [disp_vec](disp_vec.html)
* [disp_vec_nncorr](disp_vec_nncorr.html)
* [disp_vec_corr_avg](disp_vec_corr_avg.html)
* [spatial_velocity_corr](spatial_velocity_corr.html)

## References
None