## BilayerAnalyzer analysis: ```thickness_grid``` - Bilayer thickness using lipid_grid.
 
## Description
 
Estimate the bilayer thickness using a gridding procedure.

This analysis uses a lipid grid representation (lipid positions are mapped
to 2-D grids, one grid per leaflet) of the bilayer to estimate the bilayer
thickness. This is done by measuring the difference in 'z' position between
corresponding grid points  of opposing bilayer leaflets.

This protocol is identified by the analysis key: 'thickness_grid'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.BTGridProtocol'>

## Syntax

```
thickness_grid analysis-ID
```
* thickness_grid = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('thickness_grid thickness_grid_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('thickness_grid thickness_grid_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['thickness_grid', 'thickness_grid_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'thickness_grid', 'analysis_id': 'thickness_grid_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('thickness_grid_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('thickness_grid_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The P-P distance can be estimated for phospholipid by assigning positions of the phosphorous atoms to the lipid grid.  </p> 
</div> 
 
## Related analyses
* [lipid_length](lipid_length.html)
* [thickness_leaflet_distance](thickness_leaflet_distance.html)
* [thickness_peak_distance](thickness_peak_distance.html)
* [thickness_mass_weighted_std](thickness_mass_weighted_std.html)

## References

1. Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry
2. Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858
