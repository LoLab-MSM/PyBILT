## BilayerAnalyzer analysis: ```thickness_leaflet_distance``` - Distance between the leaflets' centers of mass along the normal.
 
## Description
 
Estimate the bilayer thickness using the distance between the COMs of each leaflet.

This analysis uses the COM postions from the com_frame representation. At each
frame the centers-of-mass of each leaflet is computed and the distance along
the bilayer normal is computed. This is another to estimate the bilayer thickness.

This protocol is identified by the analysis key: 'thickness_leaflet_distance'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.LeaftoLeafDistProtocol'>

## Syntax

```
thickness_leaflet_distance analysis-ID
```
* thickness_leaflet_distance = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('thickness_leaflet_distance thickness_leaflet_distance_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('thickness_leaflet_distance thickness_leaflet_distance_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['thickness_leaflet_distance', 'thickness_leaflet_distance_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'thickness_leaflet_distance', 'analysis_id': 'thickness_leaflet_distance_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('thickness_leaflet_distance_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('thickness_leaflet_distance_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The P-P distance can be estimated for phospholipids by assigning the phosphorous atoms as the reference atoms to the com_frame.  </p> 
</div> 
 
## Related analyses
* [lipid_length](lipid_length.html)
* [thickness_grid](thickness_grid.html)
* [thickness_peak_distance](thickness_peak_distance.html)
* [thickness_mass_weighted_std](thickness_mass_weighted_std.html)

## References
None