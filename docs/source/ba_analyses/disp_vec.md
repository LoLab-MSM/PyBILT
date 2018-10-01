## BilayerAnalyzer analysis: ```disp_vec``` - Displacement vectors.
 
## Description
 
Comute displacement vectors for each lipid in the specified
leaflet(s) of bilayer.

This protocol is identified by the analysis key: 'disp_vec'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.DispVecProtocol'>

## Syntax

```
disp_vec analysis-ID keyword value
```
* disp_vec = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', includes all lipid types.
    * wrapped (bool): Specify whether to use the wrapped ('True') or un-wrapped ('False') coordintes for the base of the vectors. Default: False
    * interval (int): Sets the frame interval over which to compute the displacement vectors.
    * scale (bool): Specify whether to scale the coordinates by the box dimensions of the reference frame. Default: False
    * scale_to_max (bool): Specify whether to scale the coordinates by the box dimensions of the maximum box size in the anlysis. Default: False.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('disp_vec disp_vec_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('disp_vec disp_vec_1 leaflet both')
```
 
Add by list:
```python
analyzer.add_analysis(list(['disp_vec', 'disp_vec_1', dict({'leaflet':'both'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'disp_vec', 'analysis_id': 'disp_vec_1','analysis_settings':dict({'leaflet':'both'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('disp_vec_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('disp_vec_1')
```
 
The output is type ```<class 'list'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The real frame ```interval``` setting needs to be a multiple of the frame looping interval set within the BilayerAnalyzer instance (i.e., analyzer.frame_range[2]) for this analysis to work properly.  </p> 
</div> 
 
## Related analyses
* [disp_vec_corr](disp_vec_corr.html)
* [disp_vec_nncorr](disp_vec_nncorr.html)
* [disp_vec_corr_avg](disp_vec_corr_avg.html)
* [spatial_velocity_corr](spatial_velocity_corr.html)

## References

1. Emma Falck, Tomasz Rog, Mikko Karttunen, and Ilpo Vattulainen,
Lateral Diffusion in Lipid Membranes through Collective Flows,
Journal of the American Chemical Society, 2008 130 (1), 44-45
DOI: 10.1021/ja7103558
