## BilayerAnalyzer analysis: ```msd_dist``` - Mean squared displacement as a function of distance.
 
## Description
 
Compute the MSD over a set time interval as function of nearest distance.

This protocol is identified by the analysis key: 'disp_vec'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.MSDDistProtocol'>

## Syntax

```
msd_dist analysis-ID keyword value
```
* msd_dist = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * range_inner (float): Specify the inner distance cutoff for the RDF. Default: 0.0
    * n_bins (int): Specifies the number of bins to use when estimating the RDF. Default: 25
    * interval (int): Sets the frame interval over which to compute the displacement vectors.
    * range_outer (float): Specify the outer distance cutoff for the RDF. Default: 25.0
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * wrapped (bool): Specify whether to use the wrapped ('True') or un-wrapped ('False') coordintes for the base of the vectors. Default: False
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', includes all lipid types.
    * resname_2 (str): Specify the resname of the target lipid type to include in this analysis. Special names are 'first' and 'all', which use the first and all lipid types respectively. Default: 'first', the first lipid in the list pulled from the com_frame representation.
    * resname_1 (str): Specify the resname of the reference lipid type to include in this analysis. Special names are 'first' and 'all', which use the first and all lipid types respectively. Default: 'first', the first lipid in the list pulled from the com_frame representation.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('msd_dist msd_dist_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('msd_dist msd_dist_1 n_bins 25')
```
 
Add by list:
```python
analyzer.add_analysis(list(['msd_dist', 'msd_dist_1', dict({'n_bins':25})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'msd_dist', 'analysis_id': 'msd_dist_1','analysis_settings':dict({'n_bins':25})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('msd_dist_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('msd_dist_1')
```
 
The output is type ```<type 'tuple'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The anlysis is centered on lipids of the type specified by the ```resname_1```. The real frame ```interval``` setting needs to be a multiple of the frame looping interval set within the BilayerAnalyzer instance (i.e., analyzer.frame_range[2]) for this analysis to work properly.  </p> 
</div> 
 
## Related analyses
* [disp_vec](disp_vec.html)
* [msd](msd.html)

## References
None