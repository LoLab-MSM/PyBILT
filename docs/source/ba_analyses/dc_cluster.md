## BilayerAnalyzer analysis: ```dc_cluster``` - Hiearchical clustering of lipids based on distance.
 
## Description
 
Compute lipid clusters using a hiearchical distance based method.

This analysis uses a type of hiearchical clustering where points (lipid
centers of mass) are are added to a cluster if they are within a specified
distance of any other point within the cluster.

This protocol is identified by the analysis key: 'dc_cluster'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.DCClusterProtocol'>

## Syntax

```
dc_cluster analysis-ID keyword value
```
* dc_cluster = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'all', includes all lipid types.
    * cutoff (float): The cutoff distance to use for the clustering. Default: 12.0

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('dc_cluster dc_cluster_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('dc_cluster dc_cluster_1 resname POPC')
```
 
Add by list:
```python
analyzer.add_analysis(list(['dc_cluster', 'dc_cluster_1', dict({'resname':'POPC'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'dc_cluster', 'analysis_id': 'dc_cluster_1','analysis_settings':dict({'resname':'POPC'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('dc_cluster_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('dc_cluster_1')
```
 
The output is type ```<type 'dict'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> Only finds the self clusters for a single lipid type as specified by the ```resname``` setting.  </p> 
</div> 
 
## Related analyses
* None

## References
None