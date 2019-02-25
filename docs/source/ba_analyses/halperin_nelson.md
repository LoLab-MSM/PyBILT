## BilayerAnalyzer analysis: ```halperin_nelson``` - Halperin and Nelson's rotational invariant.
 
## Description
 
Estimate the mean hexagonal packing orientation parameter.

This analysis implements Halperin and
Nelson's rotational invariant to estimate the hexagonal packing orientation
parameter. The value for lipid l is given by
phi_l = | (1/6) * sum_{j element nn(l)} exp(6i*theta_{lj}) |^2 ,
where i is complex and nn(l) are the 6 nearest neighbors of lipid
l; theta_{lj} is the angle between the vector formed by beads
representing lipid l and j and an arbitrary axis. The value is unity
for perfect hexagonal packing, and it is zero to the extent that
hexagonal packing is entirely absent. This protocol uses the 'com_frame'
representation of the bilayer.

This protocol is identified by the analysis key: 'halperin_nelson'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.HalperinNelsonProtocol'>

## Syntax

```
halperin_nelson analysis-ID keyword value
```
* halperin_nelson = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'upper'

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('halperin_nelson halperin_nelson_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('halperin_nelson halperin_nelson_1 leaflet upper')
```
 
Add by list:
```python
analyzer.add_analysis(list(['halperin_nelson', 'halperin_nelson_1', dict({'leaflet':'upper'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'halperin_nelson', 'analysis_id': 'halperin_nelson_1','analysis_settings':dict({'leaflet':'upper'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('halperin_nelson_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('halperin_nelson_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
## Related analyses
* None

## References

1. Shachi Katira, Kranthi K. Mandadapu, Suriyanarayanan
Vaikuntanathan, Berend Smit, and David Chandler, The
order-disorder transition in model lipid bilayers is a
first-order hexatic to liquid phase transition, arXiv preprint
[cond-mat.soft] 2015, arXiv:1506.04310.
https://arxiv.org/pdf/1506.04310.pdf


