## BilayerAnalyzer analysis: ```nnf``` - Lateral order nearest neighbor fraction.
 
## Description
 
Estimate the nearest neighbor fraction for one lipid type with another.

This analysis picks a specified number (n_neighbors) of nearest neighbors
centered on a lipid of reference lipid type and then counts the number of
lipids (M) of target lipid type and estimates the fraction,
nnf =  <M/n_neighbors> ,
where angle brackets denote averaging over lipids of specified by settings
resname_1. This metric is also referred to as 'fractional interations'.

This protocol is identified by the analysis key: 'nnf'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.NNFProtocol'>

## Syntax

```
nnf analysis-ID keyword value
```
* nnf = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * n_neighbors (int): Specifies the number of nearest neighbors to to include in computation. Default: 5
    * resname_1 (str): Specify the resname of the reference lipid type to include in this analysis. Default: 'first', the first lipid in the list pulled from the com_frame representation.
    * resname_2 (str): Specify the resname of the target lipid type to include in this analysis. Default: 'first', the first lipid in the list pulled from the com_frame representation.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('nnf nnf_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('nnf nnf_1 leaflet both')
```
 
Add by list:
```python
analyzer.add_analysis(list(['nnf', 'nnf_1', dict({'leaflet':'both'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'nnf', 'analysis_id': 'nnf_1','analysis_settings':dict({'leaflet':'both'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('nnf_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('nnf_1')
```
 
The output is type ```<class 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The anlysis in centered on lipids of the type specified by the ```resname_1``` setting and returns the fraction of lipids of the type specified by the ```resname_2``` setting within the ```n_neighbors``` nearest neighbors.  </p> 
</div> 
 
## Related analyses
* None

## References

1. A. H. de Vries, A. E. Mark and S. J. Marrink, J. Phys. Chem. B,
2004, 108, 2454-2463

2. M. Orsi and J. W. Essex, Faraday Discuss., 2013, 161, 249-272

3. Koldso H, Shorthouse D, He Lie Sansom MSP (2014) "Lipid
Clustering Correlates with Membrane Curvature as Revealed by
Molecular Simulations of Complex Lipid Bilayers." PloS Comput
Biol 10(10): e1003911. doi.10.1371/journal.pcbi.1003911
