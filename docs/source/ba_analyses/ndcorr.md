## BilayerAnalyzer analysis: ```ndcorr``` - Correlation between bilayer surface curvature and lipid clustering.
 
## Description
 
Correlation between bilayer surfucace curvature and the clustering of lipid molecules.

This protocol is used to estimate the cross correlation between the
normal dimension deflection of lipids and the lipid types in local
blocks of the bilayer. This serves as a measure of the correlation
between the local curvature and composition of the bilayer.

This protocol is identified by the analysis key: 'ndcorr'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.NDCorrProtocol'>

## Syntax

```
ndcorr analysis-ID
```
* ndcorr = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('ndcorr ndcorr_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('ndcorr ndcorr_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['ndcorr', 'ndcorr_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'ndcorr', 'analysis_id': 'ndcorr_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('ndcorr_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('ndcorr_1')
```
 
The output is type ```<class 'dict'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> Automatically estimates the correlations with each lipid type in the bilayer selection provided to the external BilayerAnalyzer object.  </p> 
</div> 
 
## Related analyses
* None

## References

1. Koldso H, Shorthouse D, He Lie Sansom MSP (2014) "Lipid
Clustering Correlates with Membrane Curvature as Revealed by
Molecular Simulations of Complex Lipid Bilayers." PloS Comput
Biol 10(10): e1003911. doi.10.1371/journal.pcbi.1003911
