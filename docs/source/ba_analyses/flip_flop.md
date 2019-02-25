## BilayerAnalyzer analysis: ```flip_flop``` - Count lipid flip flops.
 
## Description
 
Count any lipid flips flops between the leaflets.

This protocol is identified by the analysis key: 'flip_flop'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.FlipFlopProtocol'>

## Syntax

```
flip_flop analysis-ID
```
* flip_flop = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('flip_flop flip_flop_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('flip_flop flip_flop_1 none None')
```
 
Add by list:
```python
analyzer.add_analysis(list(['flip_flop', 'flip_flop_1', dict({'none':None})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'flip_flop', 'analysis_id': 'flip_flop_1','analysis_settings':dict({'none':None})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('flip_flop_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('flip_flop_1')
```
 
The output is type ```<type 'dict'>```
 
## Related analyses
* None

## References

1. Andrey A. Gurtovenko, and Ilpo Vattulainen, Molecular Mechanism
for Lipid Flip-Flops, The Journal of Physical Chemistry B, 2007
111 (48), 13554-13559, DOI: 10.1021/jp077094k
http://pubs.acs.org/doi/abs/10.1021/jp077094k?journalCode=jpcbfk

2. Nicolas Sapay, W. F. Drew Bennett, and D. Peter Tieleman,
Molecular Simulations of Lipid Flip-Flop in the Presence of
Model Transmembrane Helices, Biochemistry, 2010 49 (35),
7665-7673, DOI: 10.1021/bi100878q
http://pubs.acs.org/doi/abs/10.1021/bi100878q
