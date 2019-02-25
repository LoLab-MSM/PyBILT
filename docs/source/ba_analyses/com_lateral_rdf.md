## BilayerAnalyzer analysis: ```com_lateral_rdf``` - Lipid-lipid RDF in the bilayer lateral plane.
 
## Description
 
Estimate the 2-d radial pair distribution function in the bilayer lateral plane using the lipid centers of mass.

This analysis protocol uses the 'com_frame' representation.

This protocol is identified by the analysis key: 'com_lateral_rdf'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.COMLateralRDFProtocol'>

## Syntax

```
com_lateral_rdf analysis-ID keyword value
```
* com_lateral_rdf = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * n_bins (int): Specifies the number of bins to use when estimating the RDF. Default: 25
    * range_outer (float): Specify the outer distance cutoff for the RDF. Default: 25.0
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * range_inner (float): Specify the inner distance cutoff for the RDF. Default: 0.0
    * resname_2 (str): Specify the resname of the target lipid type to include in this analysis. Default: 'first', the first lipid in the list pulled from the com_frame representation.
    * resname_1 (str): Specify the resname of the reference lipid type to include in this analysis. Default: 'first', the first lipid in the list pulled from the com_frame representation.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('com_lateral_rdf com_lateral_rdf_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('com_lateral_rdf com_lateral_rdf_1 n_bins 25')
```
 
Add by list:
```python
analyzer.add_analysis(list(['com_lateral_rdf', 'com_lateral_rdf_1', dict({'n_bins':25})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'com_lateral_rdf', 'analysis_id': 'com_lateral_rdf_1','analysis_settings':dict({'n_bins':25})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('com_lateral_rdf_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('com_lateral_rdf_1')
```
 
The output is type ```<type 'tuple'>```
 
## Related analyses
* [nnf](nnf.html)

## References

1. Microsecond Molecular Dynamics Simulations of Lipid Mixing
Chunkit Hong, D. Peter Tieleman, and Yi Wang Langmuir 2014 30
(40), 11993-12001 DOI: 10.1021/la502363b
http://pubs.acs.org/doi/abs/10.1021/la502363b
