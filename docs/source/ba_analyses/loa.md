## BilayerAnalyzer analysis: ```loa``` - Lateral lipid orientation angle.
 
## Description
 
Estimate the angle for lipid orientation vectors relative to the bilayer normal.

This analysis is based on the P-N vector-normal angle analysis discussed
in Ref. 1. Using two reference atoms from the specified lipid type a vector
for each lipid is computed and the orientation relative to the bilayer
normal is computed, i.e. theta, or the angle between the two vectors.

This protocol is identified by the analysis key: 'loa'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.LateralOrientationAngleProtocol'>

## Syntax

```
loa analysis-ID keyword value
```
* loa = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer leaflet to include in the estimate. Default: 'both'
    * resname (str): Specify the resname of the lipid type to include in this analysis. Default: 'first', the first lipid type as stored in the com_frame.
    * ref_atom_2 (str): The atom name of the reference atom to use as the head of the lipid orientation vector.
    * ref_atom_1 (str): The atom name of the reference atom to use as the base of the lipid orientation vector.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('loa loa_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('loa loa_1 resname POPC')
```
 
Add by list:
```python
analyzer.add_analysis(list(['loa', 'loa_1', dict({'resname':'POPC'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'loa', 'analysis_id': 'loa_1','analysis_settings':dict({'resname':'POPC'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('loa_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('loa_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The lipid orientation vector is defined as ```ref_atom_1```->```ref_atom_2```. This analysis is only performed for one lipid type, so a new intance of this analysis must be added to the BilayerAnalyzer for each lipid type that is to be analyzed.   </p> 
</div> 
 
## Related analyses
* [lop](lop.html)
* [lipid_tilt](lipid_tilt.html)

## References

1. Zheng Li, Richard M. Venable, Laura A. Rogers, Diana Murray,
and Richard W. Pastor, "Molecular Dynamics Simulations of PIP2 and
PIP3 in Lipid Bilayers: Determination of Ring Orientation, and the
Effects of Surface Roughness on a Poisson-Boltzmann Description",
Biophys J. 2009 Jul 8; 97(1): 155-163.
doi: 10.1016/j.bpj.2009.04.037
