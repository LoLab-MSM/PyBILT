## BilayerAnalyzer analysis: ```thickness_mass_weighted_std``` - Estimates the thickness using two times the mass-weighted standard deviation of coordinates along the normal.
 
## Description
 
Estimate the thickness via the mass weighted standard deviation along the normal dimension.
This protocol is used to estimate the thickness by computing the standard
deviation of the reference atoms (e.g. phosphorous atoms) coordinates along
the bilayer normal dimension weighted by their masses:
thickness = 2xmass_weighted_std
This is same method used by MEMBPLUGIN to estimate bilayer thickness
(Blake: I had to examine the MEMBPLUGIN source code to see that this is
what it actually does.).

This protocol is identified by the analysis key: 'mass_weighted_std_thickness'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.MassWeightedStdDistProtocol'>

## Syntax

```
thickness_mass_weighted_std analysis-ID keyword value
```
* thickness_mass_weighted_std = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * selection_string (str): Provide the MDAnalysis compatible selection for the atoms to include in this analysis. Default: 'BILAYER', use all the lipids of the bilayer as recovered from the selection given to the external BilayerAnalyzer.

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```
 
Add by string - use default settings:
```python
analyzer.add_analysis('thickness_mass_weighted_std thickness_mass_weighted_std_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('thickness_mass_weighted_std thickness_mass_weighted_std_1 selection_string name P')
```
 
Add by list:
```python
analyzer.add_analysis(list(['thickness_mass_weighted_std', 'thickness_mass_weighted_std_1', dict({'selection_string':'name P'})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'thickness_mass_weighted_std', 'analysis_id': 'thickness_mass_weighted_std_1','analysis_settings':dict({'selection_string':'name P'})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('thickness_mass_weighted_std_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('thickness_mass_weighted_std_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The P-P distance can be estimated for phospholipids by assigning the phosphorous atoms as the reference atoms (e.g. with selection: name P).  </p> 
</div> 
 
## Related analyses
* None

## References

1.  R. Guixa-Gonzalez; I. Rodriguez-Espigares; J. M. Ramirez-Anguita; P.
Carrio-Gaspar; H. Martinez-Seara; T. Giorgino; J. Selent. MEMBPLUGIN: studying
membrane complexity in VMD. Bioinformatics 2014; vol. 30 (10) p. 1478-1480
doi:10.1093/bioinformatics/btu037
