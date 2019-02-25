## BilayerAnalyzer analysis: ```mass_dens``` - Mass density profile.
 
## Description
 
Estimate the mass density profile for the specified selection.

This protocol is used to estimate the 1-dimensional mass density profile
for a selection of atoms along the bilayer normal. The profile is
automatically centered on the bilayer's center of mass along the bilayer
normal.

This protocol is identified by the analysis key: 'mass_dens'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.MassDensProtocol'>

## Syntax

```
mass_dens analysis-ID keyword value
```
* mass_dens = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs 
    * n_bins (int): Set the number of bins to divide the normal dimensions into for binning. Defalt: 25
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
analyzer.add_analysis('mass_dens mass_dens_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('mass_dens mass_dens_1 n_bins 25')
```
 
Add by list:
```python
analyzer.add_analysis(list(['mass_dens', 'mass_dens_1', dict({'n_bins':25})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'mass_dens', 'analysis_id': 'mass_dens_1','analysis_settings':dict({'n_bins':25})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('mass_dens_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('mass_dens_1')
```
 
The output is type ```<type 'tuple'>```
 
## Related analyses
* None

## References
None