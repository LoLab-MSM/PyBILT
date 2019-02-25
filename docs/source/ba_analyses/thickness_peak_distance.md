## BilayerAnalyzer analysis: ```thickness_peak_distance``` - Estimate the bilayer thickness using distance between corresponding peaks in a mass density profile.
 
## Description
 
Estimate the mass density profile for the specified selection.

This protocol is used to estimate the bilayer thickness by computing the
1-dimensional mass density profile (with two peaks) for a selection of
atoms along the bilayer normal and then estimating the distance between
between the two peaks.

This protocol is identified by the analysis key: 'thickness_peak_distance'


#### Initiated by instance of:
 
    <class 'pybilt.bilayer_analyzer.analysis_protocols.PeakToPeakProtocol'>

## Syntax

```
thickness_peak_distance analysis-ID keyword value
```
* thickness_peak_distance = analysis-Key - keyword/name for this analysis.
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
analyzer.add_analysis('thickness_peak_distance thickness_peak_distance_1') 
```
 
Add by string - adjust a setting: 
```python
analyzer.add_analysis('thickness_peak_distance thickness_peak_distance_1 n_bins 25')
```
 
Add by list:
```python
analyzer.add_analysis(list(['thickness_peak_distance', 'thickness_peak_distance_1', dict({'n_bins':25})]))
```
 
Add by dict: 
```python
analyzer.add_analysis(dict({'analysis_key': 'thickness_peak_distance', 'analysis_id': 'thickness_peak_distance_1','analysis_settings':dict({'n_bins':25})}))
```
 
To remove from analyzer: 
```python
analyzer.remove_analysis('thickness_peak_distance_1')
```
 
### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('thickness_peak_distance_1')
```
 
The output is type ```<type 'numpy.ndarray'>```
 
<div class="admonition note"> 
<p class="admonition-title">Note</p> 
<p> The P-P distance can be estimated for phospholipids by assigning the phosphorous atoms as the reference atoms for the mass density profile. The mass density profile must have two symmetric peaks or this metric will yield spurious results.  </p> 
</div> 
 
## Related analyses
* None

## References
None