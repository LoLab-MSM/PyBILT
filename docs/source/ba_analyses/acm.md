## BilayerAnalyzer analysis: ```acm``` - Area compressibility modulus.

## Description

Estimate the isothermal area compressibility modulus.

This protocol is used to estimate the area compressibility modulus,

    K_A = (<A>kT) / var(A),

where A is the area in the lateral dimension of the bilayer.

This protocol is identified by the analysis key: 'acm'


#### Initiated by instance of:

    <class 'pybilt.bilayer_analyzer.analysis_protocols.AreaCompressibilityModulusProtocol'>

## Syntax

```
acm analysis-ID keyword value
```
* acm = analysis-Key - keyword/name for this analysis.
* analysis-ID = The unique name/ID being assigned to this analysis.
* keyword value = settings keyword value pairs
    * temperature (float): The absolute temperature that the simulation was run at (i.e. in Kelvin). Default: 298.15 K

### Examples
Construct analyzer:
```python
analyzer = BilayerAnalyzer(structure='name_of_structure_file',
                           trajectory='name_of_traj_file',
                           selection='resname POPC DOPC')
```

Add by string - use default settings:
```python
analyzer.add_analysis('acm acm_1')
```

Add by string - adjust a setting:
```python
analyzer.add_analysis('acm acm_1 temperature 298.15')
```

Add by list:
```python
analyzer.add_analysis(list(['acm', 'acm_1', dict({'temperature':298.15})]))
```

Add by dict:
```python
analyzer.add_analysis(dict({'analysis_key': 'acm', 'analysis_id': 'acm_1','analysis_settings':dict({'temperature':298.15})}))
```

To remove from analyzer:
```python
analyzer.remove_analysis('acm_1')
```

### Output Info:
Retrieve output after running analyses:
```python
output = analyzer.get_analysis_data('acm_1')
```

The output is type ```<class 'numpy.ndarray'>```

## Related analyses
* [ac](ac.html)
* [vcm](vcm.html)

## References

1. Christofer Hofsab, Erik Lindahl, and Olle Edholm, "Molecular
Dynamics Simulations of Phospholipid Bilayers with Cholesterol",
Biophys J. 2003 Apr; 84(4): 2192-2206.
doi:  10.1016/S0006-3495(03)75025-5
2. L. Janosi and A. A. Gorfe, J. Chem. Theory Comput. 2010, 6,
3267-3273
3. D. Aguayo, F. D. Gonzalez-Nilo, and C. Chipot, J. Chem. Theory
Comput. 2012, 8, 1765-1773
