
from pybilt.bilayer_analyzer.prefab_analysis_protocols import compressibility


def test_analysis_module_compressibility():
    sel_string = "resname POPC DOPE TLCL2"
    compressibility(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  temperature=303.15)

    return

if __name__ == '__main__':
    test_analysis_module_compressibility()
