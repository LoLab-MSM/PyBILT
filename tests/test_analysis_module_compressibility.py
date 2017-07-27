
from pybilt.analysis_modules import compressibility


def test_analysis_module_compressibility():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    compressibility(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  temperature=303.15)

    return

if __name__ == '__main__':
    test_analysis_module_compressibility()