from pybilt.analysis_modules import area_per_lipid


def test_analysis_module_area_per_lipid():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    area_per_lipid(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string)

    return

if __name__ == '__main__':
    test_analysis_module_area_per_lipid()