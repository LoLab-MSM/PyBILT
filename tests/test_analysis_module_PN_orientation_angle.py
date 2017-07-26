
from pybilt.bilayer_analyzer.analysis_modules import PN_orientational_angle


def test_analysis_module_PN_orientation_anlge():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    PN_orientational_angle(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  lipid_resnames=['POPC', 'DOPE'])

    return

if __name__ == '__main__':
    test_analysis_module_PN_orientation_anlge()
