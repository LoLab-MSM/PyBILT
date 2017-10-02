
from from pybilt.bilayer_analyzer.prefab_analysis_protocols import PN_orientational_angle


def test_analysis_module_PN_orientation_anlge():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    print("Run with default PN_name...")
    PN_orientational_angle(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  lipid_resnames=['POPC', 'DOPE'])
    print("Run with explicit PN_name...")
    PN_orientational_angle(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  lipid_resnames=['POPC', 'DOPE'],
                  PN_names={'POPC':('P', 'N'), "DOPE":('P', 'N')})

    return

if __name__ == '__main__':
    test_analysis_module_PN_orientation_anlge()
