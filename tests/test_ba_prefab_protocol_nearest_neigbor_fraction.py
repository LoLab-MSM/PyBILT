
from pybilt.bilayer_analyzer.prefab_analysis_protocols import nearest_neighbor_fraction


def test_analysis_module_nearest_neihbor_fraction():
    sel_string = "resname POPC DOPE TLCL2"
    nearest_neighbor_fraction(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  lipid_resnames=['POPC', 'DOPE', 'TLCL2'])

    return

if __name__ == '__main__':
    test_analysis_module_nearest_neihbor_fraction()
