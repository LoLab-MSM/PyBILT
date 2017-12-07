
from pybilt.bilayer_analyzer.prefab_analysis_protocols import lipid_length


def test_analysis_module_lipid_length():
    sel_string = "resname POPC DOPE TLCL2"
    ref_atoms = {'DOPE':{'start':['C218','C318'], 'end':'P'},
                 'POPC':{'start':['C218', 'C316'],'end':'P'},
                 'TLCL2':{'start':['CA18','CB18','CC18', 'CD18'],
                          'end':['P1', 'P3']}}
    lipid_length(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                       trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                       selection_string=sel_string,
                       lipid_resnames=['POPC', 'DOPE', 'TLCL2'],
                       ref_atoms=ref_atoms,
                       frame_interval=2)



    return

if __name__ == '__main__':
    test_analysis_module_lipid_length()
