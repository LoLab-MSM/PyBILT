from pybilt.bilayer_analyzer.prefab_analysis_protocols import normal_displacement_lipid_type_correlation


def test_analysis_module_normal_displacement_lipid_type_correlation():
    sel_string = "resname POPC DOPE TLCL2"
    com_sub_selection_dict =  {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    normal_displacement_lipid_type_correlation(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  com_sub_selection_dict=com_sub_selection_dict)

    return

if __name__ == '__main__':
    test_analysis_module_normal_displacement_lipid_type_correlation()
