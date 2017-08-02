from pybilt.analysis_modules import bilayer_thickness


def test_analysis_module_bilayer_thickness():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    bilayer_thickness(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  selection_string=sel_string,
                  name_dict=name_dict,
                  frame_start=2, frame_end=-2, frame_interval=2,
                  n_xbins=30, n_ybins=30)

    return

if __name__ == '__main__':
    test_analysis_module_bilayer_thickness()