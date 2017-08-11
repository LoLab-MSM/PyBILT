from pybilt.analysis_modules import position_density_maps_2d


def test_analysis_module_position_density_maps_2d():

    position_density_maps_2d(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  bilayer_selection_string='resname DOPE POPC TLCL2',
                  resnames=['TLCL2', 'POPC', 'DOPE'], frame_start=0, frame_end=-1, frame_interval=1, dump_path='./')

    return

if __name__ == '__main__':
    test_analysis_module_position_density_maps_2d()