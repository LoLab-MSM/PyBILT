
from pybilt.bilayer_analyzer.prefab_analysis_protocols import spatial_velocity_correlation_functions


def test_prefab_protocol_spatial_velocity_correlation_functions():
    sel_string = "resname POPC DOPE TLCL2"
    print("Run...")
    spatial_velocity_correlation_functions(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  bilayer_selection_string=sel_string,
                  resnames=['POPC', 'DOPE', 'TLCL2'], frame_interval=3)

    return

if __name__ == '__main__':
    test_prefab_protocol_spatial_velocity_correlation_functions()
