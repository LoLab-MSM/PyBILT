
from pybilt.bilayer_analyzer.prefab_analysis_protocols import com_lateral_rdf


def test_prefab_protocol_com_lateral_rdf():
    sel_string = "resname POPC DOPE TLCL2"
    print("Run...")
    com_lateral_rdf(structure_file='../pybilt/sample_bilayer/sample_bilayer.psf',
                  trajectory_file='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  bilayer_selection_string=sel_string,
                  resnames=['POPC', 'DOPE', 'TLCL2'], frame_interval=2)

    return

if __name__ == '__main__':
    test_prefab_protocol_com_lateral_rdf()
