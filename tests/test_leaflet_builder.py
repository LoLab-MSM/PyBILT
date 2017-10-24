from pybilt.bilayer_analyzer import BilayerAnalyzer

def test_leaflet_builder():
    sel_string = "resname POPC DOPE TLCL2"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    analyzer = BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")

    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    # Assign the leaflets using the 'avg_norm' method. (This is actually the
    # default).
    analyzer.adjust_rep_setting('leaflets', 'assign_method', 'avg_norm')
    analyzer.run_analysis()
    analyzer.reset()
    # Now redo it using the 'orientation' mehtod to assign leaflets.""
    analyzer.adjust_rep_setting('leaflets', 'assign_method', 'orientation')
    analyzer.adjust_rep_setting('leaflets', 'orientation_atoms', {'DOPE': ['C218','P'],
                                'POPC': ['C218', 'P'], 'TLCL2': ['CA18', 'P1']})
    analyzer.run_analysis()


    return

if __name__ == '__main__':
    test_leaflet_builder()
