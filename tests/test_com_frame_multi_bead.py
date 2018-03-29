import pybilt.bilayer_analyzer.bilayer_analyzer as ba
def test_com_frame_multi_bead():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")


    analyzer.remove_analysis('msd_1')
    analyzer.rep_settings['com_frame']['name_dict'] = {'DOPE':['C2'],'POPC':['C2'],'TLCL2':['C13','C32']}
    analyzer.rep_settings['com_frame']['multi_bead'] = True
    analyzer.add_analysis('halperin_nelson halperin_nelson_upper leaflet upper')
    analyzer.run_analysis()
    com_halp = analyzer.get_analysis_data('halperin_nelson_upper')
    nbeads = len(analyzer.reps['com_frame'])
    print("There are {} COM beads.".format(nbeads))
    print(com_halp)

if __name__ == '__main__':
    test_com_frame_multi_bead()
