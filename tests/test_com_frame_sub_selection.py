import vorbilt.bilayer_analyzer.bilayer_analyzer as ba
def test_com_frame_sub_selection():
    analyzer = ba.BilayerAnalyzer(psf_file='../vorbilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")


    analyzer.run_analysis()
    com_msd = analyzer.get_analysis_data('msd_1')
    #reset
    analyzer.reset()
    #use the phosphorous atoms instead of full lipid center of mass
    analyzer.com_frame_name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}

    analyzer.run_analysis()

    analyzer.get_analysis_data('msd_1')
    phospho_msd = analyzer.get_analysis_data('msd_1')

    #reset
    analyzer.reset()
    #use the phosphorous and nitrogen atoms instead of full lipid center of mass -- TLCL2 has no nitrogen
    analyzer.com_frame_name_dict = {'DOPE':['P', 'N'],'POPC':['P', 'N'],'TLCL2':['P1','P3']}
    analyzer.run_analysis()
    PN_msd = analyzer.get_analysis_data('msd_1')


    print('Center of Mass MSD:')
    print(com_msd)
    print('Phosphorous MSD:')
    print(phospho_msd)
    print('Phosphorous and Nitrogen MSD:')
    print(PN_msd)


test_com_frame_sub_selection()