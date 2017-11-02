from pybilt.bilayer_analyzer import BilayerAnalyzer
import numpy as np

def test_lipid_grid_pdb_writer():
    sel_string = "resname POPC DOPE TLCL2"
    name_dict = {'DOPE':['P'],'POPC':['P'],'TLCL2':['P1','P3']}
    analyzer = BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")
    analyzer.remove_analysis('msd_1')
    analyzer.rep_settings['com_frame']['name_dict'] = name_dict
    analyzer._analysis_protocol.use_objects['lipid_grid'] = True
    # adjust the number of bins for gridding
    nbins = 75
    analyzer.rep_settings['lipid_grid']['n_xbins'] = nbins
    analyzer.rep_settings['lipid_grid']['n_ybins'] = nbins
    analyzer.set_frame_range(0,5,1)
    i =0
    zgrid = None
    for _frame in analyzer:
        if i == 0:
            zgrid = analyzer.reps['lipid_grid'].leaf_grid['upper'].lipid_grid_z
        else:
            zgrid += analyzer.reps['lipid_grid'].leaf_grid['upper'].lipid_grid_z
        i += 1
    zgrid /= i
    analyzer.reps['lipid_grid'].write_pdb('lipid_grid_t0.pdb', leaflet='both',
                                          beta_grid_upper=zgrid)

    analyzer.reps['lipid_grid'].write_pdb('lipid_grid_t1.pdb', leaflet='upper',
                                          z_grid_upper=zgrid)
    analyzer.reps['lipid_grid'].write_pdb('lipid_grid_t2.pdb', leaflet='both',
                                          use_gaussian_filter=True,
                                          filter_sigma=5.0)                                          

    return

if __name__ == '__main__':
    test_lipid_grid_pdb_writer()
