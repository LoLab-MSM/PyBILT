import vorbilt.bilayer_analyzer.bilayer_analyzer as ba
import vorbilt.com_trajectory.COMTraj as comtraj
import MDAnalysis as mda
import numpy as np
def test_run_analysis_mp():
    analyzer = ba.BilayerAnalyzer(psf_file='../vorbilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="not resname CLA and not resname TIP3 and not resname POT")

    analyzer.add_analysis('nnf nnf_a resname_1 DOPE resname_2 POPC leaflet upper n_neighbors 6')
    #analyzer.add_analysis('nnf nnf_b resname_1 POPC resname_2 POPC leaflet upper n_neighbors 10')
    #analyzer.add_analysis('bilayer_thickness bt')
    #analyzer.add_analysis('mass_dens md')
    analyzer.set_frame_range(interval=2)
    analyzer.print_analysis_protocol()

    #run serial version first
    analyzer.run_analysis()
    msd_1_dat_s = analyzer.get_analysis_data('msd_1')
    nnf_a_dat_s = analyzer.get_analysis_data('nnf_a')

    #reset and run multiprocessing version
    analyzer.reset()
    analyzer.run_analysis_mp(nprocs=2)

    msd_1_dat_mp = analyzer.get_analysis_data('msd_1')
    nnf_a_dat_mp = analyzer.get_analysis_data('nnf_a')

    print "msd serial:"
    print(msd_1_dat_s)
    print "msd_mp:"
    print(msd_1_dat_mp)
    print("nnf serial:")
    print(nnf_a_dat_s)
    print("nnf mp:")
    print(nnf_a_dat_mp)

if __name__ == '__main__':
    test_run_analysis_mp()
