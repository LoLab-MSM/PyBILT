import pybilt.bilayer_analyzer.bilayer_analyzer as ba
import pybilt.plot_generation.plot_generation_functions as pgf
def test_ba_ap_dv_stm():
    analyzer = ba.BilayerAnalyzer(structure='../pybilt/sample_bilayer/sample_bilayer.psf',
                                  trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                                  selection="resname POPC DOPE TLCL2")


    analyzer.remove_analysis('msd_1')
    analyzer.add_analysis('disp_vec disp_vec interval 5 wrapped True leaflet upper scale_to_max True')
    analyzer.add_analysis('disp_vec disp_vec_b interval 1 wrapped True leaflet upper scale_to_max True')


    analyzer.print_analysis_protocol()

    analyzer.run_analysis()

    output = analyzer.get_analysis_data('disp_vec')
    dv_res = output[0]
    pgf.plot_step_vectors(dv_res, save=False, show=True, scaled=True, wrapped=True)
    dv_res = output[len(output)-1]
    pgf.plot_step_vectors(dv_res, save=False, show=True, scaled=True, wrapped=True)
    output = analyzer.get_analysis_data('disp_vec_b')
    pgf.plot_step_vectors_stroboscopic(output, index=0, scaled=True, wrapped=True, save=False, show=True)
if __name__ == '__main__':
    test_ba_ap_dv_stm()
