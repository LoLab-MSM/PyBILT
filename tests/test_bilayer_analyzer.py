from vorbilt.bilayer_analyzer.bilayer_analyzer import BilayerAnalyzer
from vorbilt.bilayer_analyzer.analysis_protocols import valid_analysis

#define tests

#input options
def test_input_options():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    ba = BilayerAnalyzer(
        psf_file='../vorbilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    ba.add_analysis('msd msd_2 leaflet upper resname POPC')
    ba.add_plot('msd msd_p msd_1 DOPE-U msd_2 POPC-U')
    ba.run_analysis()
    ba = BilayerAnalyzer(input_file='sample_1.in')
    ba.run_analysis()

#analsyses with default options
def test_analysis_defaults():
    sel_string = "not resname CLA and not resname TIP3 and not resname POT"
    ba = BilayerAnalyzer(
        psf_file='../vorbilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../vorbilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    for key in valid_analysis:
        analysis_in = [key, key+"_t", dict()]
        ba.add_analysis(analysis_in)
    ba.print_analysis_protocol()
    ba.run_analysis()

