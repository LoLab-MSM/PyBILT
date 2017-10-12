from pybilt.bilayer_analyzer import BilayerAnalyzer
from pybilt.bilayer_analyzer.analysis_protocols import valid_analysis

#define tests

#input options
def test_input_options():
    print("testing various input options...")
    #initialize analyzer with keyword options--and default analysis
    print("Initialize using keyword options:")
    sel_string = "segid MEMB"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    #add analysis by string input
    ba.add_analysis('msd msd_2 leaflet upper resname POPC')
    #add plot by string -- currenlty only available
    ba.add_plot('msd msd_p msd_1 DOPE-U msd_2 POPC-U')
    #add analysis by list/tuple input
    ba.add_analysis(['msd', 'msd_3', {'resname': 'POPC', 'leaflet': 'lower'}])
    #add analysis by dictionary input
    analysis_in_dict = dict()
    analysis_in_dict['analysis_key'] = 'msd'
    analysis_in_dict['analysis_id'] = 'msd_4'
    analysis_in_dict['analysis_settings'] = {'resname': 'DOPE', 'leaflet': 'lower'}
    ba.add_analysis(analysis_in_dict)
    # run the analyses
    ba.run_analysis()
    #initialize analyzer using input script
    print("Initialize with an input script:")
    ba = BilayerAnalyzer(input_file='sample_1.in')
    #run analysis
    ba.run_analysis()

    #now initialize with an input dictionary
    # define the input dictionary
    input_dict = {'structure' : '../pybilt/sample_bilayer/sample_bilayer.psf',
                 'trajectory' : '../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
                  'selection' : 'segid MEMB'
                 }

    #now initialize the analyzer
    print("Initializing using input dictionary:")
    ba = BilayerAnalyzer(input_dict=input_dict)

    ba.add_analysis('msd msd_b resname DOPE leaflet upper')
    ba.run_analysis()

    return

#analsyses with default options
def test_analysis_defaults():
    print("testing the default settings of all analysis protocols...")
    sel_string = "segid MEMB"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    for key in valid_analysis:
        analysis_in = [key, key+"_t", dict()]
        ba.add_analysis(analysis_in)
    ba.settings['print_interval'] = 1
    ba.print_analysis_protocol()
    ba.run_analysis()
    print("outputs:")
    for key in valid_analysis:
        a_id = key+"_t"
        print(a_id+":")
        print(ba.get_analysis_data(a_id))


def test_analysis_iterator():
    print("testing the analysis iterator...")
    # initialize analyzer with keyword options--and default analysis
    sel_string = "segid MEMB"
    ba = BilayerAnalyzer(
        structure='../pybilt/sample_bilayer/sample_bilayer.psf',
        trajectory='../pybilt/sample_bilayer/sample_bilayer_10frames.dcd',
        selection=sel_string,
    )
    # add analysis by string input
    ba.add_analysis('msd msd_2 leaflet upper resname POPC')
    print('Doing analysis iteration...')
    for _frame in ba:
        print(ba.reps['com_frame'])
        print(" ")


if __name__ == '__main__':
    print("Testing the various input options:")
    test_input_options()
    print("")
    print("Testing the analysis defaults:")
    test_analysis_defaults()
    print("")
    print("Testing the analysis iterator:")
    test_analysis_iterator()
