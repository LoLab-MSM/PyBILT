
"""Analysis Protocols

This is a support module that defines a set of classes used to contruct the 'analysis protocol' and the 'analysis' used
in the anlysis implemented by the pybilt.bilayer_analyzer.bilayer_anlayzer.BilayerAnalyzer class.
The AnalysisProtocol is the class used to organize and initial the the individual analysis functions/protocols. The
individual analysis functions/protocols are derived classes of AnalysisProtocol, and the availabel analysis can be
extended by simply defining a new protocol.
Example:
    # define a new analysis 'my_analysis'
    # add it to the valid_analysis list
    valid_analysis.append('my_analysis')
    #add it to analysis_obj_name_dict dictionary with the buildable object (e.g. mda_frame)  needed for the analysis.
    #The buildables are built by the bilayer analyzer and can accessed from it for the analysisanalysis.
    analysis_obj_name_dict['my_analysis'] = 'mda_frame'
    class MyAnalysisProtocol(AnalysisFunctionProtocol):
        #minimal def the __init__ and run_analysisfuncions
        def __init__(self, args):
            #define the initialization
        #The run_analysis function should always take three inputs
        #    of data that is passed in from the BilayerAnalayzer instance
        #    that owns this instance of an analysis protocol. The values are:
        #       ba_settings = BilayerAnalyzer.settings --> The adjustable
        #           settings dictionary which contains data like 'lateral', and
        #           'norm'.
        #       ba_reps = BilayerAnalyzer.reps --> The dictionary container of
        #           frame representations. e.g. 'com_frame' and 'lipid_grid'
        #       ba_mda_data = BilayerAnalyzer.mda_data --> An MDAData instance
        #           that contains data like the MDAnalysis universe, trajectory,
        #           bilayer_selection, etc.
        def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
            #.............
            #doing my analysis
            #................
            return

        #redefine any other functions from AnalysisProtocol or add new ones.


"""

# imports
import scipy.constants as scicon
import numpy as np
#import ast
try:
    import cPickle as pickle
except ImportError as error:
    import pickle
import warnings
import sys
import cmath
#range/xrange fix
if sys.version_info < (3,0):
    def range(*args, **kwargs):
        return xrange(*args, **kwargs)

# PyBILT imports
from pybilt.common.running_stats import RunningStats, binned_average
import pybilt.mda_tools.mda_density_profile as mda_dp
import pybilt.lipid_grid.lipid_grid_curv as lgc
from pybilt.common import distance_cutoff_clustering as dc_cluster
from pybilt.common.distance_cutoff_clustering import distance_euclidean_pbc
from scipy.spatial.distance import cdist
#need some containers for bookkeeping
command_protocols = {}
valid_analysis= []
analysis_obj_name_dict = {}
# obj_dict = {"com_frame":COMFrame}
#define the buildable objects that are used by the analysisfunctions
use_objects = {"mda_frame": True, "com_frame": True, "lipid_grid": False,
    "vector_frame": False}


# TO DO:
#


def word_list_to_string(word_list, delimeter=" "):
    """Creates a single string from a list of strings

    This function can be used to combine words in a list into one long sentence
    string.

    Args:
        word_list (list/tuple): A list (or other container) of strings.
        delimeter (str, Optional): A string to delimit the strings in the list
            when combining the strings.

    Returns:
        A string.
    """
    string = ""
    for word in word_list:
        string+=word+delimeter
    nchar = len(string)
    return str(string[0:nchar-1])


#def _run_analysis_alias(protocol_analyzer):
#    protocol = protocol_analyzer[0]
#    protocol.run_analysis(protocol_analyzer[1])
#    return protocol

def _run_analysis_alias(protocol_analyzer):
    """ An alias function to pass to multiprocessing threads

    This function is used internally in the BilayerAnalyzer.run_analysis_mp
    function to pass to the multiprocessing threads.
    """
    print(protocol_analyzer)
    return protocol_analyzer

# protocol for the analysis to run during the frame loop
class Analyses(object):
    """A class to facilitate analysis of the bilayers via the BilayerAnalyzer

    This object stores all the analyses that are being performed  by the
    BilayerAnalyzer class, and it provides functionality to add and remove
    analyses.

    Attributes:
        use_objects (dict): A dictionary of buildable objects that need to be
            contructed in the BilayerAnalyzer for the analysis defined in this
            protocol.
        in_commands (list): A list of the input strings for the analysis to be
            used.
        arguments (list): A list of the arguments for analysis.
        analysis_keys (list): A list of the keys assigned to analysis.
        command_protocol (dict): A dictionary of the analysis objects.
        analysis_ids (list): A list of the ids assigned to analysis.
        n_commands (int): The number of initialized analysis.

    """
    def __init__(self, analysis_commands):
        """Inits Analyses with input analysis_commands

        Args:
            analysis_commands (list): A list of the input strings for the
                analyses to be used. This is internally parsed together by
                the calling BilayerAnalyzer class object.

        """
        self.use_objects = use_objects
        self.in_commands = analysis_commands
        self.arguments = []
        self.analysis_keys = []
        self.analysis_ids = []
        self.command_protocol = dict()
        self.n_commands = 0
        # check analysis

        for command in analysis_commands:
            self.add_analysis(command)
        # object dependencies
        if self.use_objects['lipid_grid']:
            self.use_objects['com_frame'] = True
        return

    def __getitem__(self, item):
        """ Define the getitem function """
        return self.command_protocol[item]

    def __len__(self):
        """ Define the len function """
        return len(self.analysis_ids)

    def add_analysis(self, inputs):
        """Used to add a new analysis to the internal set of analyses

        Args:
            inputs (str, list, tuple, or dict): The input to parsed for the
                type of analysis and its settings to be added to set of
                analyses.
        """
        if isinstance(inputs, (str, basestring)):
            self._add_analysis_from_string(inputs)
        elif isinstance(inputs, (list, tuple)):
            self._add_analysis_from_list(inputs)
        elif isinstance(inputs, dict):
            self._add_analysis_from_dict(inputs)
        return

    def _add_analysis_from_string(self, analysis_string):
        """ Parses string inputs. """
        command = analysis_string.split()
        comp_key = command[0]
        comp_id = command[1]
        comp_args = word_list_to_string(command[1:])
        if (comp_key in valid_analysis):
            if (len(comp_args) >= 1):
                if comp_id not in self.analysis_ids:
                    comp_object = analysis_obj_name_dict[comp_key]
                    self.use_objects[comp_object] = True
                    self.arguments.append(comp_args)
                    self.analysis_ids.append(comp_id)
                    self.analysis_keys.append(comp_key)
                    self.command_protocol[comp_id] = command_protocols[
                        comp_key](comp_args)
                else:
                    raise RuntimeError("analysisid '{}' "
                                       "has already been used!".format(comp_id))
            else:
                raise RuntimeError("wrong number of arguments "
                                   "for analysis {}".format(command))
        else:
            raise RuntimeError("invalid analysisid"
                               " '{}' : {}".format(comp_key, command))
        self.n_commands += 1
        return

    def _add_analysis_from_list(self, analysis_list):
        """ Parses list/tuple inputs. """
        if len(analysis_list) == 3:
            comp_key = analysis_list[0]
            comp_id = analysis_list[1]
            comp_args = analysis_list[2]
            comp_args['analysis_id'] = comp_id
            if (comp_key in valid_analysis):
                if comp_id not in self.analysis_ids:
                    comp_object = analysis_obj_name_dict[comp_key]
                    self.use_objects[comp_object] = True
                    self.arguments.append(comp_args)
                    self.analysis_ids.append(comp_id)
                    self.analysis_keys.append(comp_key)
                    self.command_protocol[comp_id] = command_protocols[
                        comp_key](comp_args)
                else:
                    raise RuntimeError("analysis_id '{}' "
                                       "has already been used!".format(comp_id))
            else:
                raise RuntimeError("invalid analysis_key '{}'".format(comp_key))
        else:
            raise RuntimeError("wrong number of arguments in the input list/tuple. Should be three: [a_key, a_id, a_set_dict]")

        self.n_commands += 1
        return

    def _add_analysis_from_dict(self, analysis_dict):
        """ Parses dict inputs. """
        for key in ['analysis_key', 'analysis_id', "analysis_settings"]:
            if key not in analysis_dict.keys():
                raise RuntimeError("required key: {} , not in the input dictionary.".format(key))
        comp_key = analysis_dict['analysis_key']
        comp_id = analysis_dict['analysis_id']
        comp_args = analysis_dict['analysis_settings']
        comp_args['analysis_id'] = comp_id
        if (comp_key in valid_analysis):
            if comp_id not in self.analysis_ids:
                comp_object = analysis_obj_name_dict[comp_key]
                self.use_objects[comp_object] = True
                self.arguments.append(comp_args)
                self.analysis_ids.append(comp_id)
                self.analysis_keys.append(comp_key)
                self.command_protocol[comp_id] = command_protocols[
                    comp_key](comp_args)
            else:
                raise RuntimeError("analysisid '{}' "
                                   "has already been used!".format(comp_id))

        else:
            raise RuntimeError("invalid analysis_key '{}'".format(comp_key))
        self.n_commands+=1

        return

    def remove_analysis(self, analysis_id):
        """Removes the analysis with the given id from the set of analyses.

        Args:
            analysis_id (str): The string analysis_id of the analysis that is
                to be removed from the internal set of analyses.

        """
        if analysis_id in self.analysis_ids:
            del self.command_protocol[analysis_id]
            index = self.analysis_ids.index(analysis_id)
            del self.arguments[index]
            del self.analysis_keys[index]
            del self.analysis_ids[index]
            self.n_commands -= 1
        else:
            warnings.warn("no analysis with id '{}'".format(analysis_id))
        return

    def remove_all(self):
        """ Removes all of the analyses from internal set. """

        a_ids = self.analysis_ids
        for a_id in a_ids:
            self.remove_analysis(a_id)
        return

    def print_protocol(self):
        """ Prints to std out the protocol of the analyses. """

        print ('build objects:')
        for key in self.use_objects.keys():
            if self.use_objects[key]:
                print (key)
        print ("with analysis:")
        for analysis_id in self.analysis_ids:
            self.command_protocol[analysis_id].print_protocol()
        return

    def dump_data(self, path=None):
        """ Calls the individual save_data functions of each AnalysisProtocol.

        Args:
            path (str, Optional): The string path where output files should be
                dumped to disc.
        """
        print ('dumping analysis data to pickle files...')
        for analysis_id in self.analysis_ids:
            print ("analysis id: {} ---> {} ".format(
                analysis_id,
                self.command_protocol[analysis_id].save_file_name))
            self.command_protocol[analysis_id].save_data(path=path)

    def reset(self):
        """ Calls the individual rest functions for each AnalysisProtocol. """

        for analysis_id in self.analysis_ids:
            self.command_protocol[analysis_id].reset()
        return

#base class for analysis protocols
class AnalysisProtocol(object):
    """Base class for analysis protocols.

    Attributes:
        analysis_key (str): The key name of this analysis.
        analysis_id (str): The unique id assigned to this analyisis.
        save_file_name (str): The path and filename for the pickle file output of
            this analysis' results.
        settings (dict): A dict of the internal settings of the analysis.
        analysis_output (list or list like): Used to store the ouptut of this
            analysis during the frame loop.
    """
    _pickleable = True

    def __init__(self, args):
        """ Inits the AnalysisProtocol using the input args.

        Args:
            args (list): list of argument keys and values.
        """
        # required
        self._short_description = "parent analysis protocol"
        self._return_length = 1
        self.analysis_key = 'none'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        self._parse_args(args)
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        """Parses the setup arguments for this analysis.
        Args:
            args (list): List of argument keys and values.
        """
        # print args
        if isinstance(args, (basestring, str)):
            self._parse_string(args)
        elif isinstance(args, dict):
            self._parse_dict(args)
        else:
            raise RuntimeError("Invalid input structure for arguments/settings to analysis type "+self.analysis_key)
        return
        # required - a check protocol function which reports relevant settings
    def _parse_string(self, args):
        """Parses the input arguments from a string. """

        arg_dict = self._parse_str_to_dict(args)
        #type cast setttings if needed
        arg_dict = self._cast_settings(arg_dict)
        self._parse_dict(arg_dict)
        return

    def _parse_str_to_dict(self, args):
        """Converts the input argument str to dict. """
        arg_dict = dict()
        args = args.split()
        nargs = len(args)
        arg_dict['analysis_id'] = args[0]
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self._valid_settings:
                arg_dict[arg_key] =  arg_arg
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_key)
        return arg_dict

    # cast the input string values of settings to appropriate types -- should be overwritten in derived classes to
    # properly type cast their own settings.
    def _cast_settings(self, arg_dict):
        """Cast the input arguments from str to the appropriate type (e.g. int)."""

        for dummy_setting_key in arg_dict:
            pass
        return arg_dict

    def _parse_dict(self, args):
        """Parses the input arguments from a dict. """

        if 'analysis_id' not in args.keys():
            raise RuntimeError("required key \'anlaysis_id\' not assigned in input dict for analysis type: \'"+self.analysis_key+"\'")
        for arg_key in args.keys():
            arg_arg = args[arg_key]
            if arg_key in self._valid_settings:
                self.settings[arg_key] =  arg_arg
            elif arg_key == 'analysis_id':
                self.analysis_id = arg_arg
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_key)
        return

    def short_description(self):
        """Returns the protocols short description."""
        return self._short_description

    def pickleable(self):
        """Returns whether or not the protocol can be pickled."""
        return self._pickleable

    def print_protocol(self):
        """Prints to std out the internal data and settings for the protocol."""
        print ("Analysis: "+self._short_description)
        print ("  with analysis_id: {} ".format(self.analysis_id))
        print ("   and settings: ")
        for key in self.settings.keys():
            print ("    {}: {} ".format(key, self.settings[key]))
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        """Performs the analysis that this protocol represents for the current
        frame in the external BilayerAnalyzer's run_analysis function.

        Args:
            ba_settings (dict): The settings stored in the external
                BilayerAnalyzer instance.
            ba_reps (dict): The representation objects stored in the external
                BilayerAnalyzer instance.
            ba_mda_data (MDAData instance): The instance of MDAData stored in
                the external BilayerAnalyzer instance.
        """
        # do some stuff
        # get an output
        output = np.zeros(self._return_length)
        dummy_ba_settings = ba_settings
        dummy_ba_reps = ba_reps
        dummy_ba_mda_data = ba_mda_data
        # save the output
        self.analysis_output.append(output)
        return

    def save_data(self, path=None):
        """Dumps the outputs of this protocol to disc.

        Args:
            path (str, Optional): The string containing the path to the location
                that the analysis results should be dumped to on disc.
        """
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        with open(save_file, 'wb') as outfile:
            pickle.dump(np.array(self.analysis_output), outfile)

        return

    def get_data(self):
        """Returns the analysis_output of this protocol. """
        return np.array(self.analysis_output)

    def reset(self):
        """Resets the analysis by resetting the outputs and any necessary
        internal variables.
        """
        self.analysis_output = []
        return


# define a new analysis 'msd'
valid_analysis.append('msd')
analysis_obj_name_dict['msd'] = 'com_frame'


class MSDProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the mean squared displacement from a single time origin.

        The MSDProtocol is used to compute the mean squared displacement (MSD)
        of the centers of mass of the specified lipids for a single time origin
        (i.e. the first frame in the trajectory analysis) The MSD is given by
        MSD_i = <(r(t) - <r_0)**2>_i for lipid type i; the angle brackets denote
        averaging over all lipids of type i.

        This protocol is identified by the analysis key: 'msd'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', averages over all lipid types.

        References:
            1. Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical
                Properties of a Hydrated Lipid Bilayer from a Multinanosecond
                Molecular Dynamics Simulation, Biophysical Journal, Volume 81,
                Issue 5, 2001, Pages 2484-2494, ISSN 0006-3495,
                http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
                (http://www.sciencedirect.com/science/article/pii/S0006349501758948)

            2. Yoshimichi Andoh, Susumu Okazaki, Ryuichi Ueoka, Molecular
                dynamics study of lipid bilayers modeling the plasma membranes
                of normal murine thymocytes and leukemic GRSL cells, Biochimica
                et Biophysica Acta (BBA) - Biomembranes, Volume 1828, Issue 4,
                April 2013, Pages 1259-1270, ISSN 0005-2736,
                https://doi.org/10.1016/j.bbamem.2013.01.005.
                (http://www.sciencedirect.com/science/article/pii/S0005273613000096)

            3. Section 8.7,
                http://manual.gromacs.org/documentation/5.1.4/manual-5.1.4.pdf
        """
        # required
        self._short_description = "Single time origin mean squared displacement."
        self._return_length = 2
        self.analysis_key = 'msd'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        leaflet = self.settings['leaflet']
        group = self.settings['resname']
        if self._first_frame:
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            elif leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            self._indices = indices
        indices = self._indices
        n_com = len(indices)
        selcoords = np.zeros((n_com, 2))

        count = 0
        for i in indices:
            com_curr = ba_reps['com_frame'].lipidcom[i].com_unwrap[
                ba_settings['lateral']]
            selcoords[count] = com_curr[:]
            count += 1

        # initialize a numpy array to hold the msd for the selection
        msd = np.zeros(2)
        # initialize a running stats object to do the averaging over resids
        drs_stat = RunningStats()

        ref_coords = np.zeros((n_com, 2))
        if self._first_frame:
            count = 0
            for i in indices:
                com_curr = \
                    ba_reps['first_com_frame'].lipidcom[i].com_unwrap[
                        ba_settings['lateral']]
                ref_coords[count] = com_curr[:]
                count += 1
            self._ref_coords = ref_coords[:]
            self._first_frame = False
        else:
            ref_coords = self._ref_coords
            # get the current com frame list
        tc = ba_reps['com_frame'].time
        dt = tc
        dr = selcoords - ref_coords
        drs = dr * dr
        # loop over the selections for this frame
        for val in drs:
            drs_curr = val[:]
            drs_mag = drs_curr.sum()
            drs_stat.push(drs_mag)
        # get the msd for the current selection
        msdcurr = drs_stat.mean()
        msd[0] = dt
        msd[1] = msdcurr
        self.analysis_output.append(msd)
        return


# update the command_protocols dictionary
command_protocols['msd'] = MSDProtocol

# define a new analysis 'msd'
valid_analysis.append('msd_multi')
analysis_obj_name_dict['msd_multi'] = 'com_frame'


class MSDMultiProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the mean squared displacement using multiple time origins.

        The MSDMultiProtocol is used to compute the mean squared displacement
        (MSD) of the centers of mass of the specified lipids using multiple time
        origins and thus multiple time blocks.
        The MSD is given by
            MSD_i = <<(r(tau) - r_0)**2>_i>_tau
        for lipid type i; the inner angle
        brackets denote averaging over all lipids of type i and the
        outer brackets denote averaging over all time origins. The diffusion
        coefficient is estimated from the MSD using a simplified version
        Einstein's relation
            D_i ~ MSD_i/(4.0*tau)

        This protocol is identified by the analysis key: 'msd_multi'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', averages over all lipid types.
            n_tau (int): Specify the time block size in number of frames.
                Default: 50 ( 1000 picoseconds for timestep of 2 fs and frame
                output ever 100 timesteps).
            n_sigma (int): Specify the time between origins in number of frames.
                Default: 50.
        Note:
            According to results in Ref 3 the time blocks used to estimate the
                MSD should not be overlapping. Therefore, it is recommended to
                use n_sigma >= n_tau.

        References:
            1. Christofer Hofsab, Erik Lindahl, and Olle Edholm, "Molecular
                Dynamics Simulations of Phospholipid Bilayers with Cholesterol",
                Biophys J. 2003 Apr; 84(4): 2192-2206.
                doi:  10.1016/S0006-3495(03)75025-5

            2. Orsi, Mario, Julien Michel, and Jonathan W. Essex.
                "Coarse-grain modelling of DMPC and DOPC lipid bilayers."
                Journal of Physics: Condensed Matter 22.15 (2010): 155106.
                http://iopscience.iop.org/article/10.1088/0953-8984/22/15/155106/meta

            3. Gaurav Pranami and Monica H. Lamm, Estimating Error in Diffusion
                Coefficients Derived from Molecular Dynamics Simulations,
                Journal of Chemical Theory and Computation 2015 11 (10),
                4586-4592, DOI: 10.1021/acs.jctc.5b00574,
                http://pubs.acs.org/doi/full/10.1021/acs.jctc.5b00574
        """
        # required
        self._short_description = "Multiple time origin mean squared displacement."
        self._return_length = 2
        self.analysis_key = 'msd_multi'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self.settings['n_tau'] = 50
        self.settings['n_sigma'] = 50
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        self._tau_counter = 0
        self._sigma_counter = 0
        self._counters = [0]
        self._frame_counter = 0
        self._have_tau = False
        self._tau = 0.0
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'n_tau':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'n_sigma':
                    arg_dict[arg_key] = int(arg_arg)
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        leaflet = self.settings['leaflet']
        group = self.settings['resname']
        if self._first_frame:
            #determine the frame intervals
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            elif leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            self._indices = indices
            #self.analysis_output.append([])
        indices = self._indices
        n_com = len(indices)
        selcoords = np.zeros((n_com, 2))

        count = 0
        for i in indices:
            com_curr = ba_reps['com_frame'].lipidcom[i].com_unwrap[
                ba_settings['lateral']]
            selcoords[count] = com_curr[:]
            count += 1

        # initialize a running stats object to do the averaging over resids

        if self._first_frame:
            count = 0
            ref_coords = []
            for i in indices:
                com_curr = \
                    ba_reps['first_com_frame'].lipidcom[i].com_unwrap[
                        ba_settings['lateral']]
                ref_coords.append(com_curr[:])
                count += 1
            self._ref_coords = [np.array(ref_coords)]
            self.analysis_output = [0.0]
            self._first_frame = False
            # get the current com frame list
        #print(self._sigma_counter)
        #print(type(self.settings['n_sigma']))
        if self._sigma_counter < self.settings['n_sigma']:
            #print("{} {}".format(self._sigma_counter, self.settings['n_sigma']))
            #print("{}".format(self._sigma_counter < self.settings['n_sigma']))
            n_blocks = len(self.analysis_output)
            for b in range(self._tau_counter, n_blocks):
                if self._counters[b] < self.settings['n_tau']:
                    self._counters[b] += 1
                if self._counters[b] == self.settings['n_tau']:
                    drs_stat = RunningStats()
                    ref_coords = self._ref_coords[b]
                    dr = selcoords - ref_coords
                    drs = dr * dr
                    # loop over the selections for this frame
                    for val in drs:
                        drs_curr = val[:]
                        drs_mag = drs_curr.sum()
                        drs_stat.push(drs_mag)
                    # get the msd for the current selection
                    msdcurr = drs_stat.mean()
                    self.analysis_output[b] = msdcurr
                    self._tau_counter+=1
        elif self._sigma_counter == self.settings['n_sigma']:
            count = 0
            ref_coords = []
            for i in indices:
                com_curr = \
                    ba_reps['com_frame'].lipidcom[i].com_unwrap[
                        ba_settings['lateral']]
                ref_coords.append(com_curr[:])
                count += 1
            self._ref_coords.append(np.array(ref_coords))
            self.analysis_output.append(0.0)
            self._counters.append(0)
            n_blocks = len(self.analysis_output)
            for b in range(self._tau_counter, n_blocks):
                if self._counters[b] < self.settings['n_tau']:
                    self._counters[b] += 1
                if self._counters[b] == self.settings['n_tau']:
                    drs_stat = RunningStats()
                    ref_coords = self._ref_coords[b]
                    dr = selcoords - ref_coords
                    drs = dr * dr
                    # loop over the selections for this frame
                    for val in drs:
                        drs_curr = val[:]
                        drs_mag = drs_curr.sum()
                        drs_stat.push(drs_mag)
                    # get the msd for the current selection
                    msdcurr = drs_stat.mean()
                    self.analysis_output[b] = msdcurr
                    self._tau_counter+=1
            self._sigma_counter = 0

        self._sigma_counter += 1
        if not self._have_tau:
            self._frame_counter+=1
            if self._frame_counter == self.settings['n_tau']:
                self._tau = ba_reps['com_frame'].time \
                    - ba_reps['first_com_frame'].time
                self._have_tau = True
        return
    def _process_output(self, analysis_output):
        analysis_output = [val for val in analysis_output if val > 0.0]
        #print(analysis_output)
        msd_array = np.array(analysis_output)
        msd_average = msd_array.mean()
        msd_std_err = msd_array.std()/np.sqrt(len(msd_array))
        tau = self._tau
        diffusion_coeff = msd_average/(4.0*tau)
        diffusion_coeff_std_err = msd_std_err/(4.0*tau)
        analysis_output = np.array([msd_average, msd_std_err,
                                    diffusion_coeff, diffusion_coeff_std_err])
        return analysis_output

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name

        analysis_output = self._process_output(self.analysis_output)
        with open(save_file, 'wb') as outfile:
            pickle.dump(analysis_output, outfile)

        return

    def get_data(self):
        analysis_output = self._process_output(self.analysis_output)
        return analysis_output

    def reset(self):
        self.analysis_output = []
        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        self._tau_counter = 0
        self._sigma_counter = 0
        self._counters = [0]
        self._frame_counter = 0
        self._have_tau = False

        return

# update the command_protocols dictionary
command_protocols['msd_multi'] = MSDMultiProtocol

# define a new analysis 'apl_box'
valid_analysis.append('apl_box')
analysis_obj_name_dict['apl_box'] = 'mda_frame'


class APLBoxProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Estimate the area per lipid using the lateral area.

        The APLBoxProtocol is used to estimate the area per lipid (APL)
        using the lateral box dimensions. This approach is only accurate for
        homogenous lipid bilayers. If the bilayer is inhomogenous then tbis
        estimate represents a composite average of the area per lipid.

        This protocol is identified by the analysis key: 'apl_box'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            None

        References:
            1. Preston B. Moore, Carlos F. Lopez, Michael L. Klein, Dynamical Properties of a Hydrated Lipid Bilayer
                from a Multinanosecond Molecular Dynamics Simulation, Biophysical Journal, Volume 81, Issue 5, 2001,
                Pages 2484-2494, ISSN 0006-3495, http://dx.doi.org/10.1016/S0006-3495(01)75894-8.
                (http://www.sciencedirect.com/science/article/pii/S0006349501758948)
        """
        # required
        self._short_description = "Area per lipid using box dimensions."
        self._return_length = 4
        self.analysis_key = 'apl_box'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.running = RunningStats()
        self.analysis_output = []
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        box = ba_reps['current_mda_frame'].dimensions[0:3]
        plane = box[ba_settings['lateral']]
        area = plane[0] * plane[1] * 2.0
        nlipids = ba_mda_data.n_residues
        apl = area / nlipids
        time = ba_reps['current_mda_frame'].time
        self.running.push(apl)
        apl_t = np.zeros(self._return_length)
        apl_t[0] = time
        apl_t[1] = apl
        apl_t[2] = self.running.mean()
        apl_t[3] = self.running.deviation()
        self.analysis_output.append(apl_t)

        return

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        return


command_protocols['apl_box'] = APLBoxProtocol

# define a new analysis 'apl_box'
valid_analysis.append('bilayer_thickness')
analysis_obj_name_dict['bilayer_thickness'] = 'lipid_grid'


class BTGridProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Estimate the bilayer thickness using a gridding procedure.

        This protocol is identified by the analysis key: 'bilayer_thickness'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            None

        References:
            1. Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry
            2. Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858
        """
        # required
        self._short_description = "Bilayer thickness using lipid_grid."
        self._return_length = 4
        self.analysis_key = 'bilayer_thickness'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.running = RunningStats()
        self.analysis_output = []
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        current_thickness = ba_reps['lipid_grid'].average_thickness()
        # print (current_thickness)
        time = ba_reps['current_mda_frame'].time
        # print (time)
        self.running.push(current_thickness[0])
        bt_t = np.zeros(4)
        bt_t[0] = time
        bt_t[1] = current_thickness[0]
        bt_t[2] = self.running.mean()
        bt_t[3] = self.running.deviation()
        self.analysis_output.append(bt_t)

        return

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        return


command_protocols['bilayer_thickness'] = BTGridProtocol

# define a new analysis 'apl_grid'
valid_analysis.append('apl_grid')
analysis_obj_name_dict['apl_grid'] = 'lipid_grid'


class APLGridProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the indvidual area per lipid for each lipid type using a
         gridding procedure.

        This protocol is identified by the analysis key: 'apl_grid'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            None

        References:
            1. Allen et al. Vol. 30, No. 12 Journal of Computational Chemistry
            2. Gapsys et al. J Comput Aided Mol Des (2013) 27:845-858
        """
        # required
        self._short_description = "Area per lipid using lipid_grid"
        self._return_length = 4
        self.analysis_key = 'apl_grid'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        # parse input arguments
        self._parse_args(args)

        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.running = RunningStats()
        self.analysis_output = {}
        self.first_comp = True
        self.running_res = {}
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        apl_grid_out = ba_reps['lipid_grid'].area_per_lipid()
        # print (current_thickness)
        time = ba_reps['current_mda_frame'].time
        # print (time)
        self.running.push(apl_grid_out[0])
        if self.first_comp:
            self.analysis_output['composite'] = []
            for key in apl_grid_out[1].keys():
                self.running_res[key] = RunningStats()
                self.first_comp = False
                self.analysis_output[key] = []

        apl_t = np.zeros(4)
        apl_t[0] = time
        apl_t[1] = apl_grid_out[0]
        apl_t[2] = self.running.mean()
        apl_t[3] = self.running.deviation()

        self.analysis_output['composite'].append(apl_t)

        for key in apl_grid_out[1].keys():
            self.running_res[key].push(apl_grid_out[1][key][0])
            apl_res = np.zeros(4)
            apl_res[0] = time
            apl_res[1] = apl_grid_out[1][key][0]
            apl_res[2] = self.running_res[key].mean()
            apl_res[3] = self.running_res[key].deviation()
            self.analysis_output[key].append(apl_res)

        return

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name

        for key in self.analysis_output.keys():
            self.analysis_output[key] = np.array(self.analysis_output[key])
        with open(save_file, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        for key in self.analysis_output.keys():
            self.analysis_output[key] = np.array(self.analysis_output[key])
        return self.analysis_output

    def reset(self):
        self.running.reset()
        self.analysis_output = {}
        self.first_comp = True
        self.running_res = {}
        return


command_protocols['apl_grid'] = APLGridProtocol

# define a new analysis 'disp_vec'
valid_analysis.append('disp_vec')
analysis_obj_name_dict['disp_vec'] = 'com_frame'

#need to think more about box scaling (settings['scale']). currently if set True
# will scale by the box size of the reference frame
class DispVecProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Comute displacement vectors for each lipid in the specified
        leaflet(s) of bilayer.

        This protocol is identified by the analysis key: 'disp_vec'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', includes all lipid types.
            wrapped (bool): Specify whether to use the wrapped ('True') or
                un-wrapped ('False') coordintes for the base of the vectors.
                Default: False
            interval (int): Sets the frame interval over which to compute the
                    displacement vectors. f
            scale (bool): Specify whether to scale the coordinates by the box
                dimensions of the reference frame. Default: False

        References:
            1. Emma Falck, Tomasz Rog, Mikko Karttunen, and Ilpo Vattulainen,
                Lateral Diffusion in Lipid Membranes through Collective Flows,
                Journal of the American Chemical Society, 2008 130 (1), 44-45
                DOI: 10.1021/ja7103558
        """
        # required
        self._short_description = "Displacement vectors."
        self._return_length = 4
        self.analysis_key = 'disp_vec'
        self.analysis_id = 'none'

        #default settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self.settings['wrapped'] = False
        self.settings['interval'] = 5
        self.settings['scale'] = False
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.group = 'all'
        #self.wrapped = False
        #self.interval = 10
        # parse input arguments
        self._parse_args(args)

        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    arg_dict[arg_key] = arg_arg
                elif arg_key == 'scale':
                    arg_arg = arg_arg in ['True', 'true']
                    arg_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict


    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        if self.first_comp:
            self.last_com_frame = ba_reps['com_frame']
            self.first_comp = False
            self.last_frame = ba_settings['frame_range'][0]
            return
        current_frame = ba_reps['current_mda_frame'].frame
        #print(self.settings['leaflet'])
        #print(self.settings['scale'])
        interval = (current_frame) - (self.last_frame)
        #print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = ba_reps['leaflets']['upper']
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = ba_reps['leaflets']['lower']
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                warnings.warn(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            n_com = len(indices)

            # get the current frame
            curr_frame = ba_reps['com_frame']
            prev_frame = self.last_com_frame

            #get box dimensions for reference frame (i.e. prev_frame)
            box = prev_frame.box
            box_lateral = box[ba_settings['lateral']]
            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com, 4))
            # vec_ends = []
            count = 0
            resnames = []
            for i in indices:
                resname = curr_frame.lipidcom[i].type
                resnames.append(resname)
                com_i = curr_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j_w = prev_frame.lipidcom[i].com[ba_settings['lateral']]
                if self.settings['scale']:
                    #print("scaling coordinates..")
                    #print(self.settings['leaflet'])
                    #print(self.settings['scale'])
                    #print(type(self.settings['scale']))
                    #print(type(self.settings['wrapped']))
                    #quit()
                    com_i[0]/=box_lateral[0]
                    com_i[1]/=box_lateral[1]
                    com_j[0]/=box_lateral[0]
                    com_j[1]/=box_lateral[1]
                    com_j_w[0]/=box_lateral[0]
                    com_j_w[1]/=box_lateral[1]

                if self.settings['wrapped']:
                    vec_ends[count, 0] = com_j_w[0]
                    vec_ends[count, 1] = com_j_w[1]
                else:
                    vec_ends[count, 0] = com_j[0]
                    vec_ends[count, 1] = com_j[1]
                vec_ends[count, 2] = com_i[0] - com_j[0]
                vec_ends[count, 3] = com_i[1] - com_j[1]

                #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count += 1
            self.analysis_output.append([vec_ends, resnames])
            self.last_com_frame = ba_reps['com_frame']
            self.last_frame = current_frame
            return
        return

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name

        with open(save_file, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        return self.analysis_output


command_protocols['disp_vec'] = DispVecProtocol

# define a new analysis
valid_analysis.append('mass_dens')
analysis_obj_name_dict['mass_dens'] = 'mda_frame'


class MassDensProtocol(AnalysisProtocol):
    _pickleable = False
    def __init__(self, args):
        """Estimate the mass density profile for the specified selection.

        This protocol is used to estimate the 1-dimensional mass density profile
        for a selection of atoms along the bilayer normal. The profile is
        automatically centered on the bilayer's center of mass along the
        bilayer normal.

        This protocol is identified by the analysis key: 'mass_dens'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            selection_string (str): Provide the MDAnalysis compatible selection
                for the atoms to include in this analysis. Default: 'BILAYER',
                use all the lipids of the bilayer as recovered from the
                selection given to the external BilayerAnalyzer.
            n_bins (int): Set the number of bins to divide the normal dimensions
                into for binning.

        References:
            1. Needed!
        """

        # required
        self._short_description = "Mass density profile."
        self._return_length = None
        self.analysis_key = 'mass_dens'
        self.analysis_id = 'none'

        # default settings
        self.settings = dict()
        self.settings['selection_string'] = "BILAYER"
        self.settings['n_bins'] = 25
        self._valid_settings = self.settings.keys()
        #self.selection_string = 'all'
        #self.n_bins = 25

        #parse input arguments/settings
        self._parse_args(args)

        self.save_file_name = self.analysis_id + ".pickle"


        # storage for output
        self.centers = None
        self.n_frames = 0
        self.analysis_output = []
        self.first_comp = True
        self.selection = None

        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):
        #in this case the casting is taken care of in _parse_str_to_dict
        return arg_dict

    # needed to overwrite the string parser to handle the selection string
    def _parse_str_to_dict(self, args):
        # print args
        arg_dict = dict()
        args = args.split()
        nargs = len(args)
        arg_dict['analysis_id'] = args[0]
        read_sel_string = False
        n_bins_arg = False
        for i in range(1, nargs, 1):
            arg_key = args[i]

            #print ('arg_key: ', arg_key)
            if arg_key in self._valid_settings:
                if arg_key == 'n_bins':
                    arg_arg = args[i + 1]
                    arg_dict['n_bins'] = int(arg_arg)
                    read_sel_string = False
                    n_bins_arg = True
                elif arg_key == 'selection':
                    selection_words = [args[j] for j in range(i + 1, nargs) if
                                       (args[j] not in self._valid_settings)]
                    i += len(selection_words)
                    selection_string = ""
                    for word in selection_words:
                        selection_string += " " + word
                    arg_dict['selection_string'] = selection_string
                    read_sel_string = True
                    n_bins_arg = False
            elif not (read_sel_string or n_bins_arg):
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_id)
        return arg_dict


    def reset(self):
        self.centers = None
        self.n_frames = 0
        self.analysis_output = []
        self.first_comp = True
        self.selection = None

        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        first = self.first_comp
        if self.first_comp:
            print self.settings['selection_string']
            if self.settings['selection_string'] == "BILAYER":
                print("bilayer_sel")
                self.selection = ba_mda_data.bilayer_sel
            else:
                print("non-bilayer")
                self.selection = ba_mda_data.mda_universe.select_atoms(
                    self.settings['selection_string'])
            self.first_comp = False

        # print "there are ",len(indices)," members"

        norm_axis = ba_settings['normal_dimension']
        ref_sel = ba_mda_data.bilayer_sel
        # mda_density_profile
        centers_density = mda_dp.mass_density_profile(
            ba_mda_data.mda_trajectory, self.selection,
            fstart=ba_settings['frame_index'],
            fend=ba_settings['frame_index'] + 1, axis=norm_axis,
            nbins=self.settings['n_bins'], refsel=ref_sel)
        if first:
            self.centers = centers_density[0]
            self.analysis_output = np.zeros(len(self.centers))
        self.analysis_output += centers_density[1]
        self.n_frames += 1
        return

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        centers_density = (self.centers, self.analysis_output / self.n_frames)
        with open(save_file, 'wb') as outfile:
            pickle.dump(centers_density, outfile)

        return

    def get_data(self):
        return (self.centers, self.analysis_output / self.n_frames)

    def __getstate__(self):
        odict = dict()
        for key in self.__dict__.keys():
            odict[key] = self.__dict__[key]
        odict['selection'] = None
        return odict

    def __setstate__(self, in_dict):
        self.__dict__ = in_dict


command_protocols['mass_dens'] = MassDensProtocol

# define a new analysis 'nnf'
valid_analysis.append('nnf')
analysis_obj_name_dict['nnf'] = 'com_frame'

#Found in: M. Orsi and J. W. Essex, Faraday Discuss., 2013, 161, 249-272
#Originally described in: A. H. de Vries, A. E. Mark and S. J. Marrink, J. Phys. Chem. B, 2004, 108, 2454-2463
#similar method 'fractional interactions' used in: Koldso H, Shorthouse D, He Lie Sansom MSP (2014) "Lipid Clustering
# Correlates with Membrane Curvature as Revealed by Molecular Simulations of Complex Lipid Bilayers." PloS Comput Biol
# 10(10): e1003911. doi.10.1371/journal.pcbi.1003911
class NNFProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the nearest neighbor fraction for one lipid type with
        another.

        This analysis picks a specified number (n_neighbors) of nearest
        neighbors centered on a lipid of reference lipid type and then counts
        the number of lipids (M) of target lipid type and estimates the
        fraction, nnf = M/n_neighbors. This metric is also referred to as
        fractional interations.

        This protocol is identified by the analysis key: 'nnf'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname_1 (str): Specify the resname of the reference lipid type to
                include in this analysis. Default: 'first', the first lipid in
                the list pulled from the com_frame representation.
            resname_2 (str): Specify the resname of the target lipid type to
                include in this analysis. Default: 'first', the first lipid in
                the list pulled from the com_frame representation.

        References:
            1. A. H. de Vries, A. E. Mark and S. J. Marrink, J. Phys. Chem. B,
                2004, 108, 2454-2463

            2. M. Orsi and J. W. Essex, Faraday Discuss., 2013, 161, 249-272

            3. Koldso H, Shorthouse D, He Lie Sansom MSP (2014) "Lipid
                Clustering Correlates with Membrane Curvature as Revealed by
                Molecular Simulations of Complex Lipid Bilayers." PloS Comput
                Biol 10(10): e1003911. doi.10.1371/journal.pcbi.1003911
        """
        # required
        self._short_description = "Lateral order nearest neighbor fraction."

        self._return_length = 2
        self.analysis_key = 'nnf'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['n_neighbors'] = 5
        self.settings['resname_1'] = 'first'
        self.settings['resname_2'] = 'first'
        self._valid_settings = self.settings.keys()
        #parse input arguments/settings
        self._parse_args(args)

        self.save_file_name = self.analysis_id + ".pickle"

        # for outputs
        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        self.running = RunningStats()
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'n_neighbors':
                    arg_dict[arg_key] = int(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        self._first_frame = True
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        do_leaflet = []
        if self.settings['leaflet'] == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.settings['leaflet']]

        if self._first_frame:
            #pass
            # build group/resname/lipid type list
            lipid_types = []
            nlipids = 0
            for leaflet_name in do_leaflet:
                leaflet = ba_reps['leaflets'][leaflet_name]
                groups = leaflet.get_group_names()
                nlipids += len(leaflet.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.settings['resname_1'] == 'first':
                self.settings['resname_1'] =  lipid_types[0]
            if self.settings['resname_2'] == 'first':
                self.settings['resname_2'] = lipid_types[0]
            if self.settings['n_neighbors'] > nlipids:
                self.settings['n_neighbors'] = nlipids-1

            n_ltypes = len(lipid_types)
            self.lipid_types = lipid_types
            self.n_ltypes = n_ltypes

        lipid_types = self.lipid_types
        n_ltypes = self.n_ltypes
        #print lipid_types
        x_index = ba_settings['lateral'][0]
        y_index = ba_settings['lateral'][1]
        box = ba_reps['current_mda_frame'].dimensions[0:3]
        box_x = box[x_index]
        box_y = box[y_index]
        l_box = np.array([box_x, box_y])
        box_x_h = box_x / 2.0
        box_y_h = box_y / 2.0
        ltype_a = self.settings['resname_1']
        ltype_b = self.settings['resname_2']
        #print "ltype_a: ",ltype_a," ltype_b: ",ltype_b
        avg_frac = RunningStats()
        for leaflet_name in do_leaflet:
            leaflet = ba_reps['leaflets'][leaflet_name]

            if leaflet.has_group(ltype_a) and leaflet.has_group(ltype_b):
                ltype_a_indices = leaflet.get_group_indices(ltype_a)
                all_index = leaflet.get_member_indices()
                for i in ltype_a_indices:
                    neighbors = []
                    for j in all_index:
                        if i != j:
                            pos_a = ba_reps['com_frame'].lipidcom[i].com[[x_index, y_index]]
                            pos_b = ba_reps['com_frame'].lipidcom[j].com[[x_index, y_index]]
                            dist = distance_euclidean_pbc(pos_a, pos_b, l_box,
                                                          center='box_half')
                            ltype = ba_reps['com_frame'].lipidcom[j].type
                            #print "ltype: ",ltype," dist ",dist
                            neighbors.append([j, dist, ltype])
                    neighbors.sort(key=lambda x: x[1])
                    nn_neighbors = neighbors[0:self.settings['n_neighbors']]
                    #print neighbors
                    #print nn_neighbors
                    #quit()
                    n_type_b = 0.0
                    for neighbor in nn_neighbors:
                        ntype = neighbor[2]
                        if ntype == ltype_b:
                            n_type_b += 1.0
                    frac = n_type_b/self.settings['n_neighbors']
                    #print "frac ",frac," n_type_b: ",n_type_b," set_n_neighbors: ",self.settings['n_neighbors']

                    avg_frac.push(frac)

        f_current = avg_frac.mean()
        self.running.push(f_current)
        f_run = self.running.mean()
        f_run_dev = self.running.deviation()
        tc = ba_reps['com_frame'].time
        self.analysis_output.append([tc, f_current, f_run, f_run_dev])
        return

    # def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
    #
    #     if self._first_frame:
    #         pass
    #         #build group/resname/lipid type list
    #         lipid_types = []
    #         for leaflet in ba_reps['leaflets']:
    #             groups = leaflet.get_group_names()
    #             for group in groups:
    #                 if group not in lipid_types:
    #                     lipid_types.append(group)
    #         n_ltypes = len(lipid_types)
    #         self.lipid_types = lipid_types
    #         self.n_ltypes = n_ltypes
    #
    #     lipid_types = self.lipid_types
    #     n_ltypes = self.n_ltypes
    #
    #     x_index = ba_settings['lateral'][0]
    #     y_index = ba_settings['lateral'][1]
    #     box = ba_reps['current_mda_frame'].dimensions[0:3]
    #     box_x = box[x_index]
    #     box_y = box[y_index]
    #     box_x_h = box_x/2.0
    #     box_y_h = box_y/2.0
    #
    #     for leaflet in ba_reps['leaflets']:
    #         all_indices = leaflet.get_member_indices()
    #
    #         #X-X types
    #         for lipid_type in lipid_types:
    #             if leaflet.has_group(lipid_type):
    #                 ltype_indices = leaflet.get_group_indices(lipid_type)
    #                 nlip = len(ltype_indices)
    #                 neigbors = dict()
    #                 for lindex in ltype_indices:
    #                     neigbors[lindex] = []
    #                     for lindex_b in ltype_indices:
    #
    #                         if lindex_b != lindex:
    #                             pos_a = ba_reps['com_frame'].lipidcom[lindex].com
    #                             pos_b = ba_reps['com_frame'].lipidcom[lindex_b].com
    #                             dx = np.abs(pos_a[x_index] - pos_b[x_index])
    #                             dy = np.abs(pos_a[y_index] - pos_b[y_index])
    #                             #minimum image for wrapped coordinates
    #                             if dx > box_x_h:
    #                                dx = box_x - np.absolute(pos_a[x_index]-box_x_h) - np.absolute(pos_b[x_index]-box_x_h)
    #
    #                             if dy > box_yy_h:
    #                                dy = box_y - np.absolute(pos_a[y_index]-box_y_h) - np.absolute(pos_b[y_index]-box_y_h)
    #                             dist = np.sqrt(dx**2 + dy**2)
    #                             neigbors[lindex_b].append([lindex_b, dist])
    #
    #
    #
    #         #X-Y types
    #         for i in range(n_ltypes-1):
    #             for j in range(i+1, n_ltypes):
    #             ltype_i = lipid_types[i]
    #             ltype_j = lipid_types[j]
    #             if leaflet.has_group(ltype_i) and leaflet.has_group(ltype_j):
    #
    #                 ltype_indices = leaflet.get_group_indices()
    #                 for lipid_type_b in lipid_types:
    #                     if leaflet.has_group(lipid_type_b):
    #                         ltype_b_indices
    #     tc = ba_reps['com_frame'].time
    #
    #
    #     return


# update the command_protocols dictionary
command_protocols['nnf'] = NNFProtocol

# define a new analysis 'disp_vec'
valid_analysis.append('disp_vec_corr')
analysis_obj_name_dict['disp_vec_corr'] = 'com_frame'


class DispVecCorrelationProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Comute the pair-wise cross correlation matrix for the displacement
        vectors for each lipid in the specified leaflet(s) of bilayer.

        This analysis computes the displacement vectors as in DispVecProtocol,
        but then continues to compute the pair-wise cross correlation matrix
        between each vector. i.e. the cos(theta) for the angle theta between the
        vectors.

        This protocol is identified by the analysis key: 'disp_vec_corr'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', includes all lipid types.
            wrapped (bool): Specify whether to use the wrapped ('True') or
                un-wrapped ('False') coordintes for the base of the vectors.
                Default: False
            interval (int): Sets the frame interval over which to compute the
                    displacement vectors. f

        References:
            None
        """
        # required
        self._short_description = "Displacement vector correlation matrix."

        self._return_length = 4
        self.analysis_key = 'disp_vec_corr'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self.settings['wrapped'] = False
        self.settings['interval'] = 5
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    arg_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        if self.first_comp:
            self.last_com_frame = ba_reps['com_frame']
            self.first_comp = False
            return
        current_frame = ba_reps['current_mda_frame'].frame

        interval = (current_frame) - (self.last_frame)
        #print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = ba_reps['leaflets']['upper']
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = ba_reps['leaflets']['lower']
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                warnings.warn(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            n_com = len(indices)

            # get the current frame
            curr_frame = ba_reps['com_frame']
            prev_frame = self.last_com_frame
            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com, 4))
            # vec_ends = []
            count = 0
            resnames = []
            for i in indices:
                resname = curr_frame.lipidcom[i].type
                resnames.append(resname)
                com_i = curr_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j_w = prev_frame.lipidcom[i].com[ba_settings['lateral']]
                if self.settings['wrapped']:
                    vec_ends[count, 0] = com_j_w[0]
                    vec_ends[count, 1] = com_j_w[1]
                else:
                    vec_ends[count, 0] = com_j[0]
                    vec_ends[count, 1] = com_j[1]
                vec_ends[count, 2] = com_i[0] - com_j[0]
                vec_ends[count, 3] = com_i[1] - com_j[1]
                #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count += 1
            corr_mat = np.zeros((n_com, n_com))
            #loop over vectors
            for i in range(n_com):
                corr_mat[i,i] = 1.0

            for i in range(n_com-1):
                vec_end_a = vec_ends[i]
                vec_a = vec_end_a[2:4] - vec_end_a[0:2]
                for j in range(i+1, n_com):
                    vec_end_b = vec_ends[j]
                    vec_b = vec_end_b[2:4] - vec_end_b[0:2]
                    dot = np.dot(vec_a, vec_b)
                    cos_t = dot/(np.linalg.norm(vec_a)*np.linalg.norm(vec_b))
                    corr_mat[i,j] = cos_t
                    corr_mat[j,i] = cos_t



            self.analysis_output.append([corr_mat, resnames])
            self.last_com_frame = ba_reps['com_frame']
            self.last_frame = current_frame
            #return vec_ends
        return

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name

        with open(save_file, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        return self.analysis_output


command_protocols['disp_vec_corr'] = DispVecCorrelationProtocol

# define a new analysis 'disp_vec'
valid_analysis.append('disp_vec_nncorr')
analysis_obj_name_dict['disp_vec_nncorr'] = 'com_frame'


class DispVecNNCorrelationProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Comute the pair-wise cross correlations for the displacement
        vectors for each lipid in the specified leaflet(s) of bilayer and its
        nearest neighbor.

        This analysis computes the displacement vectors as in DispVecProtocol,
        but then continues to compute the pair-wise cross correlation between
        each vector and its nearest neighbor. i.e. the cos(theta) for the angle
        theta between the vectors.

        This protocol is identified by the analysis key: 'disp_vec_nncorr'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', includes all lipid types.
            wrapped (bool): Specify whether to use the wrapped ('True') or
                un-wrapped ('False') coordintes for the base of the vectors.
                Default: False
            interval (int): Sets the frame interval over which to compute the
                    displacement vectors. f

        References:
            None
        """
        # required
        self._short_description = "Displacement vector nearest neigbor correlations."
        self._return_length = 4
        self.analysis_key = 'disp_vec_nncorr'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self.settings['wrapped'] = False
        self.settings['interval'] = 5
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    arg_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        if self.first_comp:
            self.last_com_frame = ba_reps['com_frame']
            self.first_comp = False
            return
        current_frame = ba_reps['current_mda_frame'].frame

        interval = (current_frame) - (self.last_frame)
        #print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = ba_reps['leaflets']['upper']
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = ba_reps['leaflets']['lower']
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                warnings.warn(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            n_com = len(indices)

            # get the current frame
            curr_frame = ba_reps['com_frame']
            prev_frame = self.last_com_frame
            #get the vector ends

            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com, 4))
            # vec_ends = []
            count = 0
            resnames = []
            for i in indices:
                resname = curr_frame.lipidcom[i].type
                resnames.append(resname)
                com_i = curr_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j_w = prev_frame.lipidcom[i].com[ba_settings['lateral']]
                if self.settings['wrapped']:
                    vec_ends[count, 0] = com_j_w[0]
                    vec_ends[count, 1] = com_j_w[1]
                else:
                    vec_ends[count, 0] = com_j[0]
                    vec_ends[count, 1] = com_j[1]
                vec_ends[count, 2] = com_i[0] - com_j[0]
                vec_ends[count, 3] = com_i[1] - com_j[1]
                #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count += 1
            corr = np.zeros((n_com,1))
            #loop over vectors
            for i in range(n_com):
                vec_end_a = vec_ends[i]
                vec_a = vec_end_a[2:4] - vec_end_a[0:2]
                #now loop over again and find nearest neighbor at vector bases
                nn_vec = []
                nn_dist = 1000000.0
                for j in range(n_com):
                    if j != i:
                        vec_end_b = vec_ends[j]
                        vec_b = vec_end_b[2:4] - vec_end_b[0:2]
                        dist = np.linalg.norm(vec_end_b[0:2]-vec_end_a[0:2])
                        if dist < nn_dist:
                            nn_dist = dist
                            nn_vec = vec_b
                #now analysis the correlation (cos(theta))
                dot = np.dot(vec_a, nn_vec)
                cos_t = dot/(np.linalg.norm(vec_a)*np.linalg.norm(nn_vec))
                corr[i] = cos_t

            self.analysis_output.append([corr, resnames])
            self.last_com_frame = ba_reps['com_frame']
            self.last_frame = current_frame
            #return vec_ends
        return

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name

        with open(save_file, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        return self.analysis_output


command_protocols['disp_vec_nncorr'] = DispVecNNCorrelationProtocol

# define a new analysis 'apl_box'
valid_analysis.append('ndcorr')
analysis_obj_name_dict['ndcorr'] = 'com_frame'

class NDCorrProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Correlation between bilayer surfucace curvature and the clustering of
        lipid molecules.

        This protocol is used to estimate the cross correlation between the
        normal dimension deflection of lipids and the lipid types in local
        blocks of the bilayer. This serves as a measure of the correlation
        between the local curvature and composition of the bilayer.

        This protocol is identified by the analysis key: 'ndcorr'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            None

        Note:
            Automatically estimates the correlations with each lipid type in the
            bilayer selection provided to the external BilayerAnalyzer object.

        References:
            1. Koldso H, Shorthouse D, He Lie Sansom MSP (2014) "Lipid
                Clustering Correlates with Membrane Curvature as Revealed by
                Molecular Simulations of Complex Lipid Bilayers." PloS Comput
                Biol 10(10): e1003911. doi.10.1371/journal.pcbi.1003911
        """
        # required
        self._short_description = "Normal dimension displacement-lipid type cross correlation."

        self._return_length = 4
        self.analysis_key = 'ndcorr'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        # default function settings

        # parse input arguments if given
        self._parse_args(args)
        #save file name for pickle dump of results
        self.save_file_name = self.analysis_id + ".pickle"

        #for analysis and outputs
        self._first_frame = True
        # storage for output
        self.analysis_output = dict()
        self.running_stats = dict()
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        #construct the grids
        grids = lgc.LipidGrids(ba_reps['com_frame'],ba_reps['leaflets'],ba_settings['lateral'])

        if self._first_frame:
            leafs = grids.leaf_grid.keys()
            all_types = []
            for leaf in leafs:
                self.analysis_output[leaf] = dict()
                self.running_stats[leaf] = dict()
                groups = grids.leaflets[leaf].get_group_names()
                for group in groups:
                    if group not in all_types:
                        all_types.append(group)
            for leaf in leafs:
                for l_type in all_types:
                    self.analysis_output[leaf][l_type] = []
                    self.running_stats[leaf][l_type] = RunningStats()
            self._first_frame = False
        #analysis the correlations
        correlations = grids.norm_displacement_cross_correlation()
        time = ba_reps['current_mda_frame'].time
        #extract the data
        leafs = correlations.keys()
        for leaf in leafs:
            for l_type in correlations[leaf].keys():
                corr = correlations[leaf][l_type]
                self.running_stats[leaf][l_type].push(corr)
                corr_run = self.running_stats[leaf][l_type].mean()
                corr_std = self.running_stats[leaf][l_type].deviation()
                self.analysis_output[leaf][l_type].append(np.array([time, corr, corr_run, corr_std]))

        return

    def reset(self):
        self._first_frame = True
        self.analysis_output = dict()
        self.running_stats = dict()

        return

    def get_data(self):
        for key in self.analysis_output.keys():
            for tkey in self.analysis_output[key].keys():
                self.analysis_output[key][tkey] = np.array(self.analysis_output[key][tkey])
        return self.analysis_output

command_protocols['ndcorr'] = NDCorrProtocol

# define a new analysis 'apl_box'
valid_analysis.append('dc_cluster')
analysis_obj_name_dict['dc_cluster'] = 'com_frame'

#
class DCClusterProtocol(AnalysisProtocol):
    def __init__(self, args):
        """A type of hiearchical clustering where points (lipid centers of mass)
        are are added to a cluster if they are within a specified distance
        of any other point within the cluster.

        This protocol is identified by the analysis key: 'dc_cluster'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', includes all lipid types.
            cutoff (float): The cutoff distance to use for the clustering.
                Default: 12.0

        Note:
            Only finds the self clusters for a single lipid type as specified by
            the 'resname' setting.

        References:
            None
        """
        # required
        self._short_description = "Distance cutoff clustering."

        self._return_length = 4
        self.analysis_key = 'dc_cluster'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['resname'] = 'first'
        self.settings['leaflet'] = 'upper'
        self.settings['cutoff'] = 12.0
        self._valid_settings = self.settings.keys()
        # default function settings

        # parse input arguments if given
        self._parse_args(args)
        #save file name for pickle dump of results
        self.save_file_name = self.analysis_id + ".pickle"

        #for analysis and outputs
        self._first_frame = True
        self.converted = False
        # storage for output
        self.analysis_output = dict()
        self.running_stats = dict()
        self.analysis_output['clusters'] = []
        self.analysis_output['nclusters'] = []
        self.analysis_output['max_size'] = []
        self.analysis_output['min_size'] = []
        self.analysis_output['avg_size'] = []
        self.running_stats['nclusters'] = RunningStats()
        self.running_stats['max_size'] = RunningStats()
        self.running_stats['min_size'] = RunningStats()
        self.running_stats['avg_size'] = RunningStats()
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'cutoff':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        do_leaflet = []
        if self.settings['leaflet'] == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.settings['leaflet']]

        if self._first_frame:
            #pass
            # build group/resname/lipid type list
            lipid_types = []
            for leaflet_name in do_leaflet:
                leaflet = ba_reps['leaflets'][leaflet_name]
                groups = leaflet.get_group_names()
        #        nlipids += len(leaflet.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.settings['resname'] == 'first':
                self.settings['resname'] =  lipid_types[0]
            self._first_frame = False
        indices = []
        for leaflets in do_leaflet:
            curr_leaf = ba_reps['leaflets'][leaflets]
            indices += curr_leaf.get_group_indices(self.settings['resname'])
        pos = []
        for index in indices:
            pos.append(ba_reps['com_frame'].lipidcom[index].com[ba_settings['lateral']])
        pos = np.array(pos)
        box = ba_reps['com_frame'].box[ba_settings['lateral']]
        dist_func = dc_cluster.distance_euclidean_pbc
        clusters = dc_cluster.distance_cutoff_clustering(pos, self.settings['cutoff'], dist_func, 1, box, center='box_half')

        nclusters = len(clusters)
        min_size = 0
        max_size = 0
        avg_size = 0.0
        for cluster in clusters:
            n = len(cluster)
            if n > max_size:
                max_size = n
            elif n > 0 and max_size == 0:
                min_size = n
            elif max_size > 0 and n < max_size:
                min_size = n
            avg_size += n
        #protect for divide by zero
        if nclusters > 0:
            avg_size /= nclusters

        self.running_stats['nclusters'].push(nclusters)
        self.running_stats['max_size'].push(max_size)
        self.running_stats['min_size'].push(min_size)
        self.running_stats['avg_size'].push(avg_size)
        time = ba_reps['com_frame'].time
        self.analysis_output['clusters'].append(clusters)
        self.analysis_output['nclusters'].append([time,nclusters,self.running_stats['nclusters'].mean(), self.running_stats['nclusters'].deviation()])
        self.analysis_output['max_size'].append([time,max_size,self.running_stats['max_size'].mean(), self.running_stats['max_size'].deviation()])
        self.analysis_output['min_size'].append([time,min_size,self.running_stats['min_size'].mean(), self.running_stats['min_size'].deviation()])
        self.analysis_output['avg_size'].append([time,avg_size,self.running_stats['avg_size'].mean(), self.running_stats['avg_size'].deviation()])

        return

    def reset(self):
        self._first_frame = True

        for key in self.analysis_output.keys():
            self.analysis_output[key] = []
        for key in self.running_stats.keys():
            self.running_stats[key].reset()

        self.converted = False
        return

    def get_data(self):
        if not self.converted:
            for key in self.analysis_output.keys():
                if key != 'clusters':
                    self.analysis_output[key] = np.array(self.analysis_output[key])
            self.converted = True
        return self.analysis_output

command_protocols['dc_cluster'] = DCClusterProtocol

# define a new analysis 'apl_box'
valid_analysis.append('vcm')
analysis_obj_name_dict['vcm'] = 'mda_frame'

class VolumeCompressibilityModulusProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Estimate the isothermal volume compressibility modulus.

        This protocol is used to estimate the volume compressibility modulus,
            K_V = (<V>kT) / var(V),
        where V is the volume.

        This protocol is identified by the analysis key: 'vcm'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            temperature (float): The absolute temperature that the simulation
                was run at (i.e. in Kelvin). Default: 298.15 K

        References:
            1. Christofer Hofsab, Erik Lindahl, and Olle Edholm, "Molecular
                Dynamics Simulations of Phospholipid Bilayers with Cholesterol",
                Biophys J. 2003 Apr; 84(4): 2192-2206.
                doi:  10.1016/S0006-3495(03)75025-5
        """
        # required
        self._short_description = "Volume compressibility modulus."
        self._return_length = 4
        self.analysis_key = 'vcm'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['temperature'] = 298.15
        self._valid_settings = self.settings.keys()

        # parse input arguments
        self._parse_args(args)

        #output filename for pickle dump of results
        self.save_file_name = self.analysis_id + ".pickle"


        # storage for output
        self.n_frames = 0
        self.volume_run = RunningStats()
        self.analysis_output = []
        self.first_comp = True
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'temperature':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict


    def reset(self):
        self.n_frames = 0
        self.volume_run.reset()
        self.analysis_output = []
        self.first_comp = True
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        dimensions = ba_reps['current_mda_frame'].dimensions[0:3]
        volume = dimensions.prod()
        #print(area)
        self.volume_run.push(volume)

        Kv = 1.0
        if self.first_comp:
            self.first_comp = False
            self.n_frames += 1
            # return
            Kv = 0.0
        else:
            Kv = ( (self.volume_run.mean() * scicon.k *
                self.settings['temperature'])
                / self.volume_run.deviation() ** 2)
        time = ba_reps['current_mda_frame'].time
        self.analysis_output.append([time, Kv])
        self.n_frames += 1
        return

command_protocols['vcm'] = VolumeCompressibilityModulusProtocol

# define a new analysis 'apl_box'
valid_analysis.append('acm')
analysis_obj_name_dict['acm'] = 'mda_frame'

class AreaCompressibilityModulusProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Estimate the isothermal area compressibility modulus.

        This protocol is used to estimate the area compressibility modulus,
            K_A = (<A>kT) / var(A),
        where A is the area in the lateal dimension of the bilayer.

        This protocol is identified by the analysis key: 'acm'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            temperature (float): The absolute temperature that the simulation
                was run at (i.e. in Kelvin). Default: 298.15 K

        References:
            1. Christofer Hofsab, Erik Lindahl, and Olle Edholm, "Molecular
                Dynamics Simulations of Phospholipid Bilayers with Cholesterol",
                Biophys J. 2003 Apr; 84(4): 2192-2206.
                doi:  10.1016/S0006-3495(03)75025-5
            2. L. Janosi and A. A. Gorfe, J. Chem. Theory Comput. 2010, 6,
                3267-3273
            3. D. Aguayo, F. D. Gonzalez-Nilo, and C. Chipot, J. Chem. Theory
                Comput. 2012, 8, 1765-1773
        """
        # required
        self._short_description = "Area compressibility modulus."
        self._return_length = 4
        self.analysis_key = 'acm'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['temperature'] = 298.15
        self._valid_settings = self.settings.keys()

        # parse input arguments
        self._parse_args(args)

        #output filename for pickle dump of results
        self.save_file_name = self.analysis_id + ".pickle"


        # storage for output
        self.n_frames = 0
        self.area_run = RunningStats()
        self.analysis_output = []
        self.first_comp = True
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'temperature':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict


    def reset(self):
        self.n_frames = 0
        self.area_run.reset()
        self.analysis_output = []
        self.first_comp = True
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        dimensions = ba_reps['current_mda_frame'].dimensions[0:3]
        area = dimensions.prod()/dimensions[ba_settings['norm']]
        #print(area)
        self.area_run.push(area)
        Ka = 1.0
        if self.first_comp:
            self.first_comp = False
            self.n_frames+=1
            #return
            Ka = 0.0
        else:
            Ka = (self.area_run.mean()*scicon.k*self.settings['temperature'])/self.area_run.variance()
        #print "<A>: ", self.area_run.mean(), " var(A): ",self.area_run.variance()
        #print "T: ",self.settings['temperature']," k: ",scicon.k
        #conversion factor for Joules/Angstrom^2 to milliNewtons/meter
        Ka*=10.0**23

        time = ba_reps['current_mda_frame'].time
        self.analysis_output.append([time, Ka])
        self.n_frames += 1
        return

command_protocols['acm'] = AreaCompressibilityModulusProtocol

# define a new analysis 'msd'
valid_analysis.append('ald')
analysis_obj_name_dict['ald'] = 'com_frame'

#ALD From:
#Kenichiro Koshiyama, Tetsuya Kodama, Takeru Yano, Shigeo Fujikawa,
# "Molecular dynamics simulation of structural changes of lipid bilayers
# induced by shock waves: Effects of incident angles", Biochimica et Biophysica
# Acta (BBA) - Biomembranes, Volume 1778, Issue 6, June 2008, Pages 1423-1428
class ALDProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Estimate the average lateral displacement of lipids.

        This protocol is identified by the analysis key: 'ald'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', includes all lipid types.

        References:
            1. Kenichiro Koshiyama, Tetsuya Kodama, Takeru Yano, Shigeo
                Fujikawa, "Molecular dynamics simulation of structural changes
                of lipid bilayers induced by shock waves: Effects of incident
                angles", Biochimica et Biophysica Acta (BBA) - Biomembranes,
                Volume 1778, Issue 6, June 2008, Pages 1423-1428
        """
        # required
        self._short_description = "Average lateral displacement."
        self._return_length = 2
        self.analysis_key = 'ald'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'


        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.L_stat = RunningStats()
        self.analysis_output = []
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        leaflet = self.settings['leaflet']
        group = self.settings['resname']
        if self._first_frame:
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            elif leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            self._indices = indices
        indices = self._indices
        n_com = len(indices)
        selcoords = np.zeros((n_com, 2))

        count = 0
        for i in indices:
            com_curr = ba_reps['com_frame'].lipidcom[i].com_unwrap[
                ba_settings['lateral']]
            selcoords[count] = com_curr[:]
            count += 1

        # initialize a numpy array to hold the msd for the selection
        msd = np.zeros(2)
        # initialize a running stats object to do the averaging over resids
        #drs_stat = RunningStats()

        ref_coords = np.zeros((n_com, 2))
        if self._first_frame:
            count = 0
            for i in indices:
                com_curr = \
                    ba_reps['first_com_frame'].lipidcom[i].com_unwrap[
                        ba_settings['lateral']]
                ref_coords[count] = com_curr[:]
                count += 1
            self._ref_coords = ref_coords[:]
            self._first_frame = False
        else:
            ref_coords = self._ref_coords
            # get the current com frame list
        tc = ba_reps['com_frame'].time
        dt = tc
        dr = selcoords - ref_coords
        m_dr = []
        for val in dr:
            m_dr.append(np.sqrt(np.dot(val, val)))
        L = np.abs(m_dr).sum()/len(m_dr)
        self.L_stat.push(L)
        # get the msd for the current selection
        Lcurr = self.L_stat.mean()
        msd[0] = dt
        msd[1] = Lcurr
        self.analysis_output.append(msd)
        return

    def reset(self):
        self.L_stat.reset()
        self.analysis_output = []
        self._first_frame = True
        return
# update the command_protocols dictionary
command_protocols['ald'] = ALDProtocol

# define a new analysis 'apl_box'
valid_analysis.append('ac')
analysis_obj_name_dict['ac'] = 'mda_frame'
#Yoshimichi Andoha, Susumu Okazakia, Ryuichi Ueokab, "Molecular dynamics
# study of lipid bilayers modeling the plasma membranes of normal murine
# thymocytes and leukemic GRSL cells", Biochimica et Biophysica Acta (BBA)
# - Biomembranes, Volume 1828, Issue 4, April 2013, Pages 1259-1270.
# https://doi.org/10.1016/j.bbamem.2013.01.005
# Note: Area Compressibility is a defined in the reference is the inverse of
# area compressibility modulus.
class AreaCompressibilityProtocol(AnalysisProtocol):
    def __init__(self, args):

        # required
        self._short_description = "Isothermal area compressibility."
        self._return_length = 4
        self.analysis_key = 'ac'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['temperature'] = 298.15
        self._valid_settings = self.settings.keys()

        # parse input arguments
        self._parse_args(args)

        #output filename for pickle dump of results
        self.save_file_name = self.analysis_id + ".pickle"


        # storage for output
        self.n_frames = 0
        self.area_run = RunningStats()
        self.area_fluctuation = RunningStats()
        self.analysis_output = []
        self.first_comp = True
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'temperature':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict


    def reset(self):
        self.n_frames = 0
        self.area_run.reset()
        self.area_fluctuation.reset()
        self.analysis_output = []
        self.first_comp = True
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        dimensions = ba_reps['current_mda_frame'].dimensions[0:3]
        area = dimensions.prod()/dimensions[ba_settings['norm']]
        self.area_run.push(area)
        area_mean = self.area_run.mean()
        fluctuation = (area - area_mean)**2
        self.area_fluctuation.push(fluctuation)
        X_T = self.area_run.variance() / (area_mean * scicon.k * self.settings['temperature'])
        #conversion factor for Angstrom^2/Joules to meter^2/Joule
        X_T*=10.0**(-23)
        time = ba_reps['current_mda_frame'].time
        self.analysis_output.append([time, X_T])
        self.n_frames += 1
        return

command_protocols['ac'] = AreaCompressibilityProtocol

# define a new analysis
valid_analysis.append('lop')
analysis_obj_name_dict['lop'] = 'mda_frame'
#Based on P-N vector-normal angle:
#Zheng Li, Richard M. Venable, Laura A. Rogers, Diana Murray,
# and Richard W. Pastor, "Molecular Dynamics Simulations of PIP2 and PIP3
# in Lipid Bilayers: Determination of Ring Orientation, and the Effects of
# Surface Roughness on a Poisson-Boltzmann Description", Biophys J. 2009 Jul 8;
# 97(1): 155-163.
# doi:  10.1016/j.bpj.2009.04.037
class LateralOrientationParameterProtocol(AnalysisProtocol):
    def __init__(self, args):
        # required
        self._short_description = "Lateral orientation parameter."
        self._return_length = 4
        self.analysis_key = 'lop'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['resname'] = 'first'
        self.settings['leaflet'] = 'upper'
        self.settings['ref_atom_1'] = 'P'
        self.settings['ref_atom_2'] = 'N'
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.running = RunningStats()
        self.first_comp = True
        self.analysis_output = []
        self.selection = None
        self.norm_vec = np.zeros(3)
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        norm = ba_settings['norm']
        if self.first_comp:
            if self.settings['resname'] == 'first':
                groups = ba_reps['leaflets'][self.settings['leaflet']].get_group_names()
                self.settings['resname'] = groups[0]
            com_frame = ba_reps['com_frame']
            #norm_vec = np.zeros(3)
            self.norm_vec[norm] = 1.0
            indices = ba_reps['leaflets'][self.settings['leaflet']].get_group_indices(self.settings['resname'])
            for index in indices:
                res_com = com_frame.lipidcom[index]
                resid = res_com.resid

                sel_string = "resname "+self.settings['resname']+" and resid "+str(resid)
                res_sel = ba_mda_data.mda_universe.select_atoms(sel_string)
                if index == indices[0]:
                    self.selection = res_sel
                else:
                    self.selection += res_sel
            self.first_comp = False
        residues = self.selection.residues
        cos_run = RunningStats()
        for residue in residues:
            atom_1 = eval("residue.{}".format(self.settings['ref_atom_1']))
            atom_2 = eval("residue.{}".format(self.settings['ref_atom_2']))
            #atom_1 = ast.literal_eval("residue.{}".format(self.settings['ref_atom_1']))
            #atom_2 = ast.literal_eval("residue.{}".format(self.settings['ref_atom_2']))
            atom_1_i = atom_1.index
            atom_2_i = atom_2.index
            atom_1_coord = ba_reps['current_mda_frame'].positions[atom_1_i]
            atom_2_coord = ba_reps['current_mda_frame'].positions[atom_2_i]
            diff = atom_2_coord - atom_1_coord

            dist = np.sqrt(np.dot(diff, diff))
            cos_t = np.dot(diff, self.norm_vec)/dist
            cos_run.push(cos_t)
        cos_t_avg = cos_run.mean()
        self.running.push(cos_t_avg)
        time = ba_reps['current_mda_frame'].time
        #self.running.push(ap
        cos_t = np.zeros(self._return_length)
        cos_t[0] = time
        cos_t[1] = cos_t_avg
        cos_t[2] = self.running.mean()
        cos_t[3] = self.running.deviation()
        self.analysis_output.append(cos_t)

        return

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        self.first_comp = True
        return

command_protocols['lop'] = LateralOrientationParameterProtocol

# define a new analysis
valid_analysis.append('loa')
analysis_obj_name_dict['loa'] = 'mda_frame'
#Based on P-N vector-normal angle:
#Zheng Li, Richard M. Venable, Laura A. Rogers, Diana Murray,
# and Richard W. Pastor, "Molecular Dynamics Simulations of PIP2 and PIP3
# in Lipid Bilayers: Determination of Ring Orientation, and the Effects of
# Surface Roughness on a Poisson-Boltzmann Description", Biophys J. 2009 Jul 8;
# 97(1): 155-163.
# doi:  10.1016/j.bpj.2009.04.037
class LateralOrientationAngleProtocol(AnalysisProtocol):
    def __init__(self, args):
        # required
        self._short_description = "Lateral orientation angle."
        self._return_length = 4
        self.analysis_key = 'lop'
        self.analysis_id = 'none'
        #define adjustable settings
        self.settings = dict()
        self.settings['resname'] = 'first'
        self.settings['leaflet'] = 'upper'
        self.settings['ref_atom_1'] = 'P'
        self.settings['ref_atom_2'] = 'N'
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.running = RunningStats()
        self.first_comp = True
        self.analysis_output = []
        self.selection = None
        self.norm_vec = np.zeros(3)
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        norm = ba_settings['norm']
        if self.first_comp:
            if self.settings['resname'] == 'first':
                groups = ba_reps['leaflets'][self.settings['leaflet']].get_group_names()
                self.settings['resname'] = groups[0]
            com_frame = ba_reps['com_frame']
            #norm_vec = np.zeros(3)
            self.norm_vec[norm] = 1.0
            indices = ba_reps['leaflets'][self.settings['leaflet']].get_group_indices(self.settings['resname'])
            for index in indices:
                res_com = com_frame.lipidcom[index]
                resid = res_com.resid

                sel_string = "resname "+self.settings['resname']+" and resid "+str(resid)
                res_sel = ba_mda_data.mda_universe.select_atoms(sel_string)
                if index == indices[0]:
                    self.selection = res_sel
                else:
                    self.selection += res_sel
            self.first_comp = False
        residues = self.selection.residues
        cos_run = RunningStats()
        for residue in residues:
            atom_1 = eval("residue."+self.settings['ref_atom_1'])
            atom_2 = eval("residue."+self.settings['ref_atom_2'])
            atom_1_i = atom_1.index
            atom_2_i = atom_2.index
            atom_1_coord = ba_reps['current_mda_frame'].positions[atom_1_i]
            atom_2_coord = ba_reps['current_mda_frame'].positions[atom_2_i]
            diff = atom_2_coord - atom_1_coord

            dist = np.sqrt(np.dot(diff, diff))
            cos_t = np.dot(diff, self.norm_vec)/dist
            angle = 90.0 - np.arccos(cos_t)*180.0/np.pi

            cos_run.push(angle)
        cos_t_avg = cos_run.mean()
        self.running.push(cos_t_avg)
        time = ba_reps['current_mda_frame'].time
        #self.running.push(ap
        cos_t = np.zeros(self._return_length)
        cos_t[0] = time
        cos_t[1] = cos_t_avg
        cos_t[2] = self.running.mean()
        cos_t[3] = self.running.deviation()
        self.analysis_output.append(cos_t)

        return

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        self.first_comp = True
        return

command_protocols['loa'] = LateralOrientationAngleProtocol

# define a new analysis 'flip_flop'
valid_analysis.append('flip_flop')
analysis_obj_name_dict['flip_flop'] = 'com_frame'


class FlipFlopProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Count any lipid flips flops between the leaflets.

        This protocol is identified by the analysis key: 'flip_flop'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            None

        References:
            1. Andrey A. Gurtovenko, and Ilpo Vattulainen, Molecular Mechanism
                for Lipid Flip-Flops, The Journal of Physical Chemistry B, 2007
                111 (48), 13554-13559, DOI: 10.1021/jp077094k
                http://pubs.acs.org/doi/abs/10.1021/jp077094k?journalCode=jpcbfk

            2. Nicolas Sapay, W. F. Drew Bennett, and D. Peter Tieleman,
                Molecular Simulations of Lipid Flip-Flop in the Presence of
                Model Transmembrane Helices, Biochemistry, 2010 49 (35),
                7665-7673, DOI: 10.1021/bi100878q
                http://pubs.acs.org/doi/abs/10.1021/bi100878q
        """
        # required
        self._short_description = "Count lipid flip flops."
        self._return_length = 2
        self.analysis_key = 'flip_flop'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._resnames = []
        self._reference_leaf = None
        self._reference_com_frame = None
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = {}
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        if self._first_frame:
            self._reference_leaf = ba_reps['leaflets']
            self._reference_com_frame = ba_reps['com_frame']
            for leaflet in ba_reps['leaflets'].keys():
                resnames = ba_reps['leaflets'][leaflet].get_group_names()
                for resname in resnames:
                    if resname not in self._resnames:
                        self._resnames.append(resname)
            for resname in self._resnames:

                self.analysis_output[resname] = {'count': 0, 'events':[]}
            self._first_frame = False
            return
        curr_leaflets = ba_reps['leaflets']
        for leaflet in curr_leaflets.keys():
            o_leaf = [leaf for leaf in curr_leaflets.keys() if leaf != leaflet]
            if len(o_leaf) == 1:
                o_leaf = o_leaf[0]
            else:
                o_leaf = None
            curr = set(curr_leaflets[leaflet].get_member_indices())
            ref = set(self._reference_leaf[leaflet].get_member_indices())
            # o_leaf -> leaflet
            forward_diff = curr.difference(ref)
            # leaflet -> o_leaf
            backward_diff = ref.difference(curr)
            if len(forward_diff) > 0:
                for index in forward_diff:
                    resname = ba_reps['com_frame'].lipidcom[index].type
                    resid = ba_reps['com_frame'].lipidcom[index].resid
                    time = ba_reps['com_frame'].time
                    frame = ba_reps['com_frame'].mdnumber
                    self.analysis_output[resname]['count']+=1
                    self.analysis_output[resname]['events'].append([frame, time, resid,
                                                            o_leaf, leaflet])
            if len(backward_diff) > 0:
                for index in backward_diff:
                    resname = ba_reps['com_frame'].lipidcom[index].type
                    resid = ba_reps['com_frame'].lipidcom[index].resid
                    time = ba_reps['com_frame'].time
                    frame = ba_reps['com_frame'].mdnumber
                    self.analysis_output[resname]['count']+=1
                    self.analysis_output[resname]['events'].append([frame, time, resid,
                                                            leaflet, o_leaf])
            if len(forward_diff) > 0 or len(backward_diff) > 0:
                self._reference_leaf = ba_reps['leaflets']
                self._reference_com_frame = ba_reps['com_frame']
            break

        return

    def save_data(self, path=None):
        """Dumps the outputs of this protocol to disc.

        Args:
            path (str, Optional): The string containing the path to the location
                that the analysis results should be dumped to on disc.
        """
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        with open(save_file, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        """Returns the analysis_output of this protocol. """
        return self.analysis_output

    def reset(self):
        """Resets the analysis by resetting the outputs and any necessary
        internal variables.
        """
        self.analysis_output = {}
        self._first_frame = True
        self._resnames = []
        self._reference_leaf = None
        self._reference_com_frame = None
        return
# update the command_protocols dictionary
command_protocols['flip_flop'] = FlipFlopProtocol

# define a new analysis 'lipid_length'
valid_analysis.append('lipid_length')
analysis_obj_name_dict['lipid_length'] = 'vector_frame'

class LipidLengthProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the lipids length using the defined lipid vector.

        The LipidLengthProtocol is used to compute the mean lipid length using
        the vector represetation of the specified lipids.

        This protocol is identified by the analysis key: 'lipid_length'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', averages over all lipid types.

        References:
            1. Anton O. Chugunov,  Pavel E. Volynsky, Nikolay A. Krylov,
                Ivan A. Boldyrev, and Roman G. Efremov,  Liquid but Durable:
                Molecular Dynamics Simulations Explain the Unique Properties of
                Archaeal-Like Membranes, Scientific Reports, 4:7462, 2014,
                doi:10.1038/srep07462
                (https://www.nature.com/articles/srep07462)

        """
        # required
        self._short_description = "Estimate of lipid length using the lipid vectors."
        self._return_length = 2
        self.analysis_key = 'lipid_length'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        leaflet = self.settings['leaflet']
        group = self.settings['resname']
        if self._first_frame:
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            elif leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            self._indices = indices
        indices = self._indices
        n_com = len(indices)
        lengths = np.zeros((n_com, 1))
        count = 0
        for i in indices:
            vec_curr = ba_reps['vector_frame'].lipidvec[i].vector
            lengths[count] = np.sqrt(np.dot(vec_curr, vec_curr))
            count += 1

            # get the current com frame list
        tc = ba_reps['vector_frame'].time
        # get the msd for the current selection
        mean_length = lengths.mean()
        lipid_length = np.array([tc, mean_length])
        self.analysis_output.append(lipid_length)
        return


# update the command_protocols dictionary
command_protocols['lipid_length'] = LipidLengthProtocol

# define a new analysis
valid_analysis.append('lipid_tilt')
analysis_obj_name_dict['lipid_tilt'] = 'vector_frame'

class LipidTiltProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the lipids tilt angles using the defined lipid vector.

        The LipidTiltProtocol is used to compute the mean lipid tilt using
        the vector represetation of the specified lipids in reference to
        a particular axis, typically the bilayer normal.

        This protocol is identified by the analysis key: 'lipid_tilt'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'upper'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', averages over all lipid types.
            style (str: 'angle', 'order'): Specify whether to compute the
                tilt angle ('angle') or the tilt angle order parameter ('order').
                Default: 'angle'
            ref_axis (str: 'x', 'y', or 'z'): Specify the reference axis that
                should be used to estimate the tilt. This is typically the
                axis along the bilayer normal. Default: 'z'

        References:
            1. Anton O. Chugunov,  Pavel E. Volynsky, Nikolay A. Krylov,
                Ivan A. Boldyrev, and Roman G. Efremov,  Liquid but Durable:
                Molecular Dynamics Simulations Explain the Unique Properties of
                Archaeal-Like Membranes, Scientific Reports, 4:7462, 2014,
                doi:10.1038/srep07462
                (https://www.nature.com/articles/srep07462)

        """
        # required
        self._short_description = "Estimate of lipid tilt using the lipid vectors."
        self._return_length = 2
        self.analysis_key = 'lipid_tilt'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'upper'
        self.settings['resname'] = 'all'
        self.settings['style'] = 'angle'
        self.settings['ref_axis'] = 'z'
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        self._ref_axis = np.array([0.0, 0.0, 1.0])
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        leaflet = self.settings['leaflet']
        group = self.settings['resname']
        if self._first_frame:
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            elif leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_group_indices(group)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(group)
            self._indices = indices
            if self.settings['ref_axis'] != 'z':
                if self.settings['ref_axis'] == 'x':
                    self._ref_axis = np.array([1.0, 0.0, 0.0])
                elif self.settings['ref_axis'] == 'y':
                    self._ref_axis = np.array([0.0, 1.0, 0.0])

        indices = self._indices
        n_com = len(indices)
        values = np.zeros((n_com, 1))
        count = 0
        for i in indices:
            vec_curr = ba_reps['vector_frame'].lipidvec[i].vector
            if self.settings['style'] == 'angle':
                vec_length = np.sqrt(np.dot(vec_curr, vec_curr))
                angle_rad = np.arccos(np.dot(vec_curr,
                                                    self._ref_axis)/vec_length)
                values[count] = 90.0 - angle_rad*180.0/np.pi
            elif self.settings['style'] == 'order':
                vec_length = np.sqrt(np.dot(vec_curr, vec_curr))
                values[count] = np.dot(vec_curr, self._ref_axis)/vec_length

            count += 1

            # get the current com frame list
        tc = ba_reps['vector_frame'].time
        # get the msd for the current selection
        mean_value = values.mean()
        lipid_tilt = np.array([tc, mean_value])
        self.analysis_output.append(lipid_tilt)
        return


# update the command_protocols dictionary
command_protocols['lipid_tilt'] = LipidTiltProtocol

# define a new analysis
valid_analysis.append('lipid_collinearity')
analysis_obj_name_dict['lipid_collinearity'] = 'vector_frame'

class LipidCollinearityProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the lipid-lipid collinearity angles.

        The LipidCollinearityProtocol is used to compute the mean lipid-lipid
        collinearity angle (or order parameter) using the vector represetation
        of the specified lipids.

        This protocol is identified by the analysis key: 'lipid_collinearity'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'upper'
            resname_1 (str): Specify the resname of the reference lipid type to
                include in this analysis. Default: 'first', the first lipid in
                the list pulled from the com_frame representation.
            resname_2 (str): Specify the resname of the target lipid type to
                include in this analysis. Default: 'first', the first lipid in
                the list pulled from the com_frame representation.
            style (str: 'angle', 'order'): Specify whether to compute the
                tilt angle ('angle') or the tilt angle order parameter ('order').
                Default: 'angle'

        References:
            1. Anton O. Chugunov,  Pavel E. Volynsky, Nikolay A. Krylov,
                Ivan A. Boldyrev, and Roman G. Efremov,  Liquid but Durable:
                Molecular Dynamics Simulations Explain the Unique Properties of
                Archaeal-Like Membranes, Scientific Reports, 4:7462, 2014,
                doi:10.1038/srep07462
                (https://www.nature.com/articles/srep07462)

        """
        # required
        self._short_description = "Estimate of lipid-lipid collinearity."
        self._return_length = 2
        self.analysis_key = 'lipid_collinearity'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'upper'
        self.settings['resname_1'] = 'first'
        self.settings['resname_2'] = 'first'
        self.settings['style'] = 'angle'
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._ref_coords = None
        self._indices_1 = []
        self._indices_2 = []
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        if self._first_frame:
            leaflet = self.settings['leaflet']
            print("leaflet: {}".format(leaflet))
            group_1 = self.settings['resname_1']
            group_2 = self.settings['resname_2']
            lipid_types = []
            nlipids = 0
            do_leaflet = [self.settings['leaflet']]
            for leaflet_name in do_leaflet:
                leaf = ba_reps['leaflets'][leaflet_name]
                groups = leaf.get_group_names()
                nlipids += len(leaf.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.settings['resname_1'] == 'first':
                self.settings['resname_1'] =  lipid_types[0]
            if self.settings['resname_2'] == 'first':
                self.settings['resname_2'] = lipid_types[0]
            indices_1 = []
            indices_2 = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices_1 += curr_leaf.get_group_indices(group_1)
                    indices_2 += curr_leaf.get_group_indices(group_2)
            elif leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices_1 = curr_leaf.get_group_indices(group_1)
                indices_2 = curr_leaf.get_group_indices(group_2)
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices_1 = curr_leaf.get_group_indices(group_1)
                indices_2 = curr_leaf.get_group_indices(group_2)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices_1 += curr_leaf.get_group_indices(group_1)
                    indices_2 += curr_leaf.get_group_indices(group_2)
            self._indices_1 = indices_1
            self._indices_2 = indices_2
            self._first_frame = False

        indices_1 = self._indices_1
        indices_2 = self._indices_2
        n_com_1 = len(indices_1)
        n_com_2 = len(indices_2)
        values = []
        count = 0
        if self.settings['resname_1'] == self.settings['resname_2']:
            for i in range(n_com_1-1):
                ii = indices_1[i]
                vec_i = ba_reps['vector_frame'].lipidvec[ii].vector
                vec_i_l = np.sqrt(np.dot(vec_i, vec_i))
                for j in range(i+1, n_com_1):
                    jj = indices_1[j]
                    vec_j = ba_reps['vector_frame'].lipidvec[jj].vector
                    vec_j_l = np.sqrt(np.dot(vec_j, vec_j))
                    if self.settings['style'] == 'angle':
                        angle_rad = np.arccos(np.dot(vec_i, vec_j)/(vec_i_l*vec_j_l))
                        values.append(angle_rad*180.0/np.pi)
                    elif self.settings['style'] == 'order':
                        values.append( np.dot(vec_i, vec_j)/(vec_i_l*vec_j_l))
        else:
            for i in range(n_com_1):
                ii = indices_1[i]
                vec_i = ba_reps['vector_frame'].lipidvec[ii].vector
                vec_i_l = np.sqrt(np.dot(vec_i, vec_i))
                for j in range(n_com_2):
                    jj = indices_2[j]
                    vec_j = ba_reps['vector_frame'].lipidvec[jj].vector
                    vec_j_l = np.sqrt(np.dot(vec_j, vec_j))
                    if self.settings['style'] == 'angle':
                        angle_rad = np.arccos(np.dot(vec_i, vec_j)/(vec_i_l*vec_j_l))
                        values.append(angle_rad*180.0/np.pi)
                    elif self.settings['style'] == 'order':
                        values.append( np.dot(vec_i, vec_j)/(vec_i_l*vec_j_l))

        values = np.array(values)
            # get the current com frame list
        tc = ba_reps['vector_frame'].time
        # get the msd for the current selection
        mean_value = values.mean()
        lipid_collinearity = np.array([tc, mean_value])
        self.analysis_output.append(lipid_collinearity)
        return


# update the command_protocols dictionary
command_protocols['lipid_collinearity'] = LipidCollinearityProtocol


# define a new analysis 'msd'
valid_analysis.append('halperin_nelson')
analysis_obj_name_dict['halperin_nelson'] = 'com_frame'


class HalperinNelsonProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the mean Halperin and Nelson's rotational invariant.

        The HalperinNelsonProtocol is used to compute the mean Halperin and
        Nelson's rotational invariant. The value for lipid l is given by:
            phi_l = | (1/6) * sum_{j element nn(l)} exp(6i*theta_{lj}) |^2
        where i is complex and nn(l) are the 6 nearest neighbors of lipid
        l; theta_{lj} is the angle between the vector formed by beads
        representing lipid l and j and an arbitrary axis. The value is unity
        for perfect hexagonal packing, and it is zero to the extent that
        hexagonal packing is entirely absent. This protocol uses the 'com_frame'
        representation of the bilayer.

        This protocol is identified by the analysis key: 'halperin_nelson'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'upper'

        References:
            1. Shachi Katira, Kranthi K. Mandadapu, Suriyanarayanan
                Vaikuntanathan, Berend Smit, and David Chandler, The
                order-disorder transition in model lipid bilayers is a
                first-order hexatic to liquid phase transition, arXiv preprint
                [cond-mat.soft] 2015, arXiv:1506.04310.
                https://arxiv.org/pdf/1506.04310.pdf


        """
        # required
        self._short_description = "Halperin and Nelson's rotational invariant."
        self._return_length = 2
        self.analysis_key = 'halperin_nelson'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['leaflet'] = 'upper'
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._indices = []
        self._ref_axis = np.array([1.0, 0.0])
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return

    def _distance_euclidean_pbc(self, v_a, v_b, box_lengths, center='box_half'):
        if isinstance(center, (str, basestring)):
            if center == 'zero':
                center = np.zeros(len(v_a))
            elif center == 'box_half':
                center = box_lengths/2.0
        #shift center to zero for minimum image
        v_a = v_a - center
        v_b = v_b - center
        #difference
        d_v = v_a - v_b
        d_v_a = np.absolute(d_v)
        dim = len(v_a)
        #check for minimum image
        for i in range(dim):
            v_i = d_v_a[i]
            box_i = box_lengths[i]
            box_i_h = box_i/2.0
            if v_i > box_i_h:
                d_v[i] = box_i - np.absolute(v_a[i]) - np.absolute(v_b[i])
        r = np.sqrt(np.dot(d_v, d_v))
        return (r, d_v)

    def _6nn(self, indices, com_frame, lateral):
        k=6
        #initialize knn dict
        knn = {key: [] for key in indices}
        #make sure X has the right shape for the cdist function
        nX = len(indices)
        distances = [[indices[i], indices[j],\
         self._distance_euclidean_pbc(com_frame.lipidcom[indices[i]].com[lateral],\
         com_frame.lipidcom[indices[j]].com[lateral], com_frame.box[lateral], center='box_half')] \
         for i in xrange(nX-1) for j in xrange(i+1,nX)]
        #sort distances
        distances.sort(key=lambda x: x[2][0])
        #pick up the k nearest
        for d in distances:
            i = d[0]
            j = d[1]
            dist = d[2][0]
            vec = d[2][1]
            if len(knn[i]) < k:
                knn[i].append([j, dist, -vec])
            if len(knn[j]) < k:
                knn[j].append([i, dist, vec])
        return knn

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        leaflet = self.settings['leaflet']
        if self._first_frame:
            #determine the frame intervals
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "upper":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_member_indices()
            elif leaflet == "lower":
                curr_leaf = ba_reps['leaflets'][leaflet]
                indices = curr_leaf.get_member_indices()
            else:
                curr_leaf = ba_reps['leaflets']['upper']
                indices = curr_leaf.get_member_indices()
            self._indices = indices
            #self.analysis_output.append([])
        indices = self._indices
        n_com = len(indices)

        # initialize a running stats object to do the averaging over resids

        if self._first_frame:

            self._first_frame = False
        nn = self._6nn(indices, ba_reps['com_frame'],
                                    ba_settings['lateral'])
        hn = []
        for key in nn.keys():
            nn_curr = nn[key]
            nn_sum = 0.0
            for n in nn_curr:
                n_ind = n[0]
                vec = n[2]
                vec_l = n[1]
                angle = np.dot(vec, self._ref_axis)/vec_l
                nn_sum+=cmath.exp(6j*angle)
            hn.append(abs(nn_sum/6.0)**2)
        #for val in hn:
        #    print("hn: {}".format(val))
        hn_mean = np.array(hn).mean()
        self.analysis_output.append([ba_reps['com_frame'].time, hn_mean])
        return

    def reset(self):
        self.analysis_output = []
        self._first_frame = True
        self._indices = []
        return

# update the command_protocols dictionary
command_protocols['halperin_nelson'] = HalperinNelsonProtocol

# define a new analysis protocol
valid_analysis.append('area_fluctuation')
analysis_obj_name_dict['area_fluctuation'] = 'com_frame'


class AreaFluctuationProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the area fluctuation in the box along the bilayer laterals.

        This protocol is identified by the analysis key: 'area_fluctuation'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            None


        """
        # required
        self._short_description = "Bilayer lateral box area fluctuation."
        self._return_length = 2
        self.analysis_key = 'area_fluctuation'
        self.analysis_id = 'none'
        # default function settings
        # adjustable
        self.settings = dict()
        self.settings['none'] = None
        self._valid_settings = self.settings.keys()
        #self.leaflet = 'both'
        #self.resname = 'all'

        self._first_frame = True
        self._area_run = RunningStats()
        self._area_sq_run = RunningStats()
        # parse input arguments if given
        self._parse_args(args)
        self.save_file_name = self.analysis_id + ".pickle"
        # storage for output
        self.analysis_output = []
        return


    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        area = ba_reps['com_frame'].box[ba_settings['lateral']].prod()
        self._area_run.push(area)
        self._area_sq_run.push(area*area)

        if self._first_frame:
            self.analysis_output.append([ba_reps['com_frame'].time, area, 0.0])
            self._first_frame = False
        else:
            area_fluctuation = np.sqrt(self._area_sq_run.mean() - self._area_run.mean()**2)
            self.analysis_output.append([ba_reps['com_frame'].time, area, area_fluctuation])
        # print([ba_reps['com_frame'].time, area, area_fluctuation])
        return

    def reset(self):
        self._first_frame = True
        self.analysis_output = []
        self._area_run.reset()
        self._area_sq_run.reset()
        return

# update the command_protocols dictionary
command_protocols['area_fluctuation'] = AreaFluctuationProtocol

# define a new analysis 'disp_vec'
valid_analysis.append('disp_vec_corr_avg')
analysis_obj_name_dict['disp_vec_corr_avg'] = 'com_frame'


class DispVecCorrelationAverageProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Comute the pair-wise cross correlation between pairs of the
        displacement vectors for each lipid in the specified leaflet(s) of
        bilayer and do a inverse-distance weighted averaging.

        This analysis computes the displacement vectors as in DispVecProtocol,
        but then continues to compute the pair-wise cross correlations between
        each vector (i.e. the cos(theta) for the angle theta between the
        vectors) and averages the values using the inverse of the distance
        between the vector starting points as a weight.

        This protocol is identified by the analysis key: 'disp_vec_corr_avg'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname (str): Specify the resname of the lipid type to include in
                this analysis. Default: 'all', includes all lipid types.
            wrapped (bool): Specify whether to use the wrapped ('True') or
                un-wrapped ('False') coordintes for the base of the vectors.
                Default: False
            interval (int): Sets the frame interval over which to compute the
                    displacement vectors. f

        References:
            None
        """
        # required
        self._short_description = "Weighted average of displacement vector correlations."

        self._return_length = 4
        self.analysis_key = 'disp_vec_corr_avg'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self.settings['wrapped'] = False
        self.settings['interval'] = 5
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self._running = RunningStats()
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    arg_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0
        self._running.reset()

        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        if self.first_comp:
            self.last_com_frame = ba_reps['com_frame']
            self.first_comp = False
            return
        current_frame = ba_reps['current_mda_frame'].frame

        interval = (current_frame) - (self.last_frame)
        #print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = ba_reps['leaflets']['upper']
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = ba_reps['leaflets']['lower']
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                warnings.warn(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in ba_reps['leaflets']:
                    curr_leaf = ba_reps['leaflets'][leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            n_com = len(indices)
            x_index = ba_settings['lateral'][0]
            y_index = ba_settings['lateral'][1]
            box = ba_reps['current_mda_frame'].dimensions[0:3]
            box_x = box[x_index]
            box_y = box[y_index]
            box_l = np.array([box_x, box_y])
            box_x_h = box_x / 2.0
            box_y_h = box_y / 2.0
            # get the current frame
            curr_frame = ba_reps['com_frame']
            prev_frame = self.last_com_frame
            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com, 4))
            # vec_ends = []
            count = 0
            resnames = []
            for i in indices:
                resname = curr_frame.lipidcom[i].type
                resnames.append(resname)
                com_i = curr_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j_w = prev_frame.lipidcom[i].com[ba_settings['lateral']]
                if self.settings['wrapped']:
                    vec_ends[count, 0] = com_j_w[0]
                    vec_ends[count, 1] = com_j_w[1]
                else:
                    vec_ends[count, 0] = com_j[0]
                    vec_ends[count, 1] = com_j[1]
                vec_ends[count, 2] = com_i[0] - com_j[0]
                vec_ends[count, 3] = com_i[1] - com_j[1]
                #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count += 1
            total = 0.0
            weights = 0.0
            for i in range(n_com-1):
                index_i = indices[i]
                vec_end_a = vec_ends[i]
                vec_a = vec_end_a[2:4] - vec_end_a[0:2]
                com_i_w = prev_frame.lipidcom[index_i].com[ba_settings['lateral']]
                for j in range(i+1, n_com):
                    index_j = indices[j]
                    vec_end_b = vec_ends[j]
                    vec_b = vec_end_b[2:4] - vec_end_b[0:2]
                    dot = np.dot(vec_a, vec_b)
                    cos_t = dot/(np.linalg.norm(vec_a)*np.linalg.norm(vec_b))
                    com_j_w = prev_frame.lipidcom[index_j].com[ba_settings['lateral']]
                    pos_a = com_i_w[[x_index, y_index]]
                    pos_b = com_j_w[[x_index, y_index]]
                    dist = dist = distance_euclidean_pbc(pos_a, pos_b, l_box,
                                                         center='box_half')
                    total += cos_t * (1.0/dist)
                    weights += (1.0/dist)
            w_avg = total/weights
            self._running.push(w_avg)
            self.analysis_output.append([ba_reps['com_frame'].time,
                                         w_avg, self._running.mean(),
                                         self._running.deviation()])
            self.last_com_frame = ba_reps['com_frame']
            self.last_frame = current_frame
            #return vec_ends
        return


command_protocols['disp_vec_corr_avg'] = DispVecCorrelationAverageProtocol

# define a new analysis
valid_analysis.append('com_lateral_rdf')
analysis_obj_name_dict['com_lateral_rdf'] = 'com_frame'

class COMLateralRDFProtocol(AnalysisProtocol):

    def __init__(self, args):
        """Estimate the pair-wise radial distribution function in the bilayer
        lateral plane using the lipid centers of mass.

        This analysis protocol uses the 'com_frame' representation.

        This protocol is identified by the analysis key: 'com_lateral_rdf'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname_1 (str): Specify the resname of the reference lipid type to
                include in this analysis. Default: 'first', the first lipid in
                the list pulled from the com_frame representation.
            resname_2 (str): Specify the resname of the target lipid type to
                include in this analysis. Default: 'first', the first lipid in
                the list pulled from the com_frame representation.
            n_bins (int): Specifies the number of bins to use when estimating
                the RDF. Default: 25
            range_inner (float): Specify the inner distance cutoff for the RDF.
                Default: 0.0
            range_outer (float): Specify the outer distance cutoff for the RDF.
                Default: 25.0


        References:
            1.
        """
        # required
        self._short_description = "Lipid-lipid RDF in the bilayer lateral plane."

        self._return_length = 2
        self.analysis_key = 'com_lateral_rdf'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname_1'] = 'first'
        self.settings['resname_2'] = 'first'
        self.settings['n_bins'] = 25
        self.settings['range_inner'] = 0.0
        self.settings['range_outer'] = 25.0
        self._valid_settings = self.settings.keys()
        #parse input arguments/settings
        self._parse_args(args)

        self.save_file_name = self.analysis_id + ".pickle"

        # for outputs
        self._first_frame = True
        self._ref_coords = None
        self._indices = []
        self._count = None
        self._edges = None
        self._bins = None
        self._n_frames = 0
        self._area_run = RunningStats()
        self._N = 0
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'n_bins':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'range_inner':
                    arg_dict[arg_key] = float(arg_arg)
                elif arg_key == 'range_outer':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def reset(self):
        self.analysis_output = []
        self._first_frame = True
        self._count = None
        self._edges = None
        self._bins = None
        self._n_frames = 0
        self._N = 0
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):
        do_leaflet = []
        if self.settings['leaflet'] == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.settings['leaflet']]

        if self._first_frame:
            # initialize the RDF histogram
            rdf_range = [self.settings['range_inner'],
                         self.settings['range_outer']]
            count, edges = np.histogram([-1], bins=self.settings['n_bins'],
                                        range=rdf_range)
            count = count.astype(np.float64)
            count *= 0.0
            self._count = count
            self._edges = edges
            self._bins = 0.5 * (edges[:-1] + edges[1:])
            # build group/resname/lipid type list
            lipid_types = []
            nlipids = 0
            for leaflet_name in do_leaflet:
                leaflet = ba_reps['leaflets'][leaflet_name]
                groups = leaflet.get_group_names()
                nlipids += len(leaflet.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.settings['resname_1'] == 'first':
                self.settings['resname_1'] =  lipid_types[0]
            if self.settings['resname_2'] == 'first':
                self.settings['resname_2'] = lipid_types[0]

            n_ltypes = len(lipid_types)
            self.lipid_types = lipid_types
            self.n_ltypes = n_ltypes
            n_leaf = 0.0
            for leaflet_name in do_leaflet:
                self._N += len(ba_reps['leaflets'][leaflet_name].get_member_indices())
                n_leaf +=1
            self._N /= n_leaf
        lipid_types = self.lipid_types
        n_ltypes = self.n_ltypes
        #print lipid_types
        x_index = ba_settings['lateral'][0]
        y_index = ba_settings['lateral'][1]
        box = ba_reps['current_mda_frame'].dimensions[0:3]
        box_x = box[x_index]
        box_y = box[y_index]
        l_box = np.array([box_x, box_y])
        box_x_h = box_x / 2.0
        box_y_h = box_y / 2.0
        center = np.zeros(3)
        center[x_index] = box_x_h
        center[y_index] = box_y_h
        ltype_a = self.settings['resname_1']
        ltype_b = self.settings['resname_2']
        dists = []
        # COG = ba_reps['com_frame'].cog()
        # print "COG: ",COG," box_x_h: ",box_x_h
        for leaflet_name in do_leaflet:
            leaflet = ba_reps['leaflets'][leaflet_name]
            if leaflet.has_group(ltype_a) and leaflet.has_group(ltype_b):
                ltype_a_indices = leaflet.get_group_indices(ltype_a)
                ltype_b_indices = leaflet.get_group_indices(ltype_b)
                for i in ltype_a_indices:
                    for j in ltype_b_indices:
                        if i != j:
                            pos_a = ba_reps['com_frame'].lipidcom[i].com[[x_index, y_index]]
                            # print "pos_a: ",pos_a+COG
                            # print "pos_a - COG:",pos_a
                            pos_b = ba_reps['com_frame'].lipidcom[j].com[[x_index, y_index]]

                            dist = distance_euclidean_pbc(pos_a, pos_b,
                                                          l_box,
                                                          center='box_half')
                            ltype = ba_reps['com_frame'].lipidcom[j].type
                            #print "ltype: ",ltype," dist ",dist
                            dists.append(dist)

        rdf_range = [self.settings['range_inner'],
                     self.settings['range_outer']]
        count, bins = np.histogram(dists, bins=self.settings['n_bins'],
                                    range=rdf_range)
        self._count += count.astype(np.float64)
        self._area_run.push(ba_reps['com_frame'].box[ba_settings['lateral']].prod())
        self._n_frames += 1
        return

    def save_data(self, path=None):
        """Dumps the outputs of this protocol to disc.

        Args:
            path (str, Optional): The string containing the path to the location
                that the analysis results should be dumped to on disc.
        """
        # Area in each radial shell
        area = np.power(self._edges[1:], 2) - np.power(self._edges[:-1], 2)
        area *=  np.pi

        # Average number density
        box_area = self._area_run.mean()
        density = self._N / box_area

        rdf = self._count / (density * area * self._n_frames)

        output = (rdf, self._bins)
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        with open(save_file, 'wb') as outfile:
            pickle.dump(output, outfile)

        return

    def get_data(self):
        """Returns the analysis_output of this protocol. """
        # Area in each radial shell
        area = np.power(self._edges[1:], 2) - np.power(self._edges[:-1], 2)
        area *=  np.pi

        # Average number density
        box_area = self._area_run.mean()
        density = self._N / box_area

        rdf = self._count / (density * area * self._n_frames)

        return rdf, self._bins



# update the command_protocols dictionary
command_protocols['com_lateral_rdf'] = COMLateralRDFProtocol

# define a new analysis
valid_analysis.append('spatial_velocity_corr')
analysis_obj_name_dict['spatial_velocity_corr'] = 'com_frame'


class SpatialVelocityCorrelationFunctionProtocol(AnalysisProtocol):
    def __init__(self, args):
        """Comute the pair-wise cross correlation between pairs of the
        displacement vectors for each lipid in the specified leaflet(s) of
        bilayer and do a inverse-distance weighted averaging.

        This analysis computes the displacement vectors as in DispVecProtocol,
        but then continues to compute the pair-wise cross correlations between
        each vector (i.e. the cos(theta) for the angle theta between the
        vectors) and averages the values using the inverse of the distance
        between the vector starting points as a weight.

        This protocol is identified by the analysis key: 'disp_vec_corr_avg'

        Args:
            args (list): list of string keys and arguments

        Settings (parsed from args to settings dict):
            leaflet (str: 'both', 'upper', or 'lower'): Specifies the bilayer
                leaflet to include in the estimate. Default: 'both'
            resname_1 (str): Specify the resname of the reference lipid type to
                include in this analysis. Special names are 'first' and 'all',
                which use the first and all lipid types respectively. Default:
                'first', the first lipid in the list pulled from the com_frame
                representation.
            resname_2 (str): Specify the resname of the target lipid type to
                include in this analysis. Special names are 'first' and 'all',
                which use the first and all lipid types respectively. Default:
                'first', the first lipid in the list pulled from the com_frame
                representation.
            n_bins (int): Specifies the number of bins to use when estimating
                the RDF. Default: 25
            range_inner (float): Specify the inner distance cutoff for the RDF.
                Default: 0.0
            range_outer (float): Specify the outer distance cutoff for the RDF.
                Default: 25.0
            interval (int): Sets the frame interval over which to compute the
                    displacement vectors. f

        References:
            None
        Notes:
            The radial distance is centered on lipids of type resname_1 and
            averaging is taken over the pair-wise interactions of lipids of
            type resname_1 with lipids of type resname_2.
        """
        # required
        self._short_description = "Weighted average of displacement vector correlations."

        self._return_length = 4
        self.analysis_key = 'spatial_velocity_corr'
        self.analysis_id = 'none'

        # default function settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname_1'] = 'first'
        self.settings['resname_2'] = 'first'
        self.settings['n_bins'] = 25
        self.settings['range_inner'] = 0.0
        self.settings['range_outer'] = 25.0
        self.settings['interval'] = 5
        self._valid_settings = self.settings.keys()
        # parse input arguments if given
        self._parse_args(args)

        self.save_file_name = self.analysis_id + ".pickle"

        # storage for output
        self._costheta = []
        self._distances = []
        self.analysis_output = []
        self._first_frame = True
        #self._running = RunningStats()
        self.last_com_frame = None
        self.last_frame = 0
        self._processed_output = None

        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, arg_dict):

        for arg_key in arg_dict.keys():
            arg_arg = arg_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                if arg_key == 'n_bins':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'range_inner':
                    arg_dict[arg_key] = float(arg_arg)
                elif arg_key == 'range_outer':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                warnings.warn(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return arg_dict

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0
        self._costheta = []
        self._distances = []
        #self._running.reset()
        self._processed_output = None
        return

    def run_analysis(self, ba_settings, ba_reps, ba_mda_data):

        do_leaflet = []
        if self.settings['leaflet'] == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.settings['leaflet']]

        if self._first_frame:
            self.last_com_frame = ba_reps['com_frame']

            # build group/resname/lipid type list
            lipid_types = []
            nlipids = 0
            for leaflet_name in do_leaflet:
                leaflet = ba_reps['leaflets'][leaflet_name]
                groups = leaflet.get_group_names()
                nlipids += len(leaflet.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.settings['resname_1'] == 'first':
                self.settings['resname_1'] =  lipid_types[0]
            if self.settings['resname_2'] == 'first':
                self.settings['resname_2'] = lipid_types[0]
            self._first_frame = False
        current_frame = ba_reps['current_mda_frame'].frame

        interval = (current_frame) - (self.last_frame)
        #print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            for leaf in do_leaflet:
                curr_leaf = ba_reps['leaflets'][leaf]
                indices += curr_leaf.get_group_indices(self.settings['resname_1'])
                if self.settings['resname_1'] != self.settings['resname_2']:
                    indices += curr_leaf.get_group_indices(self.settings['resname_2'])
            n_com = len(indices)
            x_index = ba_settings['lateral'][0]
            y_index = ba_settings['lateral'][1]
            box = ba_reps['current_mda_frame'].dimensions[0:3]
            box_x = box[x_index]
            box_y = box[y_index]
            l_box = np.array([box_x, box_y])
            box_x_h = box_x / 2.0
            box_y_h = box_y / 2.0
            # get the current frame
            curr_frame = ba_reps['com_frame']
            prev_frame = self.last_com_frame
            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com, 4))
            # vec_ends = []
            count = 0
            resnames = []
            for i in indices:
                resname = curr_frame.lipidcom[i].type
                resnames.append(resname)
                com_i = curr_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    ba_settings['lateral']]
                com_j_w = prev_frame.lipidcom[i].com[ba_settings['lateral']]

                vec_ends[count, 0] = com_j[0]
                vec_ends[count, 1] = com_j[1]
                vec_ends[count, 2] = com_i[0] - com_j[0]
                vec_ends[count, 3] = com_i[1] - com_j[1]
                #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count += 1
            total = 0.0
            weights = 0.0
            for i in range(n_com-1):
                index_i = indices[i]
                vec_end_a = vec_ends[i]
                vec_a = vec_end_a[2:4] - vec_end_a[0:2]
                com_i_w = prev_frame.lipidcom[index_i].com[ba_settings['lateral']]
                ltype_i = prev_frame.lipidcom[index_i].type
                if (ltype_i == self.settings['resname_1']) or (self.settings['resname_1'] == 'all'):
                    for j in range(i+1, n_com):
                        index_j = indices[j]
                        ltype_j = prev_frame.lipidcom[index_j].type
                        if (ltype_j == self.settings['resname_2']) or (self.settings['resname_2'] == 'all'):
                            if index_i != index_j:
                                vec_end_b = vec_ends[j]
                                vec_b = vec_end_b[2:4] - vec_end_b[0:2]
                                dot = np.dot(vec_a, vec_b)
                                cos_t = dot/(np.linalg.norm(vec_a)*np.linalg.norm(vec_b))
                                com_j_w = prev_frame.lipidcom[index_j].com[ba_settings['lateral']]

                                dist = distance_euclidean_pbc(com_i_w[[x_index, y_index]],
                                                              com_j_w[[x_index, y_index]],
                                                              l_box, center='box_half')
                                self._costheta.append(cos_t)
                                self._distances.append(dist)

            self.last_com_frame = ba_reps['com_frame']
            self.last_frame = current_frame
            #return vec_ends
        return

    def save_data(self, path=None):
        """Dumps the outputs of this protocol to disc.

        Args:
            path (str, Optional): The string containing the path to the location
                that the analysis results should be dumped to on disc.
        """
        output = self.get_data()
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        with open(save_file, 'wb') as outfile:
            pickle.dump(output, outfile)

        return

    def _process_output(self):
        if self._processed_output is None:
            cos_t = np.array(self._costheta)
            positions = np.array(self._distances)
            pos_range = [self.settings['range_inner'], self.settings['range_outer']]
            self._processed_output = binned_average(cos_t, positions,
                                                    n_bins=self.settings['n_bins'],
                                                    position_range=pos_range)
        return

    def get_data(self):
        """Returns the analysis_output of this protocol. """
        self._process_output()
        return self._processed_output


command_protocols['spatial_velocity_corr'] = SpatialVelocityCorrelationFunctionProtocol
