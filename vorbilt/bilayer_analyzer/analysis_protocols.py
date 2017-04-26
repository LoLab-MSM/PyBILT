"""Analysis Protocols

This is a support module that defines a set of classes used to contruct the 'analysis protocol' and the 'analysiss' used
in the anlysis implemented by the vorbilt.bilayer_analyzer.bilayer_anlayzer.BilayerAnalyzer class.
The AnalysisProtocol is the class used to organize and initial the the individual analysis functions/protocols. The
individual analysis functions/protocols are derived classes of AnalysisProtocol, and the availabel analysiss can be
extended by simply defining a new protocol.
Example:
    # define a new analysis 'my_analysis'
    # add it to the valid_analysislist
    valid_analysis.append('my_analysis')
    #add it to analysis_obj_name_dict dictionary with the buildable object (e.g. mda_frame)  needed for the analysis.
    #The buildables are built by the bilayer analyzer and can accessed from it for the analysisanalysis.
    analysis_obj_name_dict['my_analysis'] = 'mda_frame'
    class MyAnalysisProtocol(AnalysisFunctionProtocol):
        #minimal def the __init__ and run_analysisfuncions
        def __init__(self, args):
            #define the initialization
        #run analysis always takes the bilayer_analyzer class object as the input
        def run_analysis(self, bilayer_analyzer):
            #.............
            #doing my analysis
            #................
            return

        #redefine any other functions from AnalysisProtocol or add new ones.


"""

# imports
import scipy.constants as scicon

try:
    import cPickle as pickle
except:
    import pickle

# VORBILT imports
from vorbilt.common.running_stats import *
import vorbilt.mda_tools.mda_density_profile as mda_dp
import vorbilt.lipid_grid.lipid_grid_curv as lgc

#need some containers for bookkeeping
command_protocols = {}
valid_analysis= []
analysis_obj_name_dict = {}
# obj_dict = {"com_frame":COMFrame}
#define the buildable objects that are used by the analysisfunctions
use_objects = {"mda_frame": True, "com_frame": True, "lipid_grid": False}


# TO DO:
# have protocols parse dictionaries as well as strings, which is more Pythonic


# protocol for the analysiss to run during the frame loop
class Analyses:
    """A class to facilitate analysis of the bilayers
    This object stores all the analyses that being performed and provides fucntionality to add and remove
    analyses.

    Attributes:
        use_objects (list of str): A list of buildable objects that need to contructed in the BilayerAnalyzer for the
            analysiss defined in this protocol.
        in_commands (list): A list of the input strings for the analysiss to be used.
        arguments (list): A list of the arguments for analysiss.
        analysis_keys (list): A list of the keys assigned to analysiss.
        command_protocol (dict): A dictionary of the analysis objects.
        analysis_ids (list): A list of the ids assigned to analysiss.
        n_commands (int): The number of initialized analysiss.

    """
    def __init__(self, analysis_commands):
        self.use_objects = use_objects
        self.in_commands = analysis_commands
        # check analysiss
        arguments = []
        analysis_keys = []
        command_protocol = {}
        analysis_ids = []
        for command in analysis_commands:
            if len(command) < 2:
                raise RuntimeError(
                    "wrong number of arguments for analysis " + str(command))
            comp_key = command[0]
            comp_id = command[1]
            comp_args = command[1:]
            if (comp_key in valid_analysis):
                if (len(comp_args) >= 1):
                    if comp_id not in analysis_ids:
                        comp_object = analysis_obj_name_dict[comp_key]
                        self.use_objects[comp_object] = True
                        arguments.append(comp_args)
                        analysis_ids.append(comp_args[0])
                        analysis_keys.append(comp_key)
                        command_protocol[comp_args[0]] = command_protocols[
                            comp_key](comp_args)
                    else:
                        raise RuntimeError(
                            "analysisid \'" + comp_id + "\' has already been used")
                else:
                    raise RuntimeError(
                        "wrong number of arguments for analysis " + str(
                            command))
            else:
                raise RuntimeError(
                    "invalid analysis \'" + comp_key + "\' :" + str(command))
        self.arguments = arguments
        self.analysis_keys = analysis_keys
        self.analysis_ids = analysis_ids
        self.command_protocol = command_protocol
        self.n_commands = len(command_protocol)
        # object dependencies
        if self.use_objects['lipid_grid']:
            self.use_objects['com_frame'] = True
        return

    def add_analysis(self, analysis_string):
        command = analysis_string.split()
        comp_key = command[0]
        comp_id = command[1]
        comp_args = command[1:]
        if (comp_key in valid_analysis):
            if (len(comp_args) >= 1):
                if comp_id not in self.analysis_ids:
                    comp_object = analysis_obj_name_dict[comp_key]
                    self.use_objects[comp_object] = True
                    self.arguments.append(comp_args)
                    self.analysis_ids.append(comp_args[0])
                    self.analysis_keys.append(comp_key)
                    self.command_protocol[comp_args[0]] = command_protocols[
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

    def remove_analysis(self, analysis_id):
        if analysis_id in self.analysis_ids:
            del self.command_protocol[analysis_id]
            index = self.analysis_ids.index(analysis_id)
            del self.arguments[index]
            del self.analysis_keys[index]
            del self.analysis_ids[index]
            self.n_commands -= 1
        else:
            raise RuntimeWarning("no analysis with id '{}'".format(analysis_id))
        return

    def print_protocol(self):
        print ('build objects:')
        for key in self.use_objects.keys():
            if self.use_objects[key]:
                print (key)
        print ("with analysis:")
        for analysis_id in self.analysis_ids:
            self.command_protocol[analysis_id].print_protocol()
        return

    def dump_data(self):
        print ('dumping analysis data to pickle files...')
        for analysis_id in self.analysis_ids:
            print ("analysis id: {} ---> {} ".format(
                analysis_id,
                self.command_protocol[analysis_id].save_file_name))
            self.command_protocol[analysis_id].save_data()

    def reset(self):
        for analysis_id in self.analysis_ids:
            self.command_protocol[analysis_id].reset()
        return

#base class for analysis protocols
class AnalysisProtocol:
    '''Base class for analysis protocols.

    Args:
        args (list): list of argument keys and values.

    Attributes:
        valid_args (list): list of the valid arguemnt keys.
        return_length (int): if applicable, length of the return vector
        analysis_key (str): the key name of this analysis.
        analysis_id (str): the unique id assigned to this analyisis.
        save_file_name (str): the path and filename for the pickle file output of
            this analysis' results.
        analysis_output (list or list like): Used to store the ouptut of this
            analyis during the frame loop.
    '''
    def __init__(self, args):

        # required
        self.valid_args = ['none']
        self.return_length = 1
        self.analysis_key = 'none'
        self.analysis_id = ''
        self.settings = dict{}
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        self._parse_args(args)
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        '''Parses the setup arguments for this analysis
        Args:
            args (list): List of argument keys and values.
        '''
        # print args
        if isinstance(args, list):
            self._parse_list(args)
        elif isinstance(args, dict):
            self._parse_dict(args)
        return
        # required - a check protocol function which reports relevant settings
    def _parse_string(self, args):
        arg_list = analysis_string.split()
        nargs = len(arg_list)
        if len(arg_list)
        self._parse_list(arg_list)
        return

    def _parse_list(self, args):
        nargs = len(args)
        self.analysis_id = args[0]
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                self.settings[arg_key] =  arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_key)
        return

    def _parse_dict(self, args):
        if 'analysis_id' not in args.keys():
            raise RuntimeError("required key \'anlaysis_id\' not assigned in input dict for analysis type: \'"+self.analysis_key+"\'")
        for arg_key in args.keys():
            arg_arg = args[arg_key]
            if arg_key in self.valid_args:
                self.settings[arg_key] =  arg_arg
            elif arg_key == 'analysis_id':
                self.analysis_id = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_key)
        return
    def print_protocol(self):
        print ("Parent protocol class for analysiss.")
        return

    def run_analysis(self, bilayer_analyzer):
        # do some stuff
        # get an output
        output = np.zeros(self.return_length)
        # save the output
        self.analysis_output.append(output)
        return

    def save_data(self):
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(np.array(self.analysis_output), outfile)

        return

    def get_data(self):
        return np.array(self.analysis_output)

    def reset(self):
        self.analysis_output = []
        return


# define a new analysis 'msd'
valid_analysis.append('msd')
analysis_obj_name_dict['msd'] = 'com_frame'


class MSDProtocol(AnalysisProtocol):
    '''Mean squared displacement
    Args:
        args (list): list of string keys and arguments
    '''
    def __init__(self, args):

        # required
        self.valid_args = ['leaflet', 'resname']
        self.return_length = 2
        self.analysis_key = 'msd'
        self.analysis_id = args[0]
        # default function settings
        # adjustable
        self.settings = dict{}
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        #self.leaflet = 'both'
        #self.resname = 'all'

        self.save_file_name = self.analysis_id + ".pickle"
        self.first_frame = True
        self.ref_coords = None
        self.indices = []
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                if arg_key == 'leaflet':
                    self.settings['leaflet'] = arg_arg
                elif arg_key == 'resname':
                    self.settings['resname'] = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for analysis \'msd\'")
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
            "\'" + self.analysis_id + "\': msd analysisof " + self.settings['resname'] + " lipids in " + self.settings['leaflet'] + " leaflet(s).")
        return

    def run_analysis(self, bilayer_analyzer):
        leaflet = self.settings['leaflet']
        group = self.settings['resname']
        if self.first_frame:
            indices = []
            # parse the leaflet and group inputs
            if leaflet == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(group)
            elif leaflet == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(group)
            elif leaflet == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(group)
            else:
                # unknown option--use default "both"
                print "!! Warning - request for unknown leaflet name \'", leaflet
                print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(group)
            self.indices = indices
        indices = self.indices
        n_com = len(indices)
        selcoords = np.zeros((n_com, 2))

        count = 0
        for i in indices:
            com_curr = bilayer_analyzer.com_frame.lipidcom[i].com_unwrap[
                bilayer_analyzer.lateral]
            selcoords[count] = com_curr[:]
            count += 1

        # initialize a numpy array to hold the msd for the selection
        msd = np.zeros(2)
        # initialize a running stats object to do the averaging over resids
        drs_stat = RunningStats()

        ref_coords = np.zeros((n_com, 2))
        if self.first_frame:
            count = 0
            for i in indices:
                com_curr = \
                    bilayer_analyzer.first_com_frame.lipidcom[i].com_unwrap[
                        bilayer_analyzer.lateral]
                ref_coords[count] = com_curr[:]
                count += 1
            self.ref_coords = ref_coords[:]
            self.first_frame = False
        else:
            ref_coords = self.ref_coords
            # get the current com frame list
        tc = bilayer_analyzer.com_frame.time
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

# define a new analysis 'apl_box'
valid_analysis.append('apl_box')
analysis_obj_name_dict['apl_box'] = 'mda_frame'


class APLBoxProtocol(AnalysisProtocol):
    def __init__(self, args):
        # required
        self.valid_args = None
        self.return_length = 4
        self.analysis_key = 'apl_box'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            raise RuntimeWarning(
                "ignoring extra arguments passed to area per lipid analysis \'" + self.analysis_id + "\'")

        # storage for output
        self.running = RunningStats()
        self.analysis_output = []
        return

    # required- function to parse the input arguments
    @staticmethod
    def _parse_args():
        pass

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print ("'{}'': area per lipid using box "
               "dimensions".format(self.analysis_id))
        return

    def run_analysis(self, bilayer_analyzer):
        box = bilayer_analyzer.current_mda_frame.dimensions[0:3]
        plane = box[bilayer_analyzer.lateral]
        area = plane[0] * plane[1] * 2.0
        nlipids = bilayer_analyzer.mda_data.n_residues
        apl = area / nlipids
        time = bilayer_analyzer.current_mda_frame.time
        self.running.push(apl)
        apl_t = np.zeros(self.return_length)
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
        # required
        self.valid_args = None
        self.return_length = 4
        self.analysis_key = 'bilayer_thickness'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            warn = "ignoring extra arguments passed to " \
                   "bilayer thickness analysis '{}''".format(self.analysis_id)
            raise RuntimeWarning(warn)

        # storage for output
        self.running = RunningStats()
        self.analysis_output = []
        return

    # required- function to parse the input arguments
    def _parse_args():
        pass

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print (
            "\'" + self.analysis_id + "\': bilayer thickness using lipid_grid")
        return

    def run_analysis(self, bilayer_analyzer):
        current_thickness = bilayer_analyzer.lipid_grid.average_thickness()
        # print (current_thickness)
        time = bilayer_analyzer.current_mda_frame.time
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

        # required
        self.valid_args = None
        self.return_length = 4
        self.analysis_key = 'apl_grid'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)

        # storage for output
        self.running = RunningStats()
        self.analysis_output = {}
        self.first_comp = True
        self.running_res = {}
        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        raise RuntimeWarning(
            "ignoring extra arguments passed to area per lipid (by grid) analysis \'" + self.analysis_id + "\'")
        return

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print ("\'" + self.analysis_id + "\': area per lipid using lipid_grid")
        return

    def run_analysis(self, bilayer_analyzer):
        apl_grid_out = bilayer_analyzer.lipid_grid.area_per_lipid()
        # print (current_thickness)
        time = bilayer_analyzer.current_mda_frame.time
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

    def save_data(self):
        for key in self.analysis_output.keys():
            self.analysis_output[key] = np.array(self.analysis_output[key])
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        for key in self.analysis_output.keys():
            self.analysis_output[key] = np.array(self.analysis_output[key])
        return self.analysis_output

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        self.first_comp = True
        self.running_res = {}
        return


command_protocols['apl_grid'] = APLGridProtocol

# define a new analysis 'disp_vec'
valid_analysis.append('disp_vec')
analysis_obj_name_dict['disp_vec'] = 'com_frame'


class DispVecProtocol(AnalysisProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['interval', 'leaflet', 'resname', "wrapped"]
        self.return_length = 4
        self.analysis_key = 'disp_vec'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        self.leaflet = 'both'
        self.group = 'all'
        self.wrapped = False
        self.interval = 10
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)

        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            print ('arg_key: ', arg_key)
            if arg_key in self.valid_args:
                if arg_key == 'interval':
                    self.interval = int(arg_arg)
                elif arg_key == 'leaflet':
                    if arg_arg in ['upper', 'lower', 'both']:
                        self.leaflet = arg_arg
                elif arg_key == 'resname':
                    self.group = arg_arg
                elif arg_key == 'wrapped':
                    if arg_arg in ['True', 'true']:
                        self.wrapped = arg_arg
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("\'" + self.analysis_id + "\': displacement vectors for " + str(
            self.interval) + " frame intervals")
        return

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_analysis(self, bilayer_analyzer):

        if self.first_comp:
            self.last_com_frame = bilayer_analyzer.com_frame
            self.first_comp = False
            return
        current_frame = bilayer_analyzer.current_mda_frame.frame

        interval = (current_frame) - (self.last_frame)
        print (interval, " ", self.interval)
        if interval == self.interval:
            indices = []
            if self.leaflet == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.group)
            elif self.leaflet == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.group)

            elif self.leaflet == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.group)
            else:
                # unknown option--use default "both"
                raise RuntimeWarning(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.leaflet = 'both'
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.group)
            n_com = len(indices)

            # print "there are ",len(indices)," members"
            xi = bilayer_analyzer.lateral[0]
            yi = bilayer_analyzer.lateral[1]
            zi = bilayer_analyzer.norm

            vec_ends_out = []
            # get the current frame
            curr_frame = bilayer_analyzer.com_frame
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
                    bilayer_analyzer.lateral]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    bilayer_analyzer.lateral]
                com_j_w = prev_frame.lipidcom[i].com[bilayer_analyzer.lateral]
                if self.wrapped:
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
            self.last_com_frame = bilayer_analyzer.com_frame
            self.last_frame = current_frame
            return vec_ends
        return

    def save_data(self):

        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(self.analysis_output, outfile)

        return

    def get_data(self):
        return self.analysis_output


command_protocols['disp_vec'] = DispVecProtocol

# define a new analysis
valid_analysis.append('mass_dens')
analysis_obj_name_dict['mass_dens'] = 'mda_frame'


class MassDensProtocol(AnalysisProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['selection', 'n_bins']
        self.return_length = None
        self.analysis_key = 'mass_dens'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        self.selection_string = 'all'
        self.n_bins = 25
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)

        # storage for output
        self.centers = None
        self.n_frames = 0
        self.analysis_output = []
        self.first_comp = True
        self.selection = None

        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        read_sel_string = False
        n_bins_arg = False
        for i in range(1, nargs, 1):
            arg_key = args[i]

            print ('arg_key: ', arg_key)
            if arg_key in self.valid_args:
                if arg_key == 'n_bins':
                    arg_arg = args[i + 1]
                    self.n_bins = int(arg_arg)
                    read_sel_string = False
                    n_bins_arg = True
                elif arg_key == 'selection':
                    selection_words = [args[j] for j in range(i + 1, nargs) if
                                       (args[j] not in self.valid_args)]
                    i += len(selection_words)
                    selection_string = ""
                    for word in selection_words:
                        selection_string += " " + word
                    self.selection_string = selection_string
                    read_sel_string = True
                    n_bins_arg = False
            elif not (read_sel_string or n_bins_arg):
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
            "\'" + self.analysis_id + "\': mass density analysis for selection \'" + self.selection_string + "\' ")
        return

    def reset(self):
        self.centers = None
        self.n_frames = 0
        self.analysis_output = []
        self.first_comp = True
        self.selection = None

        return

    def run_analysis(self, bilayer_analyzer):
        first = self.first_comp
        if self.first_comp:
            # print self.selection_string
            if self.selection_string == ' BILAYER':
                self.selection = bilayer_analyzer.mda_data.bilayer_sel
            else:
                self.selection = bilayer_analyzer.mda_data.mda_universe.select_atoms(
                    self.selection_string)
            self.first_comp = False

        # print "there are ",len(indices)," members"

        norm_axis = bilayer_analyzer.normal_dimension
        ref_sel = bilayer_analyzer.mda_data.bilayer_sel
        # mda_density_profile
        centers_density = mda_dp.mass_density_profile(
            bilayer_analyzer.mda_data.mda_trajectory, self.selection,
            fstart=bilayer_analyzer.frame_index,
            fend=bilayer_analyzer.frame_index + 1, axis=norm_axis,
            nbins=self.n_bins, refsel=ref_sel)
        if first:
            self.centers = centers_density[0]
            self.analysis_output = np.zeros(len(self.centers))
        self.analysis_output += centers_density[1]
        self.n_frames += 1
        return

    def save_data(self):
        centers_density = (self.centers, self.analysis_output / self.n_frames)
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(centers_density, outfile)

        return

    def get_data(self):
        return (self.centers, self.analysis_output / self.n_frames)


command_protocols['mass_dens'] = MassDensProtocol


# define a new analysis
valid_analysis.append('acm')
analysis_obj_name_dict['acm'] = 'mda_frame'

#L. Janosi and A. A. Gorfe, J. Chem. Theory Comput. 2010, 6, 3267-3273
#D. Aguayo, F. D. Gonzalez-Nilo, and C. Chipot, J. Chem. Theory Comput. 2012, 8, 1765-1773
class AreaCompressibilityModulusProtocol(AnalysisProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['temperature']
        self.return_length = 4
        self.analysis_key = 'acm'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)

        # storage for output
        self.temperature = 298.15
        self.n_frames = 0
        self.area = []
        self.apl = []
        self.times = []
        self.analysis_output = []
        self.first_comp = True
        self.per_leaflet = 1
        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]

            # print('arg_key: ', arg_key)
            if arg_key in self.valid_args:
                if arg_key == 'temperature':
                    arg_arg = args[i + 1]
                    self.temperature = float(arg_arg)
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print("\'" + self.analysis_id + "\': area compressibility modulus analysis for selection \'")
        return

    def reset(self):
        self.area_run.reset()
        self.area = []
        self.apl = []
        self.n_frames = 0
        self.analysis_output = []
        self.first_comp = True
        return

    def run_analysis(self, bilayer_analyzer):

        if self.first_comp:
            nlipids = bilayer_analyzer.mda_data.n_residues
            per_leaflet = nlipids / 2
            self.per_leaflet = per_leaflet
            self.first_comp = False
        lateral_indices = bilayer_analyzer.lateral
        dimensions = bilayer_analyzer.current_mda_frame.dimensions[0:3]
        lateral_dim = dimensions[lateral_indices]
        area = lateral_dim.prod()
        print(area)
        self.area.append(area)
        apl = area / self.per_leaflet
        print(apl)
        self.apl.append(apl)
        time = bilayer_analyzer.current_mda_frame.time
        self.times.append(time)
        self.n_frames += 1
        return

    def save_data(self):
        data = self.get_data()
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(data, outfile)
        return

    def get_data(self):
        area_eq = np.array(self.area).mean()
        apl = np.array(self.apl)
        apl_run = gen_running_average(apl)
        apl_minus_area = (apl - area_eq)**2
        apl_minus_area_run = gen_running_average(apl_minus_area)
        acm = scicon.k * self.temperature * apl_run[:,0] / ( self.per_leaflet * apl_minus_area_run[:,0])
        #conversion factor for Joules/Angstron^2 to milliNewtons/meter
        acm*=10.0**23
        acm_run = gen_running_average(acm)
        times = np.array(self.times)
        npoints = len(times)
        output = np.zeros((npoints, 4))
        output[:, 0] = times[:]
        output[:, 1] = acm[:]
        output[:, 2] = acm_run[:, 0]
        output[:, 3] = acm_run[:, 1]
        return output
command_protocols['acm'] = AreaCompressibilityModulusProtocol

# define a new analysis 'nnf'
valid_analysis.append('nnf')
analysis_obj_name_dict['nnf'] = 'com_frame'

#Found in: M. Orsi and J. W. Essex, Faraday Discuss., 2013, 161, 249-272
#Originally described in: A. H. de Vries, A. E. Mark and S. J. Marrink, J. Phys. Chem. B, 2004, 108, 2454-2463
class NNFProtocol(AnalysisProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['leaflet', 'resname_1', 'resname_2', 'n_neighbors']
        self.return_length = 2
        self.analysis_key = 'nnf'
        self.analysis_id = args[0]
        # default function settings
        self.leaflet = 'both'
        self.n_neighbors = 5
        self.rname_1 = 'first'
        self.rname_2 = 'first'
        self.save_file_name = self.analysis_id + ".pickle"
        self.first_frame = True
        self.ref_coords = None
        self.indices = []
        self.running = RunningStats()
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)
        else:
            raise RuntimeError('wrong number of arguments given to nnf analysis. You must specifiy ')
        # storage for output
        self.analysis_output = []
        return

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        self.first_frame = True
        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        res1 = False
        res2 = False
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                if arg_key == 'leaflet':
                    self.leaflet = arg_arg
                elif arg_key == 'resname_1':
                    self.rname_1 = arg_arg
                    res1 = True
                elif arg_key == 'resname_2':
                    self.rname_2 = arg_arg
                    res2 = True
                elif arg_key == 'n_neighbors':
                    self.n_neighbors = int(arg_arg)
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for analysis \'nnf\'")
        #if res1 and res2:
        #    return
        #else:
        #   raise RuntimeError("must specify target types by arg_keys \'resname_1\' and \'resname_2\'")
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
            "\'" + self.analysis_id + "\': lateral order nearest neighbor fraction (nnf) analysis")
        return


    def run_analysis(self, bilayer_analyzer):
        do_leaflet = []
        if self.leaflet == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.leaflet]

        if self.first_frame:
            #pass
            # build group/resname/lipid type list
            lipid_types = []
            nlipids = 0
            for leaflet_name in do_leaflet:
                leaflet = bilayer_analyzer.leaflets[leaflet_name]
                groups = leaflet.get_group_names()
                nlipids += len(leaflet.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.rname_1 == 'first':
                self.rname_1 =  lipid_types[0]
            if self.rname_2 == 'first':
                self.rname_2 = lipid_types[0]
            if self.n_neighbors > nlipids:
                self.n_neighbors = nlipids-1

            n_ltypes = len(lipid_types)
            self.lipid_types = lipid_types
            self.n_ltypes = n_ltypes

        lipid_types = self.lipid_types
        n_ltypes = self.n_ltypes
        #print lipid_types
        x_index = bilayer_analyzer.lateral[0]
        y_index = bilayer_analyzer.lateral[1]
        box = bilayer_analyzer.current_mda_frame.dimensions[0:3]
        box_x = box[x_index]
        box_y = box[y_index]
        box_x_h = box_x / 2.0
        box_y_h = box_y / 2.0
        ltype_a = self.rname_1
        ltype_b = self.rname_2
        #print "ltype_a: ",ltype_a," ltype_b: ",ltype_b
        avg_frac = RunningStats()
        for leaflet_name in do_leaflet:
            leaflet = bilayer_analyzer.leaflets[leaflet_name]

            if leaflet.has_group(ltype_a) and leaflet.has_group(ltype_b):
                ltype_a_indices = leaflet.get_group_indices(ltype_a)
                ltype_b_indices = leaflet.get_group_indices(ltype_b)
                all_index = leaflet.get_member_indices()
                for i in ltype_a_indices:
                    neighbors = []
                    for j in all_index:
                        if i != j:
                            pos_a = bilayer_analyzer.com_frame.lipidcom[i].com
                            pos_b = bilayer_analyzer.com_frame.lipidcom[j].com
                            dx = np.abs(pos_a[x_index] - pos_b[x_index])
                            dy = np.abs(pos_a[y_index] - pos_b[y_index])
                            # minimum image for wrapped coordinates
                            if dx > box_x_h:
                                dx = box_x - np.absolute(pos_a[x_index] - box_x_h) - np.absolute(
                                    pos_b[x_index] - box_x_h)

                            if dy > box_y_h:
                                dy = box_y - np.absolute(pos_a[y_index] - box_y_h) - np.absolute(
                                    pos_b[y_index] - box_y_h)
                            dist = np.sqrt(dx ** 2 + dy ** 2)
                            ltype = bilayer_analyzer.com_frame.lipidcom[j].type
                            #print "ltype: ",ltype," dist ",dist
                            neighbors.append([j, dist, ltype])
                    neighbors.sort(key=lambda x: x[1])
                    nn_neighbors = neighbors[0:self.n_neighbors]
                    #print neighbors
                    #print nn_neighbors
                    #quit()
                    n_type_b = 0.0
                    for neighbor in nn_neighbors:
                        ntype = neighbor[2]
                        if ntype == ltype_b:
                            n_type_b += 1.0
                    frac = n_type_b/self.n_neighbors
                    #print "frac ",frac
                    avg_frac.push(frac)

        f_current = avg_frac.mean()
        self.running.push(f_current)
        f_run = self.running.mean()
        f_run_dev = self.running.deviation()
        tc = bilayer_analyzer.com_frame.time
        self.analysis_output.append([tc, f_current, f_run, f_run_dev])
        return

    # def run_analysis(self, bilayer_analyzer):
    #
    #     if self.first_frame:
    #         pass
    #         #build group/resname/lipid type list
    #         lipid_types = []
    #         for leaflet in bilayer_analyzer.leaflets:
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
    #     x_index = bilayer_analyzer.lateral[0]
    #     y_index = bilayer_analyzer.lateral[1]
    #     box = bilayer_analyzer.current_mda_frame.dimensions[0:3]
    #     box_x = box[x_index]
    #     box_y = box[y_index]
    #     box_x_h = box_x/2.0
    #     box_y_h = box_y/2.0
    #
    #     for leaflet in bilayer_analyzer.leaflets:
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
    #                             pos_a = bilayer_analyzer.com_frame.lipidcom[lindex].com
    #                             pos_b = bilayer_analyzer.com_frame.lipidcom[lindex_b].com
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
    #     tc = bilayer_analyzer.com_frame.time
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

        # required
        self.valid_args = ['interval', 'leaflet', 'resname', "wrapped"]
        self.return_length = 4
        self.analysis_key = 'disp_vec'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        self.leaflet = 'both'
        self.group = 'all'
        self.wrapped = False
        self.interval = 10
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)

        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            print ('arg_key: ', arg_key)
            if arg_key in self.valid_args:
                if arg_key == 'interval':
                    self.interval = int(arg_arg)
                elif arg_key == 'leaflet':
                    if arg_arg in ['upper', 'lower', 'both']:
                        self.leaflet = arg_arg
                elif arg_key == 'resname':
                    self.group = arg_arg
                elif arg_key == 'wrapped':
                    if arg_arg in ['True', 'true']:
                        self.wrapped = arg_arg
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("\'" + self.analysis_id + "\': displacement vectors for " + str(
            self.interval) + " frame intervals")
        return

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_analysis(self, bilayer_analyzer):

        if self.first_comp:
            self.last_com_frame = bilayer_analyzer.com_frame
            self.first_comp = False
            return
        current_frame = bilayer_analyzer.current_mda_frame.frame

        interval = (current_frame) - (self.last_frame)
        print (interval, " ", self.interval)
        if interval == self.interval:
            indices = []
            if self.leaflet == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.group)
            elif self.leaflet == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.group)

            elif self.leaflet == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.group)
            else:
                # unknown option--use default "both"
                raise RuntimeWarning(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.leaflet = 'both'
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.group)
            n_com = len(indices)

            # print "there are ",len(indices)," members"
            xi = bilayer_analyzer.lateral[0]
            yi = bilayer_analyzer.lateral[1]
            zi = bilayer_analyzer.norm

            vec_ends_out = []
            # get the current frame
            curr_frame = bilayer_analyzer.com_frame
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
                    bilayer_analyzer.lateral]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    bilayer_analyzer.lateral]
                com_j_w = prev_frame.lipidcom[i].com[bilayer_analyzer.lateral]
                if self.wrapped:
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
            self.last_com_frame = bilayer_analyzer.com_frame
            self.last_frame = current_frame
            #return vec_ends
        return

    def save_data(self):

        with open(self.save_file_name, 'wb') as outfile:
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

        # required
        self.valid_args = ['interval', 'leaflet', 'resname', "wrapped"]
        self.return_length = 4
        self.analysis_key = 'disp_vec'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        self.leaflet = 'both'
        self.group = 'all'
        self.wrapped = False
        self.interval = 10
        # parse input arguments if given
        if len(args) > 1:
            self._parse_args(args)

        # storage for output
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments
    def _parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            print ('arg_key: ', arg_key)
            if arg_key in self.valid_args:
                if arg_key == 'interval':
                    self.interval = int(arg_arg)
                elif arg_key == 'leaflet':
                    if arg_arg in ['upper', 'lower', 'both']:
                        self.leaflet = arg_arg
                elif arg_key == 'resname':
                    self.group = arg_arg
                elif arg_key == 'wrapped':
                    if arg_arg in ['True', 'true']:
                        self.wrapped = arg_arg
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("\'" + self.analysis_id + "\': displacement vectors for " + str(
            self.interval) + " frame intervals")
        return

    def reset(self):
        self.analysis_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_analysis(self, bilayer_analyzer):

        if self.first_comp:
            self.last_com_frame = bilayer_analyzer.com_frame
            self.first_comp = False
            return
        current_frame = bilayer_analyzer.current_mda_frame.frame

        interval = (current_frame) - (self.last_frame)
        print (interval, " ", self.interval)
        if interval == self.interval:
            indices = []
            if self.leaflet == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.group)
            elif self.leaflet == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.group)

            elif self.leaflet == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.group)
            else:
                # unknown option--use default "both"
                raise RuntimeWarning(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.group)
            n_com = len(indices)

            # print "there are ",len(indices)," members"
            xi = bilayer_analyzer.lateral[0]
            yi = bilayer_analyzer.lateral[1]
            zi = bilayer_analyzer.norm

            vec_ends_out = []
            # get the current frame
            curr_frame = bilayer_analyzer.com_frame
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
                    bilayer_analyzer.lateral]
                com_j = prev_frame.lipidcom[i].com_unwrap[
                    bilayer_analyzer.lateral]
                com_j_w = prev_frame.lipidcom[i].com[bilayer_analyzer.lateral]
                if self.wrapped:
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
            self.last_com_frame = bilayer_analyzer.com_frame
            self.last_frame = current_frame
            #return vec_ends
        return

    def save_data(self):

        with open(self.save_file_name, 'wb') as outfile:
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
        # required
        self.valid_args = None
        self.return_length = 4
        self.analysis_key = 'ndcorr'
        self.analysis_id = args[0]
        # default function settings
        self.save_file_name = self.analysis_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            warn = "ignoring extra arguments passed to " \
                   "ndcorr analysis '{}''".format(self.analysis_id)
            raise RuntimeWarning(warn)
        self.first_frame = True
        # storage for output
        self.analysis_output = dict()
        self.running_stats = dict()
        return

    # required- function to parse the input arguments
    def _parse_args():
        pass

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print (
            "\'" + self.analysis_id + "\': Norm dimension displacement lipid type cross correlation estimate.")
        return

    def run_analysis(self, bilayer_analyzer):
        #construct the grids
        grids = lgc.LipidGrids(bilayer_analyzer.com_frame,bilayer_analyzer.leaflets,bilayer_analyzer.lateral)

        if self.first_frame:
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
                for type in all_types:
                    self.analysis_output[leaf][type] = []
                    self.running_stats[leaf][type] = RunningStats()
            self.first_frame = False
        #analysis the correlations
        correlations = grids.norm_displacement_cross_correlation()
        time = bilayer_analyzer.current_mda_frame.time
        #extract the data
        leafs = correlations.keys()
        for leaf in leafs:
            for type in correlations[leaf].keys():
                corr = correlations[leaf][type]
                self.running_stats[leaf][type].push(corr)
                corr_run = self.running_stats[leaf][type].mean()
                corr_std = self.running_stats[leaf][type].deviation()
                self.analysis_output[leaf][type].append(np.array([time, corr, corr_run, corr_std]))

        return

    def reset(self):
        self.first_frame = True
        self.analysis_output = dict()
        self.running_stats = dict()

        return


command_protocols['ndcorr'] = NDCorrProtocol
