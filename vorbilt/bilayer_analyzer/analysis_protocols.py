"""Analysis Protocols

This is a support module that defines a set of classes used to contruct the 'analysis protocol' and the 'analysis' used
in the anlysis implemented by the vorbilt.bilayer_analyzer.bilayer_anlayzer.BilayerAnalyzer class.
The AnalysisProtocol is the class used to organize and initial the the individual analysis functions/protocols. The
individual analysis functions/protocols are derived classes of AnalysisProtocol, and the availabel analysis can be
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
from vorbilt.common import distance_cutoff_clustering as dc_cluster
#need some containers for bookkeeping
command_protocols = {}
valid_analysis= []
analysis_obj_name_dict = {}
# obj_dict = {"com_frame":COMFrame}
#define the buildable objects that are used by the analysisfunctions
use_objects = {"mda_frame": True, "com_frame": True, "lipid_grid": False}


# TO DO:
# have protocols parse dictionaries as well as strings, which is more Pythonic


def word_list_to_string(word_list, delimeter=" "):
    string = ""
    for word in word_list:
        string+=word+delimeter
    nchar = len(string)
    return str(string[0:nchar-1])

# protocol for the analysis to run during the frame loop
class Analyses:
    """A class to facilitate analysis of the bilayers
    This object stores all the analyses that being performed and provides fucntionality to add and remove
    analyses.

    Attributes:
        use_objects (list of str): A list of buildable objects that need to contructed in the BilayerAnalyzer for the
            analysis defined in this protocol.
        in_commands (list): A list of the input strings for the analysis to be used.
        arguments (list): A list of the arguments for analysis.
        analysis_keys (list): A list of the keys assigned to analysis.
        command_protocol (dict): A dictionary of the analysis objects.
        analysis_ids (list): A list of the ids assigned to analysis.
        n_commands (int): The number of initialized analysis.

    """
    def __init__(self, analysis_commands):
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
        return self.command_protocol[item]

    def __len__(self):
        return len(self.analysis_ids)

    def add_analysis(self, inputs):
        if isinstance(inputs, (str, basestring)):
           self._add_analysis_from_string(inputs)
        elif isinstance(inputs, (list, tuple)):
            self._add_analysis_from_list(inputs)
        elif isinstance(inputs, dict):
           self._add_analysis_from_dict(inputs)
        return

    def _add_analysis_from_string(self, analysis_string):
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

    def dump_data(self, path=None):
        print ('dumping analysis data to pickle files...')
        for analysis_id in self.analysis_ids:
            print ("analysis id: {} ---> {} ".format(
                analysis_id,
                self.command_protocol[analysis_id].save_file_name))
            self.command_protocol[analysis_id].save_data(path=path)

    def reset(self):
        for analysis_id in self.analysis_ids:
            self.command_protocol[analysis_id].reset()
        return

#base class for analysis protocols
class AnalysisProtocol:
    """Base class for analysis protocols.

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
    """
    def __init__(self, args):

        # required
        self._short_description = "parent analysis protocol"
        self.return_length = 1
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
        """Parses the setup arguments for this analysis
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
        arg_dict = self._parse_str_to_dict(args)
        #type cast setttings if needed
        arg_dict = self._cast_settings(arg_dict)
        self._parse_dict(arg_dict)
        return

    def _parse_str_to_dict(self, args):
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
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_key)
        return arg_dict

    # cast the input string values of settings to appropriate types -- should be overwritten in derived classes to
    # properly type cast their own settings.
    def _cast_settings(self, args_dict):
        for setting_key in args_dict:
            pass
        return args_dict

    def _parse_dict(self, args):
        if 'analysis_id' not in args.keys():
            raise RuntimeError("required key \'anlaysis_id\' not assigned in input dict for analysis type: \'"+self.analysis_key+"\'")
        for arg_key in args.keys():
            arg_arg = args[arg_key]
            if arg_key in self._valid_settings:
                self.settings[arg_key] =  arg_arg
            elif arg_key == 'analysis_id':
                self.analysis_id = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_key)
        return


    def print_protocol(self):
        print ("Analysis: "+self._short_description)
        print ("  with analysis_id: {} ".format(self.analysis_id))
        print ("   and settings: ")
        for key in self.settings.keys():
            print ("    {}: {} ".format(key, self.settings[key]))
        return

    def run_analysis(self, bilayer_analyzer):
        # do some stuff
        # get an output
        output = np.zeros(self.return_length)
        # save the output
        self.analysis_output.append(output)
        return

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        with open(save_file, 'wb') as outfile:
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
    """Mean squared displacement
    Args:
        args (list): list of string keys and arguments
    """
    def __init__(self, args):

        # required
        self._short_description = "Mean squared displacement."
        self.return_length = 2
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

        self.save_file_name = self.analysis_id + ".pickle"
        self.first_frame = True
        self.ref_coords = None
        self.indices = []
        # parse input arguments if given
        self._parse_args(args)
        # storage for output
        self.analysis_output = []
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
        self._short_description = "Area per lipid using box dimensions."
        self.return_length = 4
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
        self._short_description = "Bilayer thickness using lipid_grid."
        self.return_length = 4
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
        self._short_description = "Area per lipid using lipid_grid"
        self.return_length = 4
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
        self._short_description = "Displacement vectors."
        self.return_length = 4
        self.analysis_key = 'disp_vec'
        self.analysis_id = 'none'

        #default settings
        self.settings = dict()
        self.settings['leaflet'] = 'both'
        self.settings['resname'] = 'all'
        self.settings['wrapped'] = False
        self.settings['interval'] = 5
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
    def _cast_settings(self, args_dict):

        for arg_key in args_dict.keys():
            arg_arg = args_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    args_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return args_dict


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
        print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                raise RuntimeWarning(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
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
            self.last_com_frame = bilayer_analyzer.com_frame
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
    def __init__(self, args):

        # required
        self._short_description = "Mass density."
        self.return_length = None
        self.analysis_key = 'mass_dens'
        self.analysis_id = 'none'

        # default settings
        self.settings = dict()
        self.settings['selection_string'] = 'all'
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
    def _cast_settings(self, args_dict):
        #in this case the casting is taken care of in _parse_str_to_dict
        return args_dict

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
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis " + self.analysis_id)
        return arg_dict


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
            if self.settings['selection_string'] == ' BILAYER':
                self.selection = bilayer_analyzer.mda_data.bilayer_sel
            else:
                self.selection = bilayer_analyzer.mda_data.mda_universe.select_atoms(
                    self.settings['selection_string'])
            self.first_comp = False

        # print "there are ",len(indices)," members"

        norm_axis = bilayer_analyzer.normal_dimension
        ref_sel = bilayer_analyzer.mda_data.bilayer_sel
        # mda_density_profile
        centers_density = mda_dp.mass_density_profile(
            bilayer_analyzer.mda_data.mda_trajectory, self.selection,
            fstart=bilayer_analyzer.frame_index,
            fend=bilayer_analyzer.frame_index + 1, axis=norm_axis,
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


command_protocols['mass_dens'] = MassDensProtocol


# define a new analysis
valid_analysis.append('acm')
analysis_obj_name_dict['acm'] = 'mda_frame'

#L. Janosi and A. A. Gorfe, J. Chem. Theory Comput. 2010, 6, 3267-3273
#D. Aguayo, F. D. Gonzalez-Nilo, and C. Chipot, J. Chem. Theory Comput. 2012, 8, 1765-1773
class AreaCompressibilityModulusProtocol(AnalysisProtocol):
    def __init__(self, args):

        # required
        self._short_description = "Area compressibility modulus."
        self.return_length = 4
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
        self.area = []
        self.apl = []
        self.times = []
        self.analysis_output = []
        self.first_comp = True
        self.per_leaflet = 1
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, args_dict):

        for arg_key in args_dict.keys():
            arg_arg = args_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'temperature':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return args_dict


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

    def save_data(self, path=None):
        save_file = self.save_file_name
        if path is not None:
            save_file = path+self.save_file_name
        data = self.get_data()
        with open(save_file, 'wb') as outfile:
            pickle.dump(data, outfile)
        return

    def get_data(self):
        area_eq = np.array(self.area).mean()
        apl = np.array(self.apl)
        apl_run = gen_running_average(apl)
        apl_minus_area = (apl - area_eq)**2
        apl_minus_area_run = gen_running_average(apl_minus_area)
        acm = scicon.k * self.settings['temperature'] * apl_run[:,0] / ( self.per_leaflet * apl_minus_area_run[:,0])
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
        self._short_description = "Lateral order nearest neighbor fraction."

        self.return_length = 2
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
        self.first_frame = True
        self.ref_coords = None
        self.indices = []
        self.running = RunningStats()
        # storage for output
        self.analysis_output = []
        return

    # required- function to parse the input arguments from string
    def _cast_settings(self, args_dict):

        for arg_key in args_dict.keys():
            arg_arg = args_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'n_neighbors':
                    args_dict[arg_key] = int(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return args_dict

    def reset(self):
        self.running.reset()
        self.analysis_output = []
        self.first_frame = True
        return





    def run_analysis(self, bilayer_analyzer):
        do_leaflet = []
        if self.settings['leaflet'] == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.settings['leaflet']]

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
        x_index = bilayer_analyzer.lateral[0]
        y_index = bilayer_analyzer.lateral[1]
        box = bilayer_analyzer.current_mda_frame.dimensions[0:3]
        box_x = box[x_index]
        box_y = box[y_index]
        box_x_h = box_x / 2.0
        box_y_h = box_y / 2.0
        ltype_a = self.settings['resname_1']
        ltype_b = self.settings['resname_2']
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
        self._short_description = "Displacement vector correlation matrix."

        self.return_length = 4
        self.analysis_key = 'disp_vec'
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
    def _cast_settings(self, args_dict):

        for arg_key in args_dict.keys():
            arg_arg = args_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    args_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return args_dict

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
        print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                raise RuntimeWarning(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
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

        # required
        self._short_description = "Displacement vector nearest neigbor correlations."
        self.return_length = 4
        self.analysis_key = 'disp_vec'
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
    def _cast_settings(self, args_dict):

        for arg_key in args_dict.keys():
            arg_arg = args_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'interval':
                    arg_dict[arg_key] = int(arg_arg)
                elif arg_key == 'wrapped':
                    arg_arg = arg_arg in ['True', 'true']
                    args_dict[arg_key] = arg_arg
            elif arg_key == 'analysis_id':
                pass
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return args_dict

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
        print (interval, " ", self.settings['interval'])
        if interval == self.settings['interval']:
            indices = []
            if self.settings['leaflet'] == "both":
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
            elif self.settings['leaflet'] == "upper":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.settings['resname'])

            elif self.settings['leaflet'] == "lower":
                curr_leaf = bilayer_analyzer.leaflets[leaflet]
                indices = curr_leaf.get_group_indices(self.settings['resname'])
            else:
                # unknown option--use default "both"
                raise RuntimeWarning(
                    "bad setting for \'leaflet\' in " + self.analysis_id + ". Using default \'both\'")
                self.settings['leaflet'] = 'both'
                for leaflets in bilayer_analyzer.leaflets:
                    curr_leaf = bilayer_analyzer.leaflets[leaflets]
                    indices += curr_leaf.get_group_indices(self.settings['resname'])
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
            self.last_com_frame = bilayer_analyzer.com_frame
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

#
class NDCorrProtocol(AnalysisProtocol):
    def __init__(self, args):
        # required
        self._short_description = "Normal dimension displacement-lipid type cross correlation."

        self.return_length = 4
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
        self.first_frame = True
        # storage for output
        self.analysis_output = dict()
        self.running_stats = dict()
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
        # required
        self._short_description = "Distance cutoff clustering."

        self.return_length = 4
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
        self.first_frame = True
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
    def _cast_settings(self, args_dict):

        for arg_key in args_dict.keys():
            arg_arg = args_dict[arg_key]
            if arg_key in self._valid_settings:
                if arg_key == 'cutoff':
                    arg_dict[arg_key] = float(arg_arg)
            elif arg_key == 'analysis_id':
                pass
            else:
                raise RuntimeWarning(
                    "ignoring invalid argument key " + arg_key + " for analysis" + self.analysis_id)
        return args_dict

    def run_analysis(self, bilayer_analyzer):
        do_leaflet = []
        if self.settings['leaflet'] == 'both':
            do_leaflet = ['upper', 'lower']
        else:
            do_leaflet = [self.settings['leaflet']]

        if self.first_frame:
            #pass
            # build group/resname/lipid type list
            lipid_types = []
            nlipids = 0
            for leaflet_name in do_leaflet:
                leaflet = bilayer_analyzer.leaflets[leaflet_name]
                groups = leaflet.get_group_names()
        #        nlipids += len(leaflet.get_member_indices())
                for group in groups:
                    if group not in lipid_types:
                        lipid_types.append(group)
            if self.settings['resname'] == 'first':
                self.settings['resname'] =  lipid_types[0]
            self.first_frame = False
        indices = []
        for leaflets in do_leaflet:
            curr_leaf = bilayer_analyzer.leaflets[leaflets]
            indices += curr_leaf.get_group_indices(self.settings['resname'])
        pos = []
        for index in indices:
            pos.append(bilayer_analyzer.com_frame.lipidcom[index].com[bilayer_analyzer.lateral])
        pos = np.array(pos)
        box = bilayer_analyzer.com_frame.box[bilayer_analyzer.lateral]
        cutoff = self.settings['cutoff']
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
        avg_size /= nclusters
        self.running_stats['nclusters'].push(nclusters)
        self.running_stats['max_size'].push(max_size)
        self.running_stats['min_size'].push(min_size)
        self.running_stats['avg_size'].push(avg_size)
        time = bilayer_analyzer.com_frame.time
        self.analysis_output['clusters'].append(clusters)
        self.analysis_output['nclusters'].append([time,nclusters,self.running_stats['nclusters'].mean(), self.running_stats['nclusters'].deviation()])
        self.analysis_output['max_size'].append([time,max_size,self.running_stats['max_size'].mean(), self.running_stats['max_size'].deviation()])
        self.analysis_output['min_size'].append([time,min_size,self.running_stats['min_size'].mean(), self.running_stats['min_size'].deviation()])
        self.analysis_output['avg_size'].append([time,avg_size,self.running_stats['avg_size'].mean(), self.running_stats['avg_size'].deviation()])

        return

    def reset(self):
        self.first_frame = True
        self.analysis_output = dict()
        self.running_stats = dict()
        self.analysis_output['clusters'] = []
        self.analysis_output['nclusters'] = []
        self.analysis_output['max_size'] = []
        self.analysis_output['min_size'] = []
        self.analysis_output['avg_size'] = []
        self.running_stats['nclusters'].reset()
        self.running_stats['max_size'].reset()
        self.running_stats['min_size'].reset()
        self.running_stats['avg_size'].reset()
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
