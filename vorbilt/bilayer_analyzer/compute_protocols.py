# imports
import scipy.constants as scicon

try:
    import cPickle as pickle
except:
    import pickle

# ORBILT imports
from vorbilt.common.running_stats import *
import vorbilt.mda_tools.mda_density_profile as mda_dp

command_protocols = {}
valid_computes = []
compute_obj_name_dict = {}
# obj_dict = {"com_frame":COMFrame}

use_objects = {"mda_frame": True, "com_frame": False, "lipid_grid": False}


# TO DO:
# have protocols parse dictionaries as well as strings, which is more Pythonic


# protocol for the computes to run during the frame loop
class ComputeProtocol:
    def __init__(self, compute_commands):
        self.use_objects = use_objects
        self.in_commands = compute_commands
        # check computes
        arguments = []
        compute_keys = []
        command_protocol = {}
        compute_ids = []
        for command in compute_commands:
            if len(command) < 2:
                raise RuntimeError(
                    "wrong number of arguments for compute " + str(command))
            comp_key = command[0]
            comp_id = command[1]
            comp_args = command[1:]
            if (comp_key in valid_computes):
                if (len(comp_args) >= 1):
                    if comp_id not in compute_ids:
                        comp_object = compute_obj_name_dict[comp_key]
                        self.use_objects[comp_object] = True
                        arguments.append(comp_args)
                        compute_ids.append(comp_args[0])
                        compute_keys.append(comp_key)
                        command_protocol[comp_args[0]] = command_protocols[
                            comp_key](comp_args)
                    else:
                        raise RuntimeError(
                            "compute id \'" + comp_id + "\' has already been used")
                else:
                    raise RuntimeError(
                        "wrong number of arguments for compute " + str(
                            command))
            else:
                raise RuntimeError(
                    "invalid compute \'" + comp_key + "\' :" + str(command))
        self.arguments = arguments
        self.compute_keys = compute_keys
        self.compute_ids = compute_ids
        self.command_protocol = command_protocol
        self.n_commands = len(command_protocol)
        # object dependencies
        if self.use_objects['lipid_grid']:
            self.use_objects['com_frame'] = True
        return

    def add_compute(self, compute_string):
        command = compute_string.split()
        comp_key = command[0]
        comp_id = command[1]
        comp_args = command[1:]
        if (comp_key in valid_computes):
            if (len(comp_args) >= 1):
                if comp_id not in self.compute_ids:
                    comp_object = compute_obj_name_dict[comp_key]
                    self.use_objects[comp_object] = True
                    self.arguments.append(comp_args)
                    self.compute_ids.append(comp_args[0])
                    self.compute_keys.append(comp_key)
                    self.command_protocol[comp_args[0]] = command_protocols[
                        comp_key](comp_args)
                else:
                    raise RuntimeError("compute id '{}' "
                                       "has already been used!".format(comp_id))
            else:
                raise RuntimeError("wrong number of arguments "
                                   "for compute {}".format(command))
        else:
            raise RuntimeError("invalid compute id"
                               " '{}' : {}".format(comp_key, command))
        self.n_commands += 1
        return

    def remove_compute(self, compute_id):
        if compute_id in self.compute_ids:
            del self.command_protocol[compute_id]
            index = self.compute_ids.index(compute_id)
            del self.arguments[index]
            del self.compute_keys[index]
            del self.compute_ids[index]
            self.n_commands -= 1
        else:
            raise RuntimeWarning("no compute with id '{}'".format(compute_id))
        return

    def print_protocol(self):
        print ('build objects:')
        for key in self.use_objects.keys():
            if self.use_objects[key]:
                print (key)
        print ("with computes:")
        for compute_id in self.compute_ids:
            self.command_protocol[compute_id].print_protocol()
        return

    def dump_data(self):
        print ('dumping compute data to pickle files...')
        for compute_id in self.compute_ids:
            print ("compute id: {} ---> {} ".format(
                compute_id,
                self.command_protocol[compute_id].save_file_name))
            self.command_protocol[compute_id].save_data()

    def reset(self):
        for compute_id in self.compute_ids:
            self.command_protocol[compute_id].reset()
        return


class ComputeFunctionProtocol:
    def __init__(self, args):

        # required
        self.valid_args = ['none']
        self.return_length = 1
        self.compute_key = 'none'
        self.compute_id = args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)
        # storage for output
        self.compute_output = []
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        # print args
        nargs = len(args)
        for i in range(1, nargs, 2):
            arg_key = args[i]
            arg_arg = args[i + 1]
            if arg_key in self.valid_args:
                pass
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for compute " + self.compute_key)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("Parent protocol class for computes.")
        return

    def run_compute(self, bilayer_analyzer):
        # do some stuff
        # get an output
        output = np.zeros(self.return_length)
        # save the output
        self.compute_output.append(output)
        return

    def save_data(self):
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(np.array(self.compute_output), outfile)

        return

    def get_data(self):
        return np.array(self.compute_output)

    def reset(self):
        self.compute_output = []
        return


# define a new compute 'msd'
valid_computes.append('msd')
compute_obj_name_dict['msd'] = 'com_frame'


class MSDProtocol(ComputeFunctionProtocol):
    def __init__(self, msd_args):

        # required
        self.valid_args = ['leaflet', 'resname']
        self.return_length = 2
        self.compute_key = 'msd'
        self.compute_id = msd_args[0]
        # default function settings
        self.leaflet = 'both'
        self.resname = 'all'
        self.save_file_name = self.compute_id + ".pickle"
        self.first_frame = True
        self.ref_coords = None
        self.indices = []
        # parse input arguments if given
        if len(msd_args) > 1:
            self.parse_args(msd_args)
        # storage for output
        self.compute_output = []
        return

    # required- function to parse the input arguments
    def parse_args(self, msd_args):
        # print msd_args
        nargs = len(msd_args)
        for i in range(1, nargs, 2):
            arg_key = msd_args[i]
            arg_arg = msd_args[i + 1]
            if arg_key in self.valid_args:
                if arg_key == 'leaflet':
                    self.leaflet = arg_arg
                elif arg_key == 'resname':
                    self.resname = arg_arg
            else:
                raise RuntimeError(
                    "ignoring invalid argument key " + arg_key + " for compute \'msd\'")
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
            "\'" + self.compute_id + "\': msd analysis of " + self.resname + " lipids in " + self.leaflet + " leaflet(s).")
        return

    def run_compute(self, bilayer_analyzer):
        leaflet = self.leaflet
        group = self.resname
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
        self.compute_output.append(msd)
        return


# update the command_protocols dictionary
command_protocols['msd'] = MSDProtocol

# define a new compute 'apl_box'
valid_computes.append('apl_box')
compute_obj_name_dict['apl_box'] = 'mda_frame'


class APLBoxProtocol(ComputeFunctionProtocol):
    def __init__(self, apl_args):
        # required
        self.valid_args = None
        self.return_length = 4
        self.compute_key = 'apl_box'
        self.compute_id = apl_args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        # parse input arguments if given
        if len(apl_args) > 1:
            raise RuntimeWarning(
                "ignoring extra arguments passed to area per lipid compute \'" + self.compute_id + "\'")

        # storage for output
        self.running = RunningStats()
        self.compute_output = []
        return

    # required- function to parse the input arguments
    @staticmethod
    def parse_args():
        pass

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print ("'{}'': area per lipid using box "
               "dimensions".format(self.compute_id))
        return

    def run_compute(self, bilayer_analyzer):
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
        self.compute_output.append(apl_t)

        return

    def reset(self):
        self.running.reset()
        self.compute_output = []
        return


command_protocols['apl_box'] = APLBoxProtocol

# define a new compute 'apl_box'
valid_computes.append('bilayer_thickness')
compute_obj_name_dict['bilayer_thickness'] = 'lipid_grid'


class BTGridProtocol(ComputeFunctionProtocol):
    def __init__(self, args):
        # required
        self.valid_args = None
        self.return_length = 4
        self.compute_key = 'bilayer_thickness'
        self.compute_id = args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            warn = "ignoring extra arguments passed to " \
                   "bilayer thickness compute '{}''".format(self.compute_id)
            raise RuntimeWarning(warn)

        # storage for output
        self.running = RunningStats()
        self.compute_output = []
        return

    # required- function to parse the input arguments
    def parse_args():
        pass

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print (
            "\'" + self.compute_id + "\': bilayer thickness using lipid_grid")
        return

    def run_compute(self, bilayer_analyzer):
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
        self.compute_output.append(bt_t)

        return

    def reset(self):
        self.running.reset()
        self.compute_output = []
        return


command_protocols['bilayer_thickness'] = BTGridProtocol

# define a new compute 'apl_grid'
valid_computes.append('apl_grid')
compute_obj_name_dict['apl_grid'] = 'lipid_grid'


class APLGridProtocol(ComputeFunctionProtocol):
    def __init__(self, args):

        # required
        self.valid_args = None
        self.return_length = 4
        self.compute_key = 'apl_grid'
        self.compute_id = args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)

        # storage for output
        self.running = RunningStats()
        self.compute_output = {}
        self.first_comp = True
        self.running_res = {}
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
        raise RuntimeWarning(
            "ignoring extra arguments passed to area per lipid (by grid) compute \'" + self.compute_id + "\'")
        return

    # required - a check protocol function which reports relevant settings
    def print_protocol(self):
        print ("\'" + self.compute_id + "\': area per lipid using lipid_grid")
        return

    def run_compute(self, bilayer_analyzer):
        apl_grid_out = bilayer_analyzer.lipid_grid.area_per_lipid()
        # print (current_thickness)
        time = bilayer_analyzer.current_mda_frame.time
        # print (time)
        self.running.push(apl_grid_out[0])
        if self.first_comp:
            self.compute_output['composite'] = []
            for key in apl_grid_out[1].keys():
                self.running_res[key] = RunningStats()
                self.first_comp = False
                self.compute_output[key] = []

        apl_t = np.zeros(4)
        apl_t[0] = time
        apl_t[1] = apl_grid_out[0]
        apl_t[2] = self.running.mean()
        apl_t[3] = self.running.deviation()

        self.compute_output['composite'].append(apl_t)

        for key in apl_grid_out[1].keys():
            self.running_res[key].push(apl_grid_out[1][key][0])
            apl_res = np.zeros(4)
            apl_res[0] = time
            apl_res[1] = apl_grid_out[1][key][0]
            apl_res[2] = self.running_res[key].mean()
            apl_res[3] = self.running_res[key].deviation()
            self.compute_output[key].append(apl_res)

        return

    def save_data(self):
        for key in self.compute_output.keys():
            self.compute_output[key] = np.array(self.compute_output[key])
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(self.compute_output, outfile)

        return

    def get_data(self):
        for key in self.compute_output.keys():
            self.compute_output[key] = np.array(self.compute_output[key])
        return self.compute_output

    def reset(self):
        self.running.reset()
        self.compute_output = []
        self.first_comp = True
        self.running_res = {}
        return


command_protocols['apl_grid'] = APLGridProtocol

# define a new compute 'disp_vec'
valid_computes.append('disp_vec')
compute_obj_name_dict['disp_vec'] = 'com_frame'


class DispVecProtocol(ComputeFunctionProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['interval', 'leaflet', 'resname', "wrapped"]
        self.return_length = 4
        self.compute_key = 'disp_vec'
        self.compute_id = args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        self.leaflet = 'both'
        self.group = 'all'
        self.wrapped = False
        self.interval = 10
        # parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)

        # storage for output
        self.compute_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    # required- function to parse the input arguments
    def parse_args(self, args):
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
                    "ignoring invalid argument key " + arg_key + " for compute" + self.compute_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print ("\'" + self.compute_id + "\': displacement vectors for " + str(
            self.interval) + " frame intervals")
        return

    def reset(self):
        self.compute_output = []
        self.first_comp = True
        self.last_com_frame = None
        self.last_frame = 0

        return

    def run_compute(self, bilayer_analyzer):

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
                    "bad setting for \'leaflet\' in " + self.compute_id + ". Using default \'both\'")
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
            self.compute_output.append([vec_ends, resnames])
            self.last_com_frame = bilayer_analyzer.com_frame
            self.last_frame = current_frame
            return vec_ends
        return

    def save_data(self):

        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(self.compute_output, outfile)

        return

    def get_data(self):
        return self.compute_output


command_protocols['disp_vec'] = DispVecProtocol

# define a new compute
valid_computes.append('mass_dens')
compute_obj_name_dict['mass_dens'] = 'mda_frame'


class MassDensProtocol(ComputeFunctionProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['selection', 'n_bins']
        self.return_length = None
        self.compute_key = 'mass_dens'
        self.compute_id = args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        self.selection_string = 'all'
        self.n_bins = 25
        # parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)

        # storage for output
        self.centers = None
        self.n_frames = 0
        self.compute_output = []
        self.first_comp = True
        self.selection = None

        return

    # required- function to parse the input arguments
    def parse_args(self, args):
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
                    "ignoring invalid argument key " + arg_key + " for compute " + self.compute_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print (
            "\'" + self.compute_id + "\': mass density compute for selection \'" + self.selection_string + "\' ")
        return

    def reset(self):
        self.centers = None
        self.n_frames = 0
        self.compute_output = []
        self.first_comp = True
        self.selection = None

        return

    def run_compute(self, bilayer_analyzer):
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
            self.compute_output = np.zeros(len(self.centers))
        self.compute_output += centers_density[1]
        self.n_frames += 1
        return

    def save_data(self):
        centers_density = (self.centers, self.compute_output / self.n_frames)
        with open(self.save_file_name, 'wb') as outfile:
            pickle.dump(centers_density, outfile)

        return

    def get_data(self):
        return (self.centers, self.compute_output / self.n_frames)


command_protocols['mass_dens'] = MassDensProtocol


# define a new compute
valid_computes.append('acm')
compute_obj_name_dict['acm'] = 'mda_frame'


class AreaCompressibilityModulusProtocol(ComputeFunctionProtocol):
    def __init__(self, args):

        # required
        self.valid_args = ['temperature']
        self.return_length = 4
        self.compute_key = 'acm'
        self.compute_id = args[0]
        # default function settings
        self.save_file_name = self.compute_id + ".pickle"
        # parse input arguments if given
        if len(args) > 1:
            self.parse_args(args)

        # storage for output
        self.temperature = 298.15
        self.n_frames = 0
        self.area = []
        self.apl = []
        self.times = []
        self.compute_output = []
        self.first_comp = True
        self.per_leaflet = 1
        return

    # required- function to parse the input arguments
    def parse_args(self, args):
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
                    "ignoring invalid argument key " + arg_key + " for compute " + self.compute_id)
        return
        # required - a check protocol function which reports relevant settings

    def print_protocol(self):
        print("\'" + self.compute_id + "\': area compressibility modulus compute for selection \'")
        return

    def reset(self):
        self.area_run.reset()
        self.area = []
        self.apl = []
        self.n_frames = 0
        self.compute_output = []
        self.first_comp = True
        return

    def run_compute(self, bilayer_analyzer):

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
