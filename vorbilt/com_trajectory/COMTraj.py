"""Center of Mass based bilayer analysis tools

This module defines various classes and functions used to process and analyze a lipid 
bilayer trajectory. This module assumes the structure and trajectory are initiallaly stored in MDAnalysis
objects and therefore processes MDAnalysis objects. The lipids constituting the bilayer are read in from
the MDAnalysis objects and are converted to center of mass (COM) representations. Lipids are partitioned
into an 'upper' and a 'lower' leaflet based on the z-position of the COM. The built-in analysis functions
then operate on the COM representations to compute quantities such as the lateral mean squared displacement.  
Many analysis functions allow specification of the leaflet and type of lipid to perform the the analysis on.
The primary (parent class) is the MemSys class. The analysis functions are members of the MemSys class.

Example:
    >> import MemSys as ms
    >> mem_sys = ms.MemSys(mda_universe.trajectory,mda_selection_of_bilayer_lipids)
"""

#imports
import numpy as np
import matplotlib.cm as cm
import os
import sys
import shutil
import shelve
import multiprocessing as mp
from scipy.spatial import Voronoi
from scipy.spatial import Delaunay
#import copy

#import running stats class
from orbilt.common.running_stats import *
# import the coordinate wrapping function--for unwrapping
from orbilt.mda_tools.mda_unwrap import mda_wrap_coordinates,mda_wrap_coordinates_parallel

# assumes that a 1d numpy array of floats is pass as input, but 
# does not check this
def gen_running_average(onednparray):
    """ Generates a running average 
                 
    Args:
    onednparray (numpy.array): A 1d numpy array of measurements (e.g. over time)     
    
    Returns:
    numpy.array: 2d array of dim len(onednparray)x2
        2dnparray[i][0] = running average at i
        2dnparray[i][1] = running standard deviation at i    
        for i in range(0,len(onednparray))
    """
    averager = RunningStats()
    nele = len(onednparray)
    output = np.zeros((nele,2))
    for i in xrange(nele):
        averager.push(onednparray[i])
        run_avg = averager.mean()
        run_dev = averager.deviation()
       # print run_avg, run_dev, averager.mean(), onednparray[i]
        output[i,0] = run_avg
        output[i,1] = run_dev
    return output

# This function is incomplete!
def colorize_step_vector_clusters(vectors):
    nvecs = len(vectors)
    np.zeros(nvecs,dtype=np.int)
    colors_out = np.zeros(nvecs)
    return "nothing yet!"
    
    
class LipidCOM:
    """ A lipid center of mass (COM) object.  

    This object stores the COM coordinates of a lipid (or other molecule or group
    of atoms) computed from both the wrapped and unwrapped atomic coordinates. This
    object also stores information about the type of lipid as well as the total mass
    of the lipid.            
    """
    def __init__(self):
        """LipidCOM initialization 
 
        Attributes:
            type (str): The lipid type (e.g. the lipid could be typed b resname).
            com (np.array): The length three vector holding the wrapped xyz coordinates.
            com_unwrap (np.array): The length three vector holding the unwrapped xyz coordinates.
            mass (float): The total mass of the atoms used to define LipidCOM. 
        """
        # lipid type/resname or other name    
        self.type="UNK"
        # wrapped coordinates
        self.com=np.zeros(3)
        # unwrapped coordinates
        self.com_unwrap=np.zeros(3)
        # total mass
        self.mass=1.0
        return
    # The name of this function could be changed to be more desriptive, e.g. 
    # extract_com_mda_residue  
    def extract(self, mda_residue, unwrap=False, box=None):
        """ Get the center of mass coordinates from an MDAnalysis residue  

        This function calls the MDAnalysis member function center_of_mass() of the residue
        to compute the center of mass of the atoms constituting the residue.

        Args:
            mda_residue (MDAnalysis.residue): An MDAnalysis residue object from which to extract
                a center of masss coordinates.
            unwrap (bool, optional): Define which com container to store coordiates in.
                False (default) - The COM coordinates are stored in the 
                container designated for the wrapped coordinate representation.
                True - The COM coordinates are stored in the container designated
                for the unwrapped coordinate representation.
        """
        if unwrap:
            self.com_unwrap = mda_residue.center_of_mass()

        else:
            if box is not None:
                self.com = mda_residue.center_of_mass()
                self.com_unwrap = self.com[:]
            else:
                self.com = mda_residue.center_of_mass()
                self.com_unwrap = self.com[:]
        
        self.type=mda_residue.resname
        return

# a frame object 
class Frame:
    """ A molecular dynamics style Frame object. """    
                
    # does not check that nlipids is an int
    def __init__(self, nlipids):
        """ Frame initialization.    

        Args:
            nlipids (int): The number of lipids (LipidCOM objects) that this frame will hold

        Atrributes:
            lipidcom (list of obj:LipidCOM): A list of the LipidCOM objects assigned to the Frame.
            box (np.array): A 3 element vector containing the (rectangular) xyz box edge lengths. 
            time (float): The simulation time that this Frame represents.
            number (int): The frame number of this Frame.
            mdnumber (int): The corresponding frame number in the original MD trajectory
            
        """
        # list to store the nlipids LipidCOM objects
        self.lipidcom = []
        # box dimensions -- assumes the box originates a 0,0,0
        # It might be worth adding functionality to specifiy the box origin (or center)
        # This also assumes a rectangular box
        self.box = np.zeros(3)
        # simulation time
        self.time = np.zeros(1)
        # frame number
        self.number = np.zeros(1,dtype=np.int)
        # frame number in the MD trajectory
        self.mdnumber = np.zeros(1,dtype=np.int)
        # initialize all the LipidCOM objects
        for i in xrange(nlipids):
            self.lipidcom.append(LipidCOM())

        return
    
    def set_box(self, box_lengths):
        """ Set the rectangular xyz box edge lengths.    

        Args:
            box_lengths (numpy.array): A 1d, 3 element numpy.array containing the x,y,z box sizes (or edge lengths) 

        """
        self.box = box_lengths
        return

    def set_time(self, time):
        """ Set the simulation time.    

        Args: 
            time (float): The simulation time to assign to this Frame.

        """
        self.time = time
        return

    def __len__(self):
        """ Returns the number of LipidCOM objects assigned to this Frame
        
        Returns:
            int: Number of LipidCOM objects currently assigned to this Frame
        """
        return len(self.lipidcom)

#    def COG(self,unwrapped=False):
#        cog_out = np.zeros(3)
#        for lipid in self.lipidcom:    
#            if not unwrapped:
#                cog_out+=lipid.com    
#            else:
#                cog_out+=lipid.com_unwrap
#        cog_out/=len(self)
#        return com_out
    
    def com(self, wrapped=True):
        """ Computes the center of mass (COM) for the Frame    
        
        This member function is used to compute the overall center of mass (COM) of the 
        Frame using the LipidCOM object coordinates and masses. 

        Args: 
            wrapped (bool, optional): Define which set of coordinates to use in the computation.
                True (default) - The wrapped LipidCOM coordinates are used to compute 
                the COM of the frame.
                False - The unwrapped LipidCOM coordinates are used to compute 
                the COM of the frame.
                       
        Returns:  
            np.array: A 3 element vector containing the xyz coordinates of the Frame's COM
        """
        com_out = np.zeros(3)
        total_mass = 0.0
        for lipid in self.lipidcom:    
            if wrapped:
                com_out+=lipid.com*lipid.mass
                total_mass+=lipid.mass    
            else:
                com_out+=lipid.com_unwrap*lipid.mass
                total_mass+=lipid.mass    
        com_out/=total_mass
        return com_out
    
    def write_xyz(self, xyz_name, wrapped=True):
        # Open up the file to write to
        xyz_out = open(xyz_name, "w")
        
        comment = "Memsys Frame "+str(self.number)+" MD Frame "+str(self.mdnumber)
        xyz_out.write(str(len(self.lipidcom)))
        xyz_out.write("\n")
        xyz_out.write(comment)
        xyz_out.write("\n")
        
        
        
        i=0
        for lip in self.lipidcom:
            #get the coordinates
            x = self.lipidcom[i].com[0]
            y = self.lipidcom[i].com[1]
            z = self.lipidcom[i].com_unwrap[2]
            if not wrapped:
               x = self.lipidcom[i].com_unwrap[0]
               y = self.lipidcom[i].com_unwrap[1]
               #z = self.lipidcom[i].com_unwrap[2] 
                
            #get the lipid resname
            oname = self.lipidcom[i].type
            
            #write to file
            line = str(oname)+" "+str(x)+" "+str(y)+" "+str(z)
            xyz_out.write(line)
            xyz_out.write("\n")
            i+=1
                 
        xyz_out.close()
        return

        

#frame wrapper - the name of this class may be changed. e.g. FrameShelve
class FrameShelve:
    """ Container for Frame objects    
    This class object serves as a container to store a set of Frame objects
    corresponding to a molecular dynamics trajectory. This class saves the Frame objects
    on disk using the shelve module and provides an interface to access instances of
    those saved Frames. The Frames are saved in the shelve database with integer index keys
    0 -> (nframes-1) and can accessed by indexing the instance of the frames object. 
                    
    """
    #define a non-instance type error message
    _type_error ="Instance of object MemSys.FrameShelve only excepts instances of MemSys.Frame."
    
    def __init__(self,prefix='/tmp/',save=False):
        r""" Initialization of the frames object.  

        Args:
            prefix (str, optional): The location/path to store the "shelve"d Frame data.   
                '/tmp/' (default) - The data is stored in the unix/linux tmp directory.
            save (bool, optional): Set whether to save the shelved Frame data after frames object deletion. 
                False (default) - the shelved Frame data is deleted upon calling __del__ 
                True  - the shelved Frame data is not deleted when __del__ is called 
        
        Attributes:
            nframes (int): The number Frame objects being stored.
            pid (int): The process id that the frames object is created in.
            path (str): The path to where the shelved Frame data is stored.
            save (bool): Set whether to save the shelved Frame data after frames object deletion. 
            fs_name (str): The base name of the shelve database file used to store Frame objects.
            frame_shelf (shelve.Shelf): The Shelf that will hold the Frame objects.
        """
        self.nframes = 0
        self.pid = os.getpid()
        
        if    prefix[-1] != '/':
            prefix = prefix +'/'
        path = prefix
        if    save:
            path = path+'mem_sys_frames'
        else:
            path = path+'.mem_sys_frames_'+str(self.pid)            
        self.path = path
        self.save = save
        if os.path.isdir(self.path):
            shutil.rmtree(self.path)
        os.mkdir(self.path, 0755)
        self.fs_name = self.path +'/shelf_frames.db' 
        self.frame_shelf = shelve.open(self.fs_name,flag="c", protocol=2)
        return

    def __del__(self):
        """ Non-standard implementation for the __del__ built-in for frames.  
        
        Closes the Frame shelve database file and deletes the shelved Frame
        data if the frames.save parameter is False
        """
        self.frame_shelf.close()
        if not self.save:
            if os.path.isdir(self.path):
                shutil.rmtree(self.path)
        return

    def append(self,item):
        """ Append a Frame.
    
        The new Frame is added to the shelve database with a key 'self.nframes', 
        which is then incremented (self.nframe+=1).

        Args:                
            item (obj:Frame): The instance of a Frame object to be appended.
                           
        Raises:        
            TypeError: If item is not an instance of Frame.
        """
        if isinstance(item, Frame):
            self.frame_shelf[str(self.nframes)] = item
            self.frame_shelf.sync()
            self.nframes+=1
        else:
            raise TypeError(self._type_error)
        return

    def __getitem__(self,key):
        """ Get a copy of Frame from the database at key    

        The frames object is indexed with an integer key.

        Args:
            key (int): The index of the Frame object being called.
                       
        Returns:
            
        obj:Frame : This is an instance of the Frame object stored at index key (pulled from the Shelf database).
        """
        while key < 0:
            key += self.nframes
        while key > self.nframes:
            key = self.nframes-1
        
        return self.frame_shelf[str(key)]

    def __setitem__(self,key,item):
        """ Set the Frame at key in the Shelf database   

        Args:   
            key (int): The index where the input Frame should be stored.
            item (obj:Frame): This is an instance of a Frame object to be stored at index key.               
            
        Raises:
            TypeError : If the input item is not an instance of the Frame object
            
        """
        if not isinstance(item, Frame):
            raise TypeError(self._type_error)
            return
        if key < 0:
            key+=self.nframes
        elif key >= self.nframes:
            key = self.nframes
            self.nframes+=1
        self.frame_shelf[str(key)]=item
        self.frame_shelf.sync()
        return 

    def __len__(self):
        """ Returns the number of Frame objects being stored.
        
        Returns:
            int: The number of Frame objects this instance of frames.
        """
        return self.nframes

    def __iadd__(self,item):
        """ Use += operator to append Frame objects
        Appends by calling the append() function.    
        Args:
            item (obj:Frame): An instance of Frame to be appended.               
                
        """
        self.append(item)
        return self
            
# the multiprocessor parallelized functions that get copies of this object
# still return:
#    Exception OSError: OSError(2, 'No such file or directory') in  ignored
# I'm not sure why, but it is marked as ignored and it doesn't seem to cause any problems with the Frame shelve
# database file.
class ParFrames:
    """ Read-Only version of frames object  
    This class is effectively used to generate read-only copies of the frames class, which can be passed 
    to functions that do parallelized computations over the number of frames.  Unlike frames par_frames 
    does not create a new Shelf database. It must be passed an existing Shelf of Frame objects. par_frames
    also does not have any functions defined to modify the Shelf (add/remove Frame objects). This is to
    avoid conflicts from multiple processor accesses to the Shelf database.  
            
    """
    # fs_name does not actually get used, so it should probably be removed.
    def __init__(self, nframes, fs_name, frame_shelve):
        """ Initialize the par_frames object  

        Args:
        nframes (int): The number of Frames stored in the shelve database.
        fs_name (string): The base name (prefix) of the shelve database files.
        frame_shelve (shelve.Shelf): The Shelf object containing the Frame objects. 

        """
        self.nframes = nframes            
        self.fs_name = fs_name 
        #print "par_frames instance"
        #print "self.nframes ",self.nframes
        #print "self.fs_name ",self.fs_name
        #self.frame_shelf = shelve.open(self.fs_name,flag="r", protocol=2)
        self.frame_shelf = frame_shelve
        return
    
    def __getitem__(self,key):
        """ Get a copy of the Frame object stored at key.   
        
        Args:
            key (int): The index of the Frame object being called.
                       
        Returns:
            
        obj:Frame : This is an instance of the Frame object stored at index key (pulled from the shelve database)

        """
        while key < 0:
            key += self.nframes
        while key > self.nframes:
            key = self.nframes-1
        
        return self.frame_shelf[str(key)]

    def __len__(self):
        """ Returns the number of Frame objects being stored.
        
        Returns:
            int: The number of Frame objects this instance of frames.
        """
        return self.nframes

            
# leaflet object    
class Leaflet:
    """ Create a bilayer Leaflet representation.   
    This class object is used to group lipids together according to their bilayer leaflet. It is primarily meant to 
    store the indices of LipidCOMs as they are in a Frame.lipidcom list. This class also 
    creates sub-groups within the Leaflet based on the LipidCOM.type using LipidGroup objects. Instances of Leaflet
    are created by the MemSys class.
                
    """
    def __init__(self, name):
        """Initializes an instance of a Leaflet object.    
        
        Args:
            name (str): The name of the bilayer leaflet being initialized ('upper' and 'lower' are used by the MemSys class).
        
        Attributes:
            name (str): The name of the Leaflet (e.g. 'upper' or 'lower').
            members (list of int): A list containing the integer indices associated with the LipidCOM objects within
                a Frame that are assigned to the Leaflet instance.
            groups (list of obj:LipidGroup): A list of the LipidGroup objects (uniquely named) that are created by the Leaflet instance
                as new members are added.
            group_dict (dict): A dictionary keyed according to the names of the LipidGroup objects created, which stores the 
                corresponding index of that LipidGroup in self.groups. 
  
        """
        #the name of the leaflet - e.g. 'upper' or 'lower'
        self.name = name
        #initialize a list to store the indices of lipids assigned to this leaflet
        self.members = []
        #initialize a list to hold the LipidGroup objects 
        self.groups = []
        #initialize a dictionary to store the self.groups index of LipidGroup objects
        self.group_dict = {}
        return

    def __str__(self):
        return '%s leaflet of a Center of mass trajectory with %s members and %s lipid groups' % (self.name, len(self.members), len(self.groups)) 

    def __repr__(self):
        return '%s leaflet of a Center of mass trajectory with %s members and %s lipid groups' % (self.name, len(self.members), len(self.groups))
 
    def __len__(self):
        """ Have len(Leaflet) return the number of lipids that have been added to the Leaflet instance.
        
        Returns:
            int: Number of lipids in the Leaflet.
        """
        return len(self.members)

    #consider changing var name of input 'resname' to something that doesn't conflict with LipidCOM.type  
    def add_member(self, index, resname):
        """ Add new lipids to the Leaflet.
           
        This function is meant to be used to add new lipids according to their Frame.lipidcom index
        to the Leaflet and to a LipidGroup according resname/type/name. 
        Args:  
            index (int): The index of the lipid being added to the Leaflet.
            resname (str): The resname (or LipidCOM.type) of the lipid being added.
        """
        if len(self.members) == 0:
            self.members.append([index, resname])
            self.groups.append(LipidGroup(resname))
            self.groups[0].add_member(index)
            self.group_dict.update({resname: 0})
        else:
            self.members.append([index, resname])
            addgroup = True
            group_ind = 0
            for rn in self.groups:
                if resname == rn.lg_name:
                    addgroup = False
                    break
                group_ind+=1
            if addgroup:
                self.groups.append(LipidGroup(resname))
                ng = len(self.groups)
                self.groups[ng-1].add_member(index)
                self.group_dict.update({resname: ng-1})
            else:
                self.groups[group_ind].add_member(index)
            
            #self.members=sorted(self.members,key=lambda self.members:self.members[1])  
        return

    def get_group_indices(self, group_name):
        """ Get the indices of lipids in the Leaflet belonging to a specific LipidGroup.   

        Args:
        group_name (string): The name of the LipidGroup pull LipidCOM indices from. 
            Passing the string 'all' will return indices of all the lipids assigned to
            the Leaflet instance. If the group_name is not recognised (i.e. is not in the group_dict)
            The function defaults to 'all'.
                       
        Returns:
            list of int: A list containing the integer indices of lipids in the Leaflet that
                belong to the specified LipidGroup.
        """
        indices = []
        if group_name == "all":
            for element in self.group_dict:
                gindex = self.group_dict[element]
                indices += self.groups[gindex].lg_members
        elif group_name in self.group_dict:
            gindex = self.group_dict[group_name]
            indices = self.groups[gindex].lg_members
        else:
            #unkwown group name- print warning and use the default "all"
            print "!! Warning - request for unknown Lipid Group \'",group_name,"\' from the ",self.name," leaflet"
            print "!! using the default \"all\""
            for element in self.group_dict:
                gindex = self.group_dict[element]
                indices += self.groups[gindex].lg_members

        return list(indices)

    def get_member_indices(self):
        """ Get the indices of all lipids (LipidCOM) in the Leaflet.   
        This member function Returns: the list of indices for the lipids grouped in the Leaflet instance.
                           
        Returns:
            list of int: A list of integer indices of the lipids associated with the Leaflet instance.
        """
        indices = []
        for element in self.members:
            indices.append(element[0])

        return list(indices)

    def has_group(self, group_name):
        """ Check if there is a LipidGroup with the specified name.  
            
        Args:
            group_name (str): The name to checked against the names of existing LipidGroup objects.
                           
        Returns:
            bool: True if there is a LipidGroup with name group_name, and False otherwise.
        """
        return group_name in self.group_dict.keys()

    def num_groups(self):
        """ Get the number of LipidGroups in the Leaflet.   
    
        Returns:
        int: The number of unique LipidGroups.
        """
        return len(self.groups)

    def get_group_names(self):
        """ Get the names of all the LipidGroup objects in the Leaflet  

        Returns:
        list of str: A list of the names of current LipidGroup objects.
        """
        return [group.lg_name for group in self.groups]

        
class LipidGroup:
    """ Object to group lipid indices by type/resname/name.   
        Instances of this object are created by the Leaflet class.
                       
    """
    def __init__(self, name):
        """ Initializes LipidGroup object.               

        Args:
            name (str): The name/type/resname of the lipids being grouped in this object. 

        Attributes:
            lg_members (list of int): A list to hold the indices of lipids added to this 
                this LipidGroup.
            lg_name (str): The name/type/resname of the lipids being grouped in this object.
        """
        #initialize a list to hold the member indices
        self.lg_members = []
        # the name of this lipid group
        self.lg_name = name
        return

    def add_member(self, new_mem):
        """ Add lipid index to to the LipidGroup.
   
         Args:
            new_mem (int): The index of the lipid being added to this LipidGroup. 
        """
        self.lg_members.append(new_mem)
        return

    def name(self):
        """ Get the name associated with this LipidGroup.   

        Returns:               
            str: The name of the lipid group (i.e. lg_name) 
        """
        return self.lg_name

def msd_frames(frames, fstart, fend, indices, refframe, plane):
    """ Compute the mean squared displacement for range of Frame objects in frames.

    This function allows the mean squared displacement (MSD) to be computed
    for a specified subset of the Frame objects in a frames (or par_frames) object.     
    This function was created to be called from the function MemSys.CalcMSD_parallel
    as a function to be passed to the multiprocessor threads.   
        
    Args:           
        frames (obj:frames or obj:par_frames): The object containing all the Frames of the trajectory.
        fstart (int): The index of the first frame to start the analysis on.
        fend (int): The index of the last frame to analyze.
        indices (list of int): List of integer indices of the LipidCOMs to include in the computation.
        refframe (int): The index of the frame that is to be taken as the reference for the MSD computation.
        plane (list of int): The list of the indices corresponding to the coordinate planes (x: 0,y 1,z :2) 
            to be included in the MSD computation.   
        
    Returns:
        numpy.array: This is a nx2 numpy array (of floats) containing the 
            results of the MSD computation for the specified frames.
            msd_results[i,0] = simulation time for frame f = i + fstart.
            msd_results[i,1] = the configurational average MSD over the specified LipidCOMs for frame f = i + fstart.
            for i in range( (fend-fstart)+1 ).
    """
    #initialize an array to hold the ouptut
    nfc = fend - fstart + 1
    output = np.zeros((nfc, 2))
    # number of lipids in the selection
    n_com = len(indices)
    #initialize a running stats object to do the configuration averaging
    drs_stat = RunningStats()
    # initialize an np array to hold coordinates for the selection
    # at the reference frame
    com_ref = np.zeros((n_com, 2))
    ref_frame = frames[refframe]
    count=0
    # get the coordinates
    for i in indices:
        com_i = ref_frame.lipidcom[i].com_unwrap[plane]
        com_ref[count]=com_i[:]
        count+=1
    time_ref = ref_frame.time
    #print "nframes ",len(frames)
    #print "process; fstart ",fstart," fend ",fend
    #print "process; loop range "
    #print range(fstart,(fend+1))
    # now begin loop over the frames for this process        
    for f in range(fstart, (fend+1)):
        # get the current frame
        curr_frame = frames[f]
        # get the coordinates for the selection at this frame
        com_curr = np.zeros((n_com, 2))
        count=0
        for i in indices:
            com_i = curr_frame.lipidcom[i].com_unwrap[plane]
            com_curr[count]=com_i[:]
            count+=1
        #current time
        tc = curr_frame.time
        dr = com_curr - com_ref
        drs = dr*dr
        #loop over the selections for this frame
        for    val in drs:
            drs_curr = val[:]    
            drs_mag = drs_curr.sum()
            drs_stat.push(drs_mag)
        #get the msd for the current selection
        msdcurr = drs_stat.mean()
        devcurr = drs_stat.deviation()
        drs_stat.Reset()
        findex = f-fstart
        output[findex, 0]=tc
        output[findex, 1]=msdcurr
#        output[findex,2]=devcurr
#        dt = tc - time_ref
#        DiffCon = 0.0
#        if f != 0:
#            DiffCon = msdcurr/(4.0*dt)
#        output[findex,3]=DiffCon
    #    print "msdcurr ",msdcurr," DiffCon ",DiffCon
    return output

#function to compute the thickness of the membrane (in the normal direction). The algorithm is based on  
# the GridMAT-MD bilayer thickness calculation (except without the gridding procedure) 
def thickness_frames(frames, fstart, fend, leaflets, nlipids, plane, norm):
    """ Compute the bilayer thickness for a range of Frame objects in frames.
    Computes the thickness of the bilayer (along the normal direction). The algorithm is based on  
    the GridMAT-MD bilayer thickness calculation.  
    This function was created to be called used in MemSys.CalcThickness_parallel
    as a function to be passed to the multiprocessor threads.   
    
    Args:   
        frames (obj:frames or obj:par_frames): Object containing all the Frame objects of the trajectory.
        fstart (int): The index of the first frame to start the analysis on.
        fend (int): The index of the last frame to analyze.
        leaflets (dict of Leaflet): A dict containing the Leaflet instances used to define the upper and lower
            bilayer leaflets for this calculation. This dict should contain the two keys, 'upper' and 'lower', corresponding
            to instances of the Leaflet class.
        nlipids (int): The total number of LipidCOMs (or lipids) in the two Leaflet instances.
        plane (list of int): A list of the integer indices corresponding to the bilayer lateral 
            coordinate planes (0 for x,1 for y, and 2 for z)
        norm (int): The integer index corresponding to the bilayer normal coordinate plane 
            (0 for x,1 for y, or 2 for z).
        
    Returns:
        tuple of numpy.array: This is a two element tuple containing numpy arrays of the computation results.
            tuple[0] => thickness: A nx3 numpy array containing the 
            results of the thickness computation for the specified frames. Specifically:
                thickness[i,0] = simulation time for frame f = i + fstart.
                thickness[i,1] = the configurational average thickness for frame f = i + fstart.
                thickness[i,2] = the standard deviation of the configurational average thickness for frame f = i + fstart.
                For i in xrange( (fend-fstar) + 1).
            tuple[1] => thickness_map: A nxNx6 numpy array containing the thickness data that can be 
            used to generate a 3d thickness map/plot. Specifically:
                thickness[i,j,0] = simulation time for frame f = i + fstart and lipid j.
                thickness[i,j,1] = the average x position for lipid j and its cross leaflet partner at frame f = i + fstart.
                thickness[i,j,2] = the average y position for lipid j and its cross leaflet partner at frame f = i + fstart.
                thickness[i,j,3] = the lower z position for lipid j and its cross leaflet partner at frame f = i + fstart.
                thickness[i,j,4] = the upper z position for lipid j and its cross leaflet partner at frame f = i + fstart.
                thickness[i,j,5] = the difference between the upper and lower z positions 
                for lipid j and its cross leaflet partner at frame f = i + fstart.
                For i in xrange((fend-fstar) + 1) and For j in xrange(nlipids).
    """
    #upper_match = []
    #lower_match = []
    xi = plane[0]
    yi = plane[1]
    zi = norm
    comcup = np.zeros(3)
    comclo = np.zeros(3)
    dcom = np.zeros(3)
    nfc = fend - fstart + 1
    nlc = nlipids
    zdists = np.zeros((nfc, nlc, 1))
    zmaps = np.zeros((nfc, nlc, 6))
    #dcoms = np.zeros(3)
    f=0
    times = np.zeros(nfc)
    
    for    f in range(fstart,(fend+1)):
        n=0
        fr = frames[f]
        boxc = fr.box
        boxc_xh = boxc[xi]/2.0
        boxc_yh = boxc[yi]/2.0
        dt = fr.time
        findex = f-fstart
        times[findex]=dt
        for memu in leaflets['upper'].members:
            idu = memu[0]
            comcup = fr.lipidcom[idu].com
            distxy = 10000.0
            distz = 0.0
            mindex = 0
            zlom = 0.0
            zhim = 0.0
            xavgm = 0.0
            yavgm = 0.0
            for meml in leaflets['lower'].members:
                idl = meml[0]
                comclo = fr.lipidcom[idl].com
                dcom = comcup-comclo
                dx = dcom[xi]
                dy = dcom[yi]
                dz = dcom[zi]
                #Minimum image -- coordinates must be pre-wrapped 
                if np.absolute(dx) > boxc_xh:
                    dx = boxc[xi] - np.absolute(comcup[xi]-boxc_xh) - np.absolute(comclo[xi]-boxc_xh)
                if np.absolute(dy) > boxc_yh:
                    dy = boxc[yi] - np.absolute(comcup[yi]-boxc_yh) - np.absolute(comclo[yi]-boxc_yh)
                rxy = np.sqrt(dx**2+dy**2)
                #get 4d map values
                comavg = (comcup+comclo)/2.0
                xavg = comavg[xi]
                yavg = comavg[yi]
                zlo = comclo[zi]
                zhi = comcup[zi]
                if    rxy<distxy:
                    distxy=rxy
                    distz = np.absolute(dz)
                    mindex=meml
                    xavgm = xavg
                    yavgm = yavg
                    zlom = zlo
                    zhim = zhi
                    
            #upper_match.append([mindex,distz])
            #print "n ",n," xvg ", xavgm," yvg ", yavgm
            
            zdists[findex,n]=distz
            #maps
            zmaps[findex,n,0]=dt
            zmaps[findex,n,1]=xavgm
            zmaps[findex,n,2]=yavgm
            zmaps[findex,n,3]=zlom
            zmaps[findex,n,4]=zhim
            zmaps[findex,n,5]=distz
            
            n+=1
        for meml in leaflets['lower'].members:
            idl = meml[0]
            comclo = fr.lipidcom[idl].com
            distxy = 10000.0
            distz = 0.0
            mindex = 0
            zlom = 0.0
            zhim = 0.0
            xavgm = 0.0
            yavgm = 0.0
            for memu in leaflets['upper'].members:
                idu = memu[0]
                comcup = fr.lipidcom[idu].com
                dcom = comclo-comcup
                dx = dcom[xi]
                dy = dcom[yi]
                dz = dcom[zi]
                #Minimum image -- coordinates must be pre-wrapped 
                if np.absolute(dx) > boxc_xh:
                    dx = boxc[xi] - np.absolute(comclo[xi]-boxc_xh) - np.absolute(comcup[xi]-boxc_xh)
                if np.absolute(dy) > boxc_yh:
                    dy = boxc[yi] - np.absolute(comclo[yi]-boxc_yh) - np.absolute(comcup[yi]-boxc_yh)
                rxy = np.sqrt(dx**2+dy**2)
                #get 4d map values
                comavg = (comcup+comclo)/2.0
                xavg = comavg[xi]
                yavg = comavg[yi]
                zlo = comclo[zi]
                zhi = comcup[zi]
                if    rxy<distxy:
                    distxy=rxy
                    distz = np.absolute(dz)
                    mindex=meml
                    xavgm = xavg
                    yavgm = yavg
                    zlom = zlo
                    zhim = zhi
            #upper_match.append([mindex,distz])
            #print "n ",n," xvg ", xavgm," yvg ", yavgm
            zdists[findex,n]=distz
            #maps
            zmaps[findex,n,0]=dt
            zmaps[findex,n,1]=xavgm
            zmaps[findex,n,2]=yavgm
            zmaps[findex,n,3]=zlom
            zmaps[findex,n,4]=zhim
            zmaps[findex,n,5]=distz
            n+=1
        
        #break
    zavgs = np.zeros((nfc, 3))
    zdtstat = RunningStats()    
    for fr in xrange(nfc):
        currtime = times[fr]
        dt = currtime 
        curr = zdists[fr,:]
        zavgcurr = curr.mean()            
        zdevcurr = curr.std()
#        zdtstat.push(zavgcurr)
#        zdtcurr = zdtstat.mean()
#        zdtdcurr = zdtstat.deviation()
        zavgs[fr,0]=dt
        zavgs[fr,1]=zavgcurr
        zavgs[fr,2]=zdevcurr
#        zavgs[fr,3]=zdtcurr
#        zavgs[fr,4]=zdtdcurr
    out = (zavgs, zmaps)
    return out
    #return zavgs
    #return zmaps
        

## this is the main class - the Center of mass trajectory (MemSys) object
# The lipids are assumed to be separate residues, so the input selection (mem_sel)
# is split according to residues with each (full) residue assigned as a LipidCOM. However,
# it would nice to make it possible to something like a list of selections. That way 
# you could for example use selections of the lipid headgroups as the LipidCOMs representation.
class COMTraj:
    """ This is the main class object.
    An instance of this class reads in the trajectory and a selection (both MDAnalysis objects) 
    and creates/reduces the lipids to center of mass (COM) representations; i.e. a center of mass trajectory. There are several member
    functions to perform various types of analyses based on the COM representations.    
                       
    """
    # pass the mda anaylis trajectory object and a selection with the membrane (i.e. w/o water and ions)
    # optional - specify the plane that the membrane is in - default is xy with normal in z
    def __init__(self, mda_traj, mem_sel, plane="xy",fstart=0, fend=-1, fskip=1,frame_path='Default',frame_save=False,nprocs=1):
        """ COMTraj initialization.    

        Each lipid in the system is assumed to be its own residue. It is also assumed that the coordinates 
        in the MDAnalysis trajectory are wrapped, so when this object is initialized it computes COMs using
        the wrapped coordinates. It then unwraps the raw coordinates and removes the system COM motion before 
        recomputing another set of lipid COMs (now based on the unwrapped coordinates).
        Both representations are stored and are used in the various analyses. 
         

        Args:
            mda_traj (MDAnalysi.Universe.trajectory): The MDAnalysis trajectory object containg the 
                the molecular dynamics trajectory data for the system.
            mem_sel (MDAnalysis.AtomGroup): The MDAnalysis atom selection (AtomGroup) containing
                the bilayer lipids.
            plane (str, optional): Defines the lateral plane of the bilayer. The default is 'xy'. 
                The other options are 'yz' or 'xz', or their equivalent 'zy' or 'zx'.
            fstart (int, optional): The index of the initial frame to use from the MDAnalysis trajectory.
                The default is 0.
            fstart (int, optional): The index of the final frame to use from the MDAnalysis trajectory.
                The default is -1 (i.e. the last frame).
            fskip (int, optional): Use every fskip frame in mda_traj. 
                The default is 1 (i.e. use all frames in mda_traj).
            frame_path (str, optional): Path in which to store the Shelf database of the Frame objects
                created for the COM representations. The default is 'Default' which uses the '/tmp/' dir.
            frame_save (bool, optional): Preserve the Shelf database files of Frame objects for this
                MemSys after deletion. The default is False (i.e. the Shelf database files are deleted). 
                If True then the Shelf database files are not removed when this MemSys object is deleted.
            nprocs (int, optional): Define the number of processors to use when unwrapping the coordinates.
                The default is 1 (i.e. no paralellization). Greater than 1 uses a parallelized version
                of the unwrap function when uwrapping coordinates. nprocs should not exceed the 
                the number of cores on your machine.   
              

        Atrributes:
           plane (list of int): A two element list of the integer indices corresponding to the bilayer lateral 
               coordinate planes (0 for x,1 for y, and 2 for z). Determined by the value of Args:plane. 
           norm (int): The integer index corresponding to the bilayer normal coordinate plane 
               (0 for x,1 for y, or 2 for z). Determined by the value of Args:plane.
           leaflets (dict of Leaflet): A dict containing the Leaflet instances used to define the upper and lower
               bilayer leaflets for this system. This dict will contain the two keys, 'upper' and 'lower', upper
               and lower Leaflet objects respectively.
           com_leaflet (list of str): A list of leaflet name strings that maps the lipid index
               to the Leaflet ('upper' or 'lower') it has been assigned to. 
               E.g. 
               >>> com_leaf[0]
               >>> 'upper'
               means lipid index 0 was assigned to the 'upper' leaflet.  
           nlipids (int): The number of lipids in the system. Each MDAnalysis residue in Args:mem_sel is 
               treated as a lipid. 
           clusters (list of list of list int): This is initialized as an empty list. It is used by the function
               CheckClustering() as a persistent container to store the list of clusters for each frame. clusters is
                overwritten after each call to CheckClustering().
           frame (obj:frames): An instance of the Frame container object (frames) to hold the COM trajectory data
               for this system.
           nframes (int): The total number of frames processed and stored for the MemSys instance.  
            
            
        """
        #defaults - xy plane with z normal
        ii=0
        jj=1
        kk=2    
        if    plane=="yz" or plane=="zy":
            ii=1
            jj=2
            kk=0
        if    plane=="xz" or plane=="zx":
            ii=0
            jj=2
            kk=1
        #parallelize loading -- currently just applies to unwrapping
        parallel=False
        if nprocs>1:
            parallel=True
        #store the indices of the plane directions
        self.plane = [ii, jj]
        # store the index of the normal direction
        self.norm = kk 
            
        #initialize leaflet objects
        self.leaflets = {'upper':Leaflet('upper'),'lower':Leaflet('lower')}
        self.com_leaflet = []
    
        
        #get the number of lipids (residues)
        self.nlipids=mem_sel.n_residues
        #initialize an empty cluster list - used to store the clusters built in the last call of 'CheckClustering'
        self.clusters = [] # after 'CheckClustering' is called, the outersize len(self.clusters) should equal self.nframes
        #initialize empty frame list
        #self.frame=[]
        if frame_path == 'Default':
            frame_path = '/tmp/'
        self.frame = FrameShelve(prefix=frame_path,save=frame_save)
        #adjust for slicing index
        len_mda_traj = len(mda_traj)
        #adjust for negative indexing
        while fstart < 0:
            fstart+=len_mda_traj 
        while fend < 0:
            fend+=len_mda_traj 

        #adjust endpoint for slicing
        fend+=1
         
        #loop over the frames
        f=0
        #fskip-=1
        # you can slice the MDAnalysis trajectory in a loop
        for frame in mda_traj[fstart:fend:fskip]:
            mdframe = frame.frame
            print "doing full trajectory frame ", mdframe, ", COM trajectory frame ",f 
            #add the frame object for this frame
            cframe = Frame(self.nlipids)
            # set the box dimensions and the time for this frame
            dimensions = frame.dimensions[0:3]
            if (dimensions == 0).any() and f>0:
                dimensions = self.frame[f-1].box[:]
            cframe.set_box(dimensions)
            cframe.set_time(frame.time)
            print "time ",mda_traj.time
            cframe.number = f
            cframe.mdnumber = mdframe
            # loop over the residues (lipids) and get the centers of mass
            r=0            
            for res in mem_sel.residues:
                cframe.lipidcom[r].extract(res, box=dimensions)
                cframe.lipidcom[r].mass = res.total_mass()
                r+=1
            #append the frame
            self.frame.append(cframe)
            f+=1
        #get the number of frames from the trajectory
        self.nframes = f
        #now we need to unwrap the coordinates
        natoms = len(mem_sel)
        oldcoord = np.zeros((natoms,3))
        currcoord = np.zeros((natoms,3))
        wrapcoord = np.zeros((natoms,3))
        index = mem_sel.indices

        firstframe = True
        # loop over the trajectory again to get unwrapped coordinates
        # unwrap the raw residue coordinates - then get the COMs
        f=0
        for frame in mda_traj[fstart:fend:fskip]:    
            #first we unwrapp
            #print "unwrapping frame ",frame.frame
            currcoord = frame._pos[index]
            if firstframe:
                oldcoord = np.copy(currcoord)
                firstframe = False
            else:
                abc = frame.dimensions[0:3]
                if parallel:
                    wrapcoord = mda_wrap_coordinates_parallel(abc, currcoord, oldcoord,nprocs=nprocs)
                else:
                    wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
                frame._pos[index] = wrapcoord[:]
                oldcoord = np.copy(wrapcoord)
            #now we need to adjust for the center of mass motion of the membrane -- for simplicity set all frames to (0,0,0)
            # to remove center of mass motion of the membrane
            mem_com = mem_sel.center_of_mass()
            frame._pos[index] -= mem_com
            r=0    
            cframe = self.frame[f]        
            for res in mem_sel.residues:
                cframe.lipidcom[r].extract(res, unwrap=True)
                r+=1
            self.frame[f]=cframe
            f+=1            

        # now we can assign the lipids to the leaflets 
        # NOTE: Lipids are only assigned to leaflets once based on the 
        #       first frame of the trajectory
        
        #first- compute the average position along the normal direction
        zstat = RunningStats()
        for lipcom in self.frame[0].lipidcom:
            zstat.push(lipcom.com_unwrap[self.norm])
        zavg = zstat.mean()
        # now loop over the lipids
        l = 0        
        for lipcom in self.frame[0].lipidcom:
            pos = ""
            # decide which leaflet
            if lipcom.com_unwrap[self.norm]>zavg:
                pos = 'upper'                
            elif lipcom.com_unwrap[self.norm]<zavg:
                pos = 'lower'
            #add to the chosen leaflet
            self.com_leaflet.append(pos)
            self.leaflets[pos].add_member(l, lipcom.type)    
            l+=1
        #complete
        return
    
    def __str__(self):
        return 'Center of mass trajectory with %s frames and %s lipids/components' % (self.nframes, self.nlipids) 
    def __repr__(self):
        return 'Center of mass trajectory with %s frames and %s lipids/components' % (self.nframes, self.nlipids)

    def number_of_unique_groups(self):
        """ Get the number of uniquely named LipidGroups within both Leaflet objects.   
    
        Returns:
        int: The number of uniquely named LipidGroup objects.
        """
        resnames = []
        for leaflet in self.leaflets:
            for group in leaflet.groups:
                gname = group.name()
                if gname not in resnames:
                    resnames.append(gname)
        return len(resnames)
    #def LeafletCOM(leaflet_name,frame_num):
        

    # function to compute the mean squared displace (msd) along with the diffusion constant of a group 
    # Possibly add functionality to specify the range of trajectory frames to include int computaton.
    def calc_msd(self, leaflet="both",group="all"):
        """ Compute the configurational average mean squared displacement for select lipids in a select leaflet(s).

        This function allows the mean squared displacement (MSD) to be computed
        for a specified leaflet ('upper', 'lower', or 'both') and for a specified LipidGroup within the 
        chosen Leaflet. This calculation is over the whole trajectory and assumes the first frame is the 
        reference frame for computing the displacement.  
            
        Args:           
            leaflet (str): A string designating which Leaflet to include in the computation.
            group (str): A string with the name of a specific LipidGroup to in the computation. 
            
        Returns:
            numpy.array: This is a nframesx2 numpy array (of floats) containing the 
                results of the MSD computation across all frames.
                msd_results[i,0] = simulation time for frame i.
                msd_results[i,1] = the configurational average MSD over the specified LipidCOMs for frame f = i.
                For i in range( nframes ).
        """
        # initialize a list to hold the indices of LipidCOMs to be included in this computaton.
        indices = []
        #diffusion dimension - assume lateral so, dim=2 -- Although this could made an optional parameter.
        dim=2
        # parse the leaflet and group inputs
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        n_com = len(indices)
        #store the coordinates of the selected LipidCOMs in a single numpy array
        selcoords = np.zeros((self.nframes,n_com,2))
        
        for f in xrange(self.nframes): 
            count=0
            for i in indices:
                com_curr = self.frame[f].lipidcom[i].com_unwrap[self.plane]
                selcoords[f,count]=com_curr[:]
                count+=1
        
        #initialize a numpy array to hold the msd for the selection        
        msd = np.zeros((self.nframes, 2))
        #initialize a running stats object to do the averaging
        drs_stat = RunningStats()
        #initialize a running stats object for the diffusion constant (frame/time average)
        diff_stat = RunningStats()
        #running stats object for time averaging
        msd_stat = RunningStats()
        #loop over the frames starting at index 1
        #print comlist
        #print len(comlist)
        coml0 = selcoords[0,:,:]
        t0 = self.frame[0].time
        #print coml0
        for i in xrange(1, self.nframes):
            # get the current com frame list
            tc = self.frame[i].time
            dt = tc
            comlcurr = selcoords[i,:,:]
            dr = comlcurr - coml0
            drs = dr*dr
            #loop over the selections for this frame
            for    val in drs:
                drs_curr = val[:]    
                drs_mag = drs_curr.sum()
                drs_stat.push(drs_mag)
            #get the msd for the current selection
            msdcurr = drs_stat.mean()
            devcurr = drs_stat.deviation()
            drs_stat.Reset()
            msd_stat.push(msdcurr)
            msd_tavg = msd_stat.mean()
            msd_dev = msd_stat.deviation()            
            #dt = times[i]-times[0]
            DiffCon = msd_tavg/(2.0*dim*dt)
            diff_stat.push(DiffCon)
            #print "msdcurr ",msdcurr
            #push to the msd array
            
            msd[i,0]=dt
            msd[i,1]=msdcurr
#            msd[i,2]=msd_tavg
#            msd[i,3]=msd_dev
#            msd[i,4]=DiffCon
#            msd[i,5]=diff_stat.mean()
#            msd[i,6]=diff_stat.deviation()
        #return msd array
        return msd 

    #function to compute the thickness of the membrane (in the normal direction). The algorithm is based on  
    # the GridMAT-MD bilayer thickness calculation (except without the gridding procedure) 
    def calc_membrane_thickness(self):
        """ Compute the bilayer thickness across the trajectory.
        Computes the thickness of the bilayer (along the normal direction). The algorithm is based on  
        the GridMAT-MD bilayer thickness calculation.     
            
        Returns:
            tuple of numpy.array: This is a two element tuple containing numpy arrays of the computation results.
                tuple[0] => thickness: A nx5 numpy array containing the 
                results of the thickness computation for the specified frames. Specifically:
                    thickness[i,0] = simulation time for frame f = i + fstart.
                    thickness[i,1] = the configurational average thickness for frame f = i + fstart.
                    thickness[i,2] = the standard deviation of the configurational average thickness for frame f = i + fstart.
                    thickness[i,3] = the running time average of the configurational average thickness for frame f = i + fstart.
                    thickness[i,4] = the running standard deviation of the time averaged configurational average thickness for frame f = i + fstart.
                    For i in xrange( (fend-fstar) + 1).
                tuple[1] => thickness_map: A nxNx6 numpy array containing the thickness data that can be 
                used to generate a 3d thickness map/plot. Specifically:
                    thickness[i,j,0] = simulation time for frame f = i + fstart and lipid j.
                    thickness[i,j,1] = the average x position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,2] = the average y position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,3] = the lower z position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,4] = the upper z position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,5] = the difference between the upper and lower z positions 
                    for lipid j and its cross leaflet partner at frame f = i + fstart.
                    For i in xrange((fend-fstar) + 1) and For j in xrange(nlipids).
        """
        #upper_match = []
        #lower_match = []
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        comcup = np.zeros(3)
        comclo = np.zeros(3)
        dcom = np.zeros(3)
        zdists = np.zeros((self.nframes, self.nlipids, 1))
        zmaps = np.zeros((self.nframes, self.nlipids, 6))
        #dcoms = np.zeros(3)
        f=0
        
        for f in xrange(self.nframes): 
            n=0
            fr = self.frame[f]
            boxc = fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            dt = fr.time
            for memu in self.leaflets['upper'].members:
                idu = memu[0]
                comcup = fr.lipidcom[idu].com
                distxy = 10000.0
                distz = 0.0
                mindex = 0
                zlom = 0.0
                zhim = 0.0
                xavgm = 0.0
                yavgm = 0.0
                for meml in self.leaflets['lower'].members:
                    idl = meml[0]
                    comclo = fr.lipidcom[idl].com
                    dcom = comcup-comclo
                    dx = dcom[xi]
                    dy = dcom[yi]
                    dz = dcom[zi]
                    #Minimum image -- coordinates must be pre-wrapped 
                    if np.absolute(dx) > boxc_xh:
                        dx = boxc[xi] - np.absolute(comcup[xi]-boxc_xh) - np.absolute(comclo[xi]-boxc_xh)
                    if np.absolute(dy) > boxc_yh:
                        dy = boxc[yi] - np.absolute(comcup[yi]-boxc_yh) - np.absolute(comclo[yi]-boxc_yh)
                    rxy = np.sqrt(dx**2+dy**2)
                    #get 4d map values
                    comavg = (comcup+comclo)/2.0
                    xavg = comavg[xi]
                    yavg = comavg[yi]
                    zlo = comclo[zi]
                    zhi = comcup[zi]
                    if    rxy<distxy:
                        distxy=rxy
                        distz = np.absolute(dz)
                        mindex=meml
                        xavgm = xavg
                        yavgm = yavg
                        zlom = zlo
                        zhim = zhi
                        
                #upper_match.append([mindex,distz])
                #print "n ",n," xvg ", xavgm," yvg ", yavgm
                zdists[f,n]=distz
                #maps
                zmaps[f,n,0]=dt
                zmaps[f,n,1]=xavgm
                zmaps[f,n,2]=yavgm
                zmaps[f,n,3]=zlom
                zmaps[f,n,4]=zhim
                zmaps[f,n,5]=distz
                
                n+=1
            for meml in self.leaflets['lower'].members:
                idl = meml[0]
                comclo = fr.lipidcom[idl].com
                distxy = 10000.0
                distz = 0.0
                mindex = 0
                zlom = 0.0
                zhim = 0.0
                xavgm = 0.0
                yavgm = 0.0
                for memu in self.leaflets['upper'].members:
                    idu = memu[0]
                    comcup = fr.lipidcom[idu].com
                    dcom = comclo-comcup
                    dx = dcom[xi]
                    dy = dcom[yi]
                    dz = dcom[zi]
                    #Minimum image -- coordinates must be pre-wrapped 
                    if np.absolute(dx) > boxc_xh:
                        dx = boxc[xi] - np.absolute(comclo[xi]-boxc_xh) - np.absolute(comcup[xi]-boxc_xh)
                    if np.absolute(dy) > boxc_yh:
                        dy = boxc[yi] - np.absolute(comclo[yi]-boxc_yh) - np.absolute(comcup[yi]-boxc_yh)
                    rxy = np.sqrt(dx**2+dy**2)
                    #get 4d map values
                    comavg = (comcup+comclo)/2.0
                    xavg = comavg[xi]
                    yavg = comavg[yi]
                    zlo = comclo[zi]
                    zhi = comcup[zi]
                    if    rxy<distxy:
                        distxy=rxy
                        distz = np.absolute(dz)
                        mindex=meml
                        xavgm = xavg
                        yavgm = yavg
                        zlom = zlo
                        zhim = zhi
                #upper_match.append([mindex,distz])
                #print "n ",n," xvg ", xavgm," yvg ", yavgm
                zdists[f,n]=distz
                #maps
                zmaps[f,n,0]=dt
                zmaps[f,n,1]=xavgm
                zmaps[f,n,2]=yavgm
                zmaps[f,n,3]=zlom
                zmaps[f,n,4]=zhim
                zmaps[f,n,5]=distz
                n+=1
            
            #break
        zavgs = np.zeros((self.nframes, 5))
        zdtstat = RunningStats()    
        for fr in xrange(self.nframes):
            currtime = self.frame[fr].time
            dt = currtime 
            curr = zdists[fr,:]
            zavgcurr = curr.mean()            
            zdevcurr = curr.std()
            zdtstat.push(zavgcurr)
            zdtcurr = zdtstat.mean()
            zdtdcurr = zdtstat.deviation()
            zavgs[fr,0]=dt
            zavgs[fr,1]=zavgcurr   
            zavgs[fr,2]=zdevcurr
            zavgs[fr,3]=zdtcurr
            zavgs[fr,4]=zdtdcurr

        return (zavgs, zmaps)
        #return zmaps

    # a simple cluster/chain analysis routine
    def check_clustering(self, leaflet="both",group="all", dist=15.0):
        """ Determine physical cluster for select lipid COMs in a select leaflet(s).

        This function determines physical clusters based on Cartesian distance criteria
        in a specified leaflet ('upper', 'lower', or 'both') and for a specified LipidGroup within the 
        chosen Leaflet. The routine determines the physical clusters for each frame in the trajectory
        and computes properties of these clusters, e.g. number of clusters, which are time averaged over the trajectory. 
        The information on the clusters is stored in MemSys.clusters which can then be used by other functions 
        (like MemSys.ExportClustersForPlotting).    
            
        Args:           
            leaflet (str, optional): A string designating which Leaflet to include in the computation.
            group (str, optional): A string with the name of a specific LipidGroup to in the computation.
            dist (float, optional): The distance cutoff for determining clusters.
            
        Returns:
            numpy.array: This is a nframesx5 numpy array (of floats) containing the 
                results of this computation across all frames.
                cluster_results[i,0] = simulation time for frame i.
                cluster_results[i,1] = the current running time averaged number of clusters at frame f = i.
                cluster_results[i,2] = the current standard deviation of the running time averaged number of clusters at frame f = i.
                cluster_results[i,3] = the current running time averaged configurational average lipids per cluster at frame f = i.
                cluster_results[i,4] = the current standard deviation of the running time averaged 
                    configurational average lipids per cluster at frame f = i.
                For i in range( nframes ).
        """
        indices = []
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        #reset the system cluster list
        self.clusters = []
        # numpy array to store output for return
        outdata = np.zeros((self.nframes,5))
        #stats objects - time averages
        ncstat = RunningStats() #number of clusters
        asstat = RunningStats() # average cluster size
        misstat = RunningStats() # minimum cluster size
        masstat = RunningStats() # maximum cluster size                
        #loop over frames        
        
        for f in xrange(self.nframes):
            fr = self.frame[f]
            ctime = fr.time
            clusters = []
#            masterlistf = []
#            masterlistf += masterlist
            #rebuild the master list each frame
            masterlistf = list()
            for i in indices:
                masterlistf.append([i, False])
#            print "master ",masterlistf
            boxc=fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            #print boxc
            clustind = 0
            neighborlist = []
            while len(masterlistf)>0:
                #print "master ",masterlistf
                start = masterlistf[0][0]
                masterlistf[0][1]=True
            #    print 
                # reset the neighborlist
                neighborlist = []
                #seed the neighborlist with the start
                neighborlist.append(start)
                #now loop over the neighborlist and build neighbors and neighbors of neigbors for this cluster
                i=0
                while i < len(neighborlist):
                    ele = neighborlist[i]
                    startn = ele
                    coms = fr.lipidcom[startn].com        
                    #get neighbors of the start
                    #mindex=0
                    for j in xrange(len(masterlistf)):
                    #for elem in masterlistf:
                        elem = masterlistf[j]
                        incluster = elem[1]
                    #    print "second incluster ",incluster
                        if not incluster:
                            ci = elem[0]
                            comc = fr.lipidcom[ci].com
                            #dcom = comc-coms
                            dx = comc[xi]-coms[xi]
                            dy = comc[yi]-coms[yi]
                            #rxy = np.sqrt(dx*dx+dy*dy)
                            #print dx," ",dy," ",rxy
                            
                            #Minimum image -- coordinates must be pre-wrapped 
                            if np.absolute(dx) > boxc_xh:
                                dx = boxc[xi] - np.absolute(comc[xi]-boxc_xh) - np.absolute(coms[xi]-boxc_xh)
                            if np.absolute(dy) > boxc_yh:
                                dy = boxc[yi] - np.absolute(comc[yi]-boxc_yh) - np.absolute(coms[yi]-boxc_yh)
                            rxy = np.sqrt(dx*dx+dy*dy)
                            #print "rxy ",rxy," dx ",dx," dy ",dy
                            if    rxy <= dist:
                                #print "passed! adding ",masterlistf[mindex][0]," to the neighborlist"
                                neighborlist.append(masterlistf[j][0])
                                masterlistf[j][1]=True
                        #mindex+=1
                    i+=1
                #filter the masterlistf
            #    print "neighlist", neighborlist
                masterlistf=list([v for v in masterlistf if v[1] == False])
                if len(neighborlist) > 1:
                    clusters.append([])
                    clusters[clustind]=list(neighborlist)
                    #print "clustind clusters[clustind]"
                    #print clustind, " ",clusters
                    clustind+=1
                        
            #print masterlistf
            #filter out single points
            #clusters = [v for v in clusters if len(v) > 1]
            nclusters = len(clusters)
            clsizestat = RunningStats()
            mini = 100000000
            maxi = -1000000
            for cluster in clusters:
                size = len(cluster)
                clsizestat.push(size)
                if    size>maxi:
                    maxi=size
                if size < mini:
                    mini=size
            avgsize = clsizestat.mean()
            #store instantaneous values
            outdata[f,0] = ctime
#            outdata[f,1]= nclusters
#            outdata[f,2] = avgsize
#            outdata[f,3] = mini
#            outdata[f,4] = maxi
            #push to the time averages
            ncstat.push(nclusters)
            asstat.push(avgsize)
            misstat.push(mini)
            masstat.push(maxi)
            #store current time averages
#            outdata[f,5] = ncstat.mean()
#            outdata[f,6] = ncstat.deviation()
#            outdata[f,7] = asstat.mean()
#            outdata[f,8] = asstat.deviation()
#            outdata[f,9] = misstat.mean()
#            outdata[f,10] = misstat.deviation()
#            outdata[f,11] = masstat.mean()
#            outdata[f,12] = masstat.deviation()
            outdata[f,1] = ncstat.mean()
            outdata[f,2] = ncstat.deviation()
            outdata[f,3] = asstat.mean()
            outdata[f,4] = asstat.deviation()
#            outdata[f,5] = misstat.mean()
#            outdata[f,6] = misstat.deviation()
#            outdata[f,7] = masstat.mean()
#            outdata[f,8] = masstat.deviation()
            # now add cluster list to the system storage
            self.clusters.append(list(clusters))
            #print clusters
           # print "Frame ",f
           # print "There are ",nclusters," clusters with an average size of ",avgsize
           # print "the largest cluster was ",maxi," and the smallest was ",mini

        return outdata    
    #takes the cluster lists from self.clusters and gets the plane coordinates 
    # need to call the 'CheckClustering' function before calling this one
    def export_clusters_for_plotting(self):
        """ Determine physical cluster for select lipid COMs in a select leaflet(s).

        This function determines physical clusters based on Cartesian distance criteria
        in a specified leaflet ('upper', 'lower', or 'both') and for a specified LipidGroup within the 
        chosen Leaflet. The routine determines the physical clusters for each frame in the trajectory
        and computes properties of these clusters, e.g. number of clusters, which are time averaged over the trajectory. 
        The information on the clusters is stored in MemSys.clusters which can then be used by other functions 
        (like MemSys.ExportClustersForPlotting).    
            
        Args:           
            leaflet (str, optional): A string designating which Leaflet to include in the computation.
            group (str, optional): A string with the name of a specific LipidGroup to in the computation.
            dist (float, optional): The distance cutoff for determining clusters.
            
        Returns:
            numpy.array: This is a nframesx5 numpy array (of floats) containing the 
                results of this computation across all frames.
                cluster_results[i,0] = simulation time for frame i.
                cluster_results[i,1] = the current running time averaged number of clusters at frame f = i.
                cluster_results[i,2] = the current standard deviation of the running time averaged number of clusters at frame f = i.
                cluster_results[i,3] = the current running time averaged configurational average lipids per cluster at frame f = i.
                cluster_results[i,4] = the current standard deviation of the running time averaged 
                    configurational average lipids per cluster at frame f = i.
                For i in range( nframes ).
        """
        if len(self.clusters) == 0:
            print "Warning!! - call to \'ExportClustersForPlotting\' of a MemSys object with no cluster lists"
            print "      the \'CheckClustering\' function needs to be called first!"
            return 
        xi = self.plane[0]
        yi = self.plane[1]
        #get the maximum number of clusters from any of the frames
        maxsize = 0
        for f in xrange(len(self.clusters)):
            nclust = len(self.clusters[f])
            if nclust>maxsize:
                maxsize=nclust
        #generate a color array
        colors = cm.rainbow(np.linspace(0, 1, maxsize))
        output = []
        for f in xrange(len(self.clusters)):
            frame_clusters = self.clusters[f]
            frame_data = []
            nclust = len(frame_clusters)
            #print len(frame_clusters)
            #print len(colors)
            c = 0
            xcoord = []
            #xm1 = []
            #xp1 = []
            ycoord = []
            #ym1 = []
            #yp1 =[]
            coord_color = []
            for cluster in frame_clusters:
                for index in cluster:
                    xc = self.frame[f].lipidcom[index].com[xi]
                    #xcm1 = self.frame[f].lipidcom[index].com[xi]-self.frame[f].box[xi]
                    #xcp1 = self.frame[f].lipidcom[index].com[xi]+self.frame[f].box[xi]
                    yc = self.frame[f].lipidcom[index].com[yi]
                    #ycm1 = self.frame[f].lipidcom[index].com[yi]-self.frame[f].box[yi]
                    #ycp1 = self.frame[f].lipidcom[index].com[yi]+self.frame[f].box[yi]
                    xcoord.append(xc)
                    #xm1.append(xcm1)
                    #xp1.append(xcp1)
                    ycoord.append(yc)
                    #ym1.append(ycm1)
                    #yp1.append(ycp1)
                    #print c," ",colors[c]
                    coord_color.append(colors[c])
                c+=1    
            #output.append([xm1,xcoord,xp1,ym1,ycoord,yp1,coord_color])
            output.append([xcoord,ycoord,coord_color])
        return output

    # function to compute an approximation of the area per lipid of a group using 
    # closest neighbor circles 
    def calc_area_per_lipid_closest_neighbor_circle(self, leaflet="both",group="all"):
        """ Approximate the area per lipid for select lipid COMs in a select leaflet(s).

        This function computes an approximation of the area per lipid of a selection of lipid COMs
        by locating the closest neighbor and computing the non-overlapping circular area the two 
        equiradius circles formed by the vector between the two COM points. At each frame this is averaged 
        over the configuration.    
            
        Args:           
            leaflet (str, optional): A string designating which Leaflet to include in the computation.
            group (str, optional): A string with the name of a specific LipidGroup to in the computation.
            
        Returns:
            numpy.array: This is a nframesx4 numpy array (of floats) containing the 
                results of this computation across all frames.
                apl_results[i,0] = simulation time for frame i.
                apl_results[i,1] = configurational average area per lipid at frame f = i.
                apl_results[i,2] = the current running time average (of the configurational average) at frame f = i.
                apl_results[i,3] = the current standard deviation of the running time average at frame f = i.
                For i in range( nframes ).
        """
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        do_leaflet = []
        nlip = 0
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
            nlip=self.nlipids
        
        elif leaflet == "upper" or leaflet == "lower":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        sub_fact = (2.0*np.pi/3.0 - np.sqrt(3.0)/2.0)
        #initialize a numpy array to hold the msd for the selection        
        areas = np.zeros((self.nframes, 4))
        #initialize a running stats object to do the averaging
        area_stat = RunningStats()
        n_leaflet = len(do_leaflet)
        #build the index lists
        indices_leaflet = {}
        all_mem_leaflet = {}
        for leaflets in do_leaflet:
            indices = list()
            curr_leaf = self.leaflets[leaflets]
            indices+=curr_leaf.get_group_indices(group)
            n_com = len(indices)
            all_mem = list(self.leaflets[leaflets].get_member_indices())
            all_mem_leaflet[leaflets] = list(all_mem)
            indices_leaflet[leaflets]=list(indices)
            
        
        #loop over the frames
        for f in xrange(self.nframes): 
            fr = self.frame[f]
            dt = fr.time
            boxc=fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            lat_area = boxc_xh*boxc_yh*4.0
            if leaflet == 'both':
                lat_area*=2.0
            
            area_stat_config = RunningStats()
            #loop over the leaflets
            for leaflets in do_leaflet:
                indices = indices_leaflet[leaflets]
                all_mem = all_mem_leaflet[leaflets]
                #loop over the group indices in this leaflet
                for    index in indices:
                    comc = fr.lipidcom[index].com[:]
                    rdist_min = 10000.0
                    #loop over the COMs of non group 
                    #get all the leaflet members
                    
                    for a in all_mem:
                        #print "a ",a
                        if a != index:
                            comn = fr.lipidcom[a].com[:]
                            dx = comc[xi]-comn[xi]
                            dy = comc[yi]-comn[yi]
                        
                            #Minimum image -- coordinates must be pre-wrapped 
                            if np.absolute(dx) > boxc_xh:
                                dx = boxc[xi] - np.absolute(comc[xi]-boxc_xh) - np.absolute(comn[xi]-boxc_xh)
                            if np.absolute(dy) > boxc_yh:
                                dy = boxc[yi] - np.absolute(comc[yi]-boxc_yh) - np.absolute(comn[yi]-boxc_yh)
                            rxy = np.sqrt(dx*dx+dy*dy)
                            #print "rxy ",rxy," dx ",dx," dy ",dy
                            if    rxy < rdist_min:
                                rdist_min = rxy
                    #got the min dist, now compute area
                    #print "rdist_min ",rdist_min
                    area = np.pi*rdist_min*rdist_min - (rdist_min*rdist_min)*sub_fact
                    area_stat_config.push(area)
            area_conf_avg = area_stat_config.mean()
            area_stat.push(area_conf_avg)
            area_time_run = area_stat.mean()
            area_time_run_dev = area_stat.deviation()
            #print "time ",dt
            areas[f][0]=dt
            areas[f][1]=area_conf_avg
            areas[f][2]=area_time_run
            areas[f][3]=area_time_run_dev
            #areas[f][4]=lat_area/nlip
        return areas

    # function to compute the area per lipid using the lateral box sizes and numbers of lipids:
    def calc_area_per_lipid_box(self, leaflet="both"):
        """ Approximate the area per lipid for select lipid COMs in a select leaflet(s).

        This function computes a composite approximation of the area per lipid by simply 
        dividing lateral area of the simulation box by the number lipids in a specific 
        leaflet. 
            
        Args:           
            leaflet (str, optional): A string designating which Leaflet to include in the computation.
                Default is 'both': the area per lipid is averaged over the two leaflets.
            group (str, optional): A string with the name of a specific LipidGroup to in the computation.
            
        Returns:
            numpy.array: This is a nframesx4 numpy array (of floats) containing the 
                results of this computation across all frames.
                apl_results[i,0] = simulation time for frame i.
                apl_results[i,1] = area per lipid at frame f = i.
                apl_results[i,2] = the current running time average (area per lipid) at frame f = i.
                apl_results[i,3] = the current standard deviation of the running time average at frame f = i.
                For i in range( nframes ).
        """
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        do_leaflet = []
        nlip = 0
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
            nlip = []
            for leaflets in do_leaflet:
                nlip.append(float(len(self.leaflets[leaflets])))
        
        elif leaflet == "upper":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        elif leaflet == "lower":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        
        #initialize a numpy array to hold the msd for the selection        
        areas = np.zeros((self.nframes, 4))
        #initialize a running stats object to do the averaging
        area_stat = RunningStats()
        n_leaflet = len(do_leaflet)
            
        
        #loop over the frames
        for f in xrange(self.nframes): 
            fr = self.frame[f]
            dt = fr.time
            boxc=fr.box
            boxc_xh = boxc[xi]/2.0
            boxc_yh = boxc[yi]/2.0
            lat_area = boxc_xh*boxc_yh*4.0
            area_per_lip = lat_area/nlip
            if leaflet == 'both':
                area_per_lip = (lat_area/2.0)*( (nlip[0]+nlip[1])/(nlip[0]*nlip[1]))
            
            area_stat.push(area_per_lip)
            area_time_run = area_stat.mean()
            area_time_run_dev = area_stat.deviation()
            areas[f][0]=dt
            areas[f][1]=area_per_lip
            areas[f][2]=area_time_run
            areas[f][3]=area_time_run_dev
        return areas

    # do Voronoi tesselation using the COMs as generators
    def voronoi_tesselate(self, leaflet="both",group="all"):
        indices = []
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        out_tess = []
        for    f in xrange(self.nframes):
            # get the current frame
            curr_frame = self.frame[f]
            # get the coordinates for the selection at this frame
            com_curr = np.zeros((n_com,2))
            count=0
            for i in indices:
                com_i = curr_frame.lipidcom[i].com_unwrap[self.plane]
                com_curr[count]=com_i[:]
                count+=1
            vor = Voronoi(com_curr)
            #out_tess.append([com_curr[:,0],com_curr[:,1],vor])
            out_tess.append(vor)
        return out_tess

    # do Delauny tesselation using the COMs as generators
    def delaunay_tesselate(self, leaflet="both",group="all"):
        indices = []
        
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        out_tess = []
        for    f in xrange(self.nframes):
            # get the current frame
            curr_frame = self.frame[f]
            # get the coordinates for the selection at this frame
            com_curr = np.zeros((n_com,2))
            count=0
            for i in indices:
                com_i = curr_frame.lipidcom[i].com_unwrap[self.plane]
                com_curr[count]=com_i[:]
                count+=1
            tri = Delaunay(com_curr)
            out_tess.append([com_curr[:,0],com_curr[:,1],tri])
        return out_tess


    # generate the step vectors of the center of mass--in the lateral dimensions
    def step_vector(self, leaflet="both",group="all",fstart=0,fend=-1,fstep=1000,wrapped=False):
        
        indices = []
        if fstart<0:
            fstart+=self.nframes
        if fend < 0:
            fend+=self.nframes
        if fstep == 'single':
            fstep = fend-fstart
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
    
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
            
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        n_com = len(indices)

        #print "there are ",len(indices)," members"
        xi = self.plane[0]
        yi = self.plane[1]
        zi = self.norm
        
        vec_ends_out = []
        for    f in xrange(fstart+fstep,fend+1,fstep):
            fprev = f-fstep
            # get the current frame
            curr_frame = self.frame[f]
            prev_frame = self.frame[fprev]
            # get the coordinates for the selection at this frame
            vec_ends = np.zeros((n_com,4))
            #vec_ends = []
            count=0
            for i in indices:
                com_i = curr_frame.lipidcom[i].com_unwrap[self.plane]
                com_j = prev_frame.lipidcom[i].com_unwrap[self.plane]
                com_j_w = prev_frame.lipidcom[i].com[self.plane]
                if wrapped: 
                    vec_ends[count,0]=com_j_w[0]
                    vec_ends[count,1]=com_j_w[1]
                else:
                    vec_ends[count,0]=com_j[0]
                    vec_ends[count,1]=com_j[1]
                vec_ends[count,2]=com_i[0] - com_j[0]
                vec_ends[count,3]=com_i[1] - com_j[1]
            #    vec_ends.append([com_j[0],com_j[0],com_i[0]-com_j[0],com_i[1]-com_j[1]])
                count+=1
            vec_ends_out.append(vec_ends)
            
        return vec_ends_out

    # return the MemSys frame numbers associated with step vectors calculation
    def step_vector_frames(self,fstart=0,fend=-1,fstep=1000):
        
        
        if fstart<0:
            fstart+=self.nframes
        if fend < 0:
            fend+=self.nframes
        if fstep == 'single':
            fstep = fend-fstart
        output = []
        
        for f in xrange(fstart+fstep,fend+1,fstep):
            fprev = f-fstep
            output.append([fprev, f])
        return np.array(output, dtype=np.int)
    
    # generate the step vectors of the center of mass
    def step_vector_colors(self, leaflet="both",group="all"):
        indices = []            
        ngroups = 1
        group_names = []
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
                curr_group_names = curr_leaf.get_group_names()
                if group == 'all':
                    for gname in curr_group_names:
                        if gname not in group_names:
                            group_names.append(gname)
                else:
                    group_names.append(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
            curr_group_names = curr_leaf.get_group_names()
            if group == 'all':
                for gname in curr_group_names:
                    if gname not in group_names:
                        group_names.append(gname)
            else:
                group_names.append(group)
                
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
            curr_group_names = curr_leaf.get_group_names()
            if group == 'all':
                for gname in curr_group_names:
                    if gname not in group_names:
                        group_names.append(gname)
            else:
                group_names.append(group)
                
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
                curr_group_names = curr_leaf.get_group_names()
                if group == 'all':
                    for gname in curr_group_names:
                        if gname not in group_names:
                            group_names.append(gname)
                else:
                    group_names.append(group)
                
        n_com = len(indices)
        ngroups = len(group_names)
        colors = cm.rainbow(np.linspace(0, 1, ngroups))
        #build color map
        cmap = {}
        n = 0
        for name in group_names:
            cmap[name] = colors[n]
            n+=1
        #pick a frame-just use first frame    
        curr_frame = self.frame[0]
        colors_out = np.zeros( (n_com,4))
        count=0
        for i in indices:
            name_i = curr_frame.lipidcom[i].type
            colors_out[count] = cmap[name_i]
            count+=1
            
        return (colors_out, cmap)

    def remove_leaflet_com_motion(self,leaflet="both"):
        """ Remove the independent center of mass (COM) motion of the Leaflets.  
            
        Args:           
            leaflet (str): A string designating which Leaflet to include.
            
        """
        do_leaflet = []
        nlip = 0
        if leaflet == "both":
            do_leaflet.append('upper')
            do_leaflet.append('lower')
            nlip = []
            for leaflets in do_leaflet:
                nlip.append(float(len(self.leaflets[leaflets])))
        
        elif leaflet == "upper":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        elif leaflet == "lower":
            do_leaflet.append(leaflet)
            nlip = len(self.leaflets[leaflet])
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            do_leaflet.append('upper')
            do_leaflet.append('lower')
        leaf_indices = {}

        for leaf in do_leaflet:
            leaf_indices[leaf]=list(self.leaflets[leaf].get_member_indices())
        
        
        for f in xrange(self.nframes):
            fr = self.frame[f]
            
            for leaf in do_leaflet:
                indices=leaf_indices[leaf]
                #get the leaflet COM
                lcom = np.zeros(3)
                masst = 0.0
                for i in indices:
                    lcom+=(fr.lipidcom[i].com_unwrap*fr.lipidcom[i].mass)
                    masst+=fr.lipidcom[i].mass
                lcom/=masst
                lcom[2]=0.0
                for i in indices:
                    fr.lipidcom[i].com_unwrap-=lcom
            self.frame[f]=fr
        return

    


    ############### multiprocessor parallelized versions of calculation member functions

    # parallelized version of calc_msd- using the multiprocessing module 
    def calc_msd_parallel(self, leaflet="both",group="all",nprocs=2):            
        """ Parallelized version of calc_msd

        This function allows the mean squared displacement (MSD) to be computed
        for a specified leaflet ('upper', 'lower', or 'both') and for a specified LipidGroup within the 
        chosen Leaflet. This calculation is over the whole trajectory and assumes the first frame is the 
        reference frame for computing the displacement. The calculation is parallized over the 
        number of frames in the trejectory using the multiprocessor module. 
            
        Args:           
            leaflet (str): A string designating which Leaflet to include in the computation.
            group (str): A string with the name of a specific LipidGroup to in the computation. 
            
        Returns:
            numpy.array: This is a nframesx2 numpy array (of floats) containing the 
                results of the MSD computation across all frames.
                msd_results[i,0] = simulation time for frame i.
                msd_results[i,1] = the configurational average MSD over the specified LipidGroup for frame f = i.
                For i in range( nframes ).
        """
        indices = []
        #diffusion dimension - assume lateral so, dim=2
        dim=2
        if leaflet == "both":
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        elif leaflet == "upper":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        elif leaflet == "lower":
            curr_leaf = self.leaflets[leaflet]
            indices=curr_leaf.get_group_indices(group)
        else:
            #unknown option--use default "both"
            print "!! Warning - request for unknown leaflet name \'",leaflet,"\' from the ",self.name," leaflet"
            print "!! the options are \"upper\", \"lower\", or \"both\"--using the default \"both\""
            for leaflets in self.leaflets:
                curr_leaf = self.leaflets[leaflets]
                indices+=curr_leaf.get_group_indices(group)
        n_com = len(indices)

        frame_ranges = []
        total_frames = self.nframes
        frames_per_proc_base = total_frames/nprocs
        left_over = total_frames % (frames_per_proc_base * nprocs)
#        print "total frames ",total_frames
#        print "frames per proc ",frames_per_proc_base
#        print "left over ",left_over
        #assign base ranges
        for i in xrange(nprocs):
            fs = i*frames_per_proc_base
            fe = fs + frames_per_proc_base - 1
            frame_ranges.append([fs,fe])
#        print "frame_ranges (pre-adjust):"    
#        print frame_ranges
        #now adjust for leftovers - divide them "equally" over the processes
        lo = left_over
        while lo > 0:
            for i in xrange(nprocs):
                frame_ranges[i][1]+=1
                for j in xrange(i+1,nprocs):
                    frame_ranges[j][0]+=1
                    frame_ranges[j][1]+=1
                lo-=1
                if lo == 0:
                    break
        
#        print "nprocs ",nprocs
#        print "frame_ranges (post adjust): "
#        print frame_ranges
        #initialize a numpy array to hold the msd for the selection        
        msd = np.zeros((self.nframes, 2))
        #
        msd_frames = msd_frames
        #frames_local = getattr(self, 'frame')
        #shelf_local = shelve.open(self.frame.fs_name,flag="r", protocol=2)
        frames_local = par_frames(self.frame.nframes,self.frame.fs_name,self.frame.frame_shelf)
        #frames_local = par_frames(self.frame.nframes,self.frame.fs_name)
        #frames_local = par_frames(self.frame.nframes,self.frame.fs_name,shelf_local)
        plane_local = self.plane
        #create process pool
        pool = mp.Pool(processes=nprocs)
        results = [pool.apply_async(msd_frames,args=(frames_local,frame_ranges[i][0],frame_ranges[i][1],indices,0,plane_local)) for i in range(0,nprocs)]
    #    print "results:"
    #    print results
        results_ordered = [p.get() for p in results]
    #    print "results ordered: "
    #    print results_ordered
#        #collect results  into single array for return
        i = 0
    #    print "len(results_ordered) ",len(results_ordered)
        for p in results_ordered:
            fs = frame_ranges[i][0]
            fe = frame_ranges[i][1]
            #print fs, fe
            #print msd[fs:(fe+1)].shape
            #print p[:].shape
            msd[fs:(fe+1)] = p[:]
            i+=1
        pool.close()
        pool.join()        

        return msd
    #function to compute the thickness of the membrane (in the normal direction). The algorithm is based on  
    # the GridMAT-MD bilayer thickness calculation (except without the gridding procedure) 
    def calc_membrane_thickness_parallel(self,nprocs=2):
        """ Parallelized version of MemSys.CalcMembraneThickness 

        Computes the thickness of the bilayer (along the normal direction) across the trajectory.
        The algorithm is based on the GridMAT-MD bilayer thickness calculation. The calculation is 
        parallized over the number of frames in the trejectory using the multiprocessor module.  
        This function passes the thickness_frames fucntion to the multiprocessor threads.   

        Args:   
            nprocs (int): The number processors (threads) to use in the computation.
            
        Returns:
            tuple of numpy.array: This is a two element tuple containing numpy arrays of the computation results.
                tuple[0] => thickness: A nx5 numpy array containing the 
                results of the thickness computation for the specified frames. Specifically:
                    thickness[i,0] = simulation time for frame f = i + fstart.
                    thickness[i,1] = the configurational average thickness for frame f = i + fstart.
                    thickness[i,2] = the standard deviation of the configurational average thickness for frame f = i + fstart.
                    thickness[i,3] = the running time average of the configurational average thickness for frame f = i + fstart.
                    thickness[i,4] = the running standard deviation of the time averaged configurational average thickness for frame f = i + fstart.
                    For i in xrange( (fend-fstar) + 1).
                tuple[1] => thickness_map: A nxNx6 numpy array containing the thickness data that can be 
                used to generate a 3d thickness map/plot. Specifically:
                    thickness[i,j,0] = simulation time for frame f = i + fstart and lipid j.
                    thickness[i,j,1] = the average x position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,2] = the average y position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,3] = the lower z position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,4] = the upper z position for lipid j and its cross leaflet partner at frame f = i + fstart.
                    thickness[i,j,5] = the difference between the upper and lower z positions 
                    for lipid j and its cross leaflet partner at frame f = i + fstart.
                    For i in xrange((fend-fstar) + 1) and For j in xrange(nlipids).
        """
        nlip = self.nlipids
        comcup = np.zeros(3)
        comclo = np.zeros(3)
        dcom = np.zeros(3)
        zdists = np.zeros((self.nframes, 3))
        zmaps = np.zeros((self.nframes, self.nlipids, 6))
        frame_ranges = []
        total_frames = self.nframes
        frames_per_proc_base = total_frames/nprocs
        left_over = total_frames % (frames_per_proc_base * nprocs)
#        print "total frames ",total_frames
#        print "frames per proc ",frames_per_proc_base
#        print "left over ",left_over
        #assign base ranges
        for i in xrange(nprocs):
            fs = i*frames_per_proc_base
            fe = fs + frames_per_proc_base - 1
            frame_ranges.append([fs,fe])
#        print "frame_ranges (pre-adjust):"    
#        print frame_ranges
        #now adjust for leftovers - divide them "equally" over the processes
        lo = left_over
        while lo > 0:
            for i in xrange(nprocs):
                frame_ranges[i][1]+=1
                for j in xrange(i+1,nprocs):
                    frame_ranges[j][0]+=1
                    frame_ranges[j][1]+=1
                lo-=1
                if lo == 0:
                    break
        
#        print "nprocs ",nprocs
#        print "frame_ranges (post adjust): "
#        print frame_ranges
#        
        thick_frames = thickness_frames
        frames_local = par_frames(self.frame.nframes,self.frame.fs_name,self.frame.frame_shelf)
        plane_local = self.plane
        norm_local = self.norm
        #create process pool
        pool = mp.Pool(processes=nprocs)
        results = [pool.apply_async(thick_frames,args=(frames_local,frame_ranges[i][0],frame_ranges[i][1],self.leaflets,nlip,plane_local,norm_local)) for i in range(0,nprocs)]
        #print "results:"
    #    print results
        #print "len(results) ",len(results)
        results_ordered = [p.get() for p in results]
        #print "results ordered: "
    #    print results_ordered
#        #collect results  into single array for return
        i = 0
        #print "len(results_ordered) ",len(results_ordered)
        for p in results_ordered:
            fs = frame_ranges[i][0]
            fe = frame_ranges[i][1]
            print fs, fe
            #print msd[fs:(fe+1)].shape
            #print p[:].shape
            zdistf = p[0]
            zmapf = p[1]
            #print zdistf.shape," ",zmapf.shape
            zdists[fs:(fe+1)] = zdistf[:]
            zmaps[fs:(fe+1)] = zmapf[:]
            #zdists[fs:(fe+1)] = pg[:]
            i+=1
        pool.close()
        pool.join()
                
        #regenerate the container
        zdist_tavg = np.zeros((self.nframes, 5))
        # get the running time average
        tavg_dz = gen_running_average(zdists[:,1])
        #slice together the values
        zdist_tavg[:,0:3]=zdists[:,:]
        zdist_tavg[:,3:5]=tavg_dz[:,:]
            
            
        #shelf_local.close()
        return (zdist_tavg, zmaps)
        #return zdist_tavg


