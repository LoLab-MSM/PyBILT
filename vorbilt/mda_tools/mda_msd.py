#we are going to use the MDAnalysis to read in topo and traj
#numpy
import numpy as np
#import my running stats class
from RunningStats import *
# import the coordinate wrapping function--for unwrapping
from pUnwrap import mda_wrap_coordinates
'''
 function to compute the mean square displacement (MSD) and diffusion constant
 of a list of MDAnalysis atom selections (atom_sel_list). The list of atom selections 
 are averaged at each timestep.
 Returns 2d numpy array with len(atom_sel_list)X6 elements: 
 [:,0]=dt [:,1]=msd [:,2]=msd_dev [:,3]=diff_con_instantaneous 
 [:,4]=diff_con_running_average [:,5]=diff_con_running_dev
 
 Long time mean squared displacement: 
	MSD = lim_(t->inf) <||r_i(t) - r_i(0)||**2>_(nsels) = 2*dim*D*t
'''
def mda_msd (trajectory, atom_sel_list, lateral=False, plane="xy", unwrap=True, verbose=False):
	dim=3
	plane_index = [0,1,2]
	if	lateral:
		dim=2		
		ii=0
		jj=1
		if	plane=="yz" or plane=="zy":
			ii=1
			jj=2
		if	plane=="xz" or plane=="zx":
			ii=0
			jj=2
		plane_index = [ii, jj]
	naxes = len(plane_index)
	#get the number of frames from the trajectory
	nframes = len(trajectory)
	#get the number of atomselections
	nsels = len(atom_sel_list)
	#initialize a numpy array to hold the center of mass vectors
	comlist = np.zeros((nsels, nframes, 3))
	#print "l comlist ",len(comlist)
	times = np.zeros(nframes)
	#index counter for the frame number
	comit = 0

	#combine all the selections into one (for wrapping)
	msel = atom_sel_list[0]
	for s in xrange(1, nsels):
		 msel+=atom_sel_list[s]
	
	natoms = len(msel)
	oldcoord = np.zeros((natoms,naxes))
	index = msel.indices
	firstframe = True
	# loop over the trajectory
	for ts in trajectory:
		time=ts.time
		if	verbose:
			print " "
			print "frame ",ts.frame 
	#unwrap coordinates -- currently unwraps all the coordinates
		if	unwrap:	
			if	verbose:
				print "unwrapping frame ",ts.frame 		
			currcoord = ts._pos[index]
			if firstframe:
				oldcoord = currcoord
				firstframe = False
			else:
				abc = ts.dimensions[0:3]
				wrapcoord = mda_wrap_coordinates(abc, currcoord, oldcoord)
				ts._pos[index] = wrapcoord[:]

		#loop over the selections
		for i in xrange(nsels):
			if	verbose:
				print "frame ",ts.frame," getting com of selection ",atom_sel_list[i] 
			#compute the center of mass of the current selection and current frame
			com = atom_sel_list[i].center_of_mass()
			#print "com ",com
			#add to the numpy array
			comlist[i,comit]=com
			#print comlist
		times[comit]=time
		comit+=1
	#initialize a numpy array to hold the msd for each selection		
	msd = np.zeros((nframes, 6))
	#initialize a running stats object to do the averaging
	drs_stat = RunningStats()
	#initialize a running stats object for the diffusion constant (frame/time average)
	diff_stat = RunningStats()
	#loop over the frames starting at index 1
	#print comlist
	#print len(comlist)
	coml0 = comlist[:,0,plane_index]
	#print coml0
	for i in xrange(1, nframes):
		# get the current com frame list
		comlcurr = comlist[:,i,plane_index]
		dr = comlcurr - coml0
		drs = dr*dr
		#loop over the selections for this frame
		for	j in xrange(nsels):
			drs_curr = drs[j,:]	
			drs_mag = drs_curr.sum()
			drs_stat.Push(drs_mag)
		#get the msd for the current selection
		msdcurr = drs_stat.Mean()
		devcurr = drs_stat.Deviation()
		dt = times[i]-times[0]
		DiffCon = msdcurr/(2.0*dim*dt)
		diff_stat.Push(DiffCon)
		#print "msdcurr ",msdcurr
		#push to the msd array
		msd[i,0]=dt
		msd[i,1]=msdcurr
		msd[i,2]=devcurr
		msd[i,3]=DiffCon
		msd[i,4]=diff_stat.Mean()
		msd[i,5]=diff_stat.Deviation()
		if	verbose:
				print "selection number ",i," has msd ",msdcurr," with deviation ",devcurr 
		#reset the running stats object--prepare for next selection
		drs_stat.Reset()
	#return msd array
	return msd 

