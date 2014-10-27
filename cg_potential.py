#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.22"
parser = argparse.ArgumentParser(prog = 'cg_potential', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/cg_potential
DOI: 
**********************************************

[ DESCRIPTION ]
 
This script calculates the electrosatic potential along z in CG systems - i.e. taking
into account the force shift and cutoff.

A file containing the charged particles can also be supplied to calculate the density
of charges.
 
The script follows a 2-step process:
 1. calculate charge density along z for each slice
 2. sum the corresponding shifted potentials over all the slices

density calculation
-------------------
  You can specify which particles to take into account for the calculation of the total
  charge density by supplying a file via the --charges option. Each line of this file
  should follow the format (without quotation marks):
   -> 'label,value,MDAnalysis selection string'

  The absolute charge for each group will be plotted on the charge density profile. The
  group colour must be specified for each charge.
  
  By default the charged are defined as follows:
   -> Na+,1,name Na+
   -> Cl-,-1,name CL-
   -> PO4,-1,name PO4
   -> NC3,1,name NC3
   -> NH3,1,name NH3
  (note that the MDAnalysis selection string should not contain any commas)


electrostatic potential
-----------------------
  The shifted electrostatic potential is calculated as per Gromacs manual for each point charge.
  In particular when rs = 0, the potential created by a charge q at a distance r is:
   -> phi(r) = q/(4*pi*e0*er)*(1/r - 5/(3*rc) + 5*r^3/(3*rc^4) - r^4/rc^5)


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - numpy
 - scipy
 - networkX (if option --algorithm is set to 'min' or 'cog')
 - sklearn (if option --algorithm is set to 'density')


[ NOTES ]

1. The density is calculated with respect to the z axis, not the bilayer normal. So the
   more your system deforms the noiser the less meaningful the results get.

 
[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		10	: process every t-frames
 
Density profile options
-----------------------------------------------------
--charges		: definition of charged particles, see 'DESCRIPTION'  [TO DO]
--slices	[200] 	: number of slizes along z

Electrostatic potential
-----------------------------------------------------
--er		[15]	: dielectric constant
--rs		[0] 	: distance from charge where the electrostatic starts to be shifted (Angstrom)
--rc		[12]	: distance from charge where the electrostatic should reach 0 (Angstrom)

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[10], type=int, help=argparse.SUPPRESS)

#density profile options
parser.add_argument('--charges', nargs=1, dest='chargesfilename', default=['mine'], help=argparse.SUPPRESS)
parser.add_argument('--slices', nargs=1, dest='slices', default=[200], type=int, help=argparse.SUPPRESS)

#electrostatic potential
parser.add_argument('--er', nargs=1, dest='er', default=[15], type=float, help=argparse.SUPPRESS)
parser.add_argument('--rs', nargs=1, dest='rs', default=[0], type=float, help=argparse.SUPPRESS)
parser.add_argument('--rc', nargs=1, dest='rc', default=[12], type=float, help=argparse.SUPPRESS)

#lipids identification options
parser.add_argument('--beads', nargs=1, dest='beadsfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]
#density profile options
args.chargesfilename = args.chargesfilename[0]
args.slices = args.slices[0]
##electrostatic potential
args.er = args.er[0]
args.rs = args.rs[0]/float(10)				#conversion to nm for unit consistency
args.rc = args.rc[0]/float(10)				#conversion to nm for unit consistency

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.chargesfilename != "no" and args.chargesfilename != "mine" and not os.path.isfile(args.chargesfilename):
	print "Error: file " + str(args.chargesfilename) + " not found."
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)

if args.xtcfilename == "no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "cg_potential_" + args.grofilename[:-4]
	else:
		args.output_folder = "cg_potential_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	filename_log = os.getcwd() + '/' + str(args.output_folder) + '/cg_potential.log'
	output_log = open(filename_log, 'w')		
	output_log.write("[cg_potential v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python cg_potential.py"
	for c in sys.argv[1:]:
		tmp_log += " " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	
	#copy input files
	if args.chargesfilename != "no" and args.chargesfilename != "mine":
		shutil.copy2(args.chargesfilename,args.output_folder + "/")

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

global f_factor
global bins
global vols
global thicks
global potential
global charge_density
global bins_stats
global potential_stats
global charge_density_stats

f_factor = 138.935485/float(args.er)
bins_stats = {}
bins_stats["avg"] = np.zeros(args.slices+1)
bins_stats["std"] = np.zeros(args.slices+1)
potential_stats = {}
potential_stats["avg"] = np.zeros(args.slices)
potential_stats["std"] = np.zeros(args.slices)
charge_density_stats = {}
charge_density_stats["avg"] = np.zeros(args.slices)
charge_density_stats["std"] = np.zeros(args.slices)

#=========================================================================================
# data loading
#=========================================================================================

def set_charges():
			
	global charges_groups
	charges_groups = {}
	
	charges_groups["Na+"] = {}
	charges_groups["Na+"]["value"] = 1
	charges_groups["Na+"]["sele_string"] = "name NA+"
	
	charges_groups["Cl-"] = {}
	charges_groups["Cl-"]["value"] = 1
	charges_groups["Cl-"]["sele_string"] = "name Cl-"

	charges_groups["PO4"] = {}
	charges_groups["PO4"]["value"] = -1
	charges_groups["PO4"]["sele_string"] = "name PO4"

	charges_groups["NH3"] = {}
	charges_groups["NH3"]["value"] = 1
	charges_groups["NH3"]["sele_string"] = "name NH3"

	charges_groups["NC3"] = {}
	charges_groups["NC3"]["value"] = 1
	charges_groups["NC3"]["sele_string"] = "name NC3"
		
	return
def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global nb_frames_to_process
	global f_start
	global f_end
	global residues_list
	global water_pres
	global water_sele
	global bilayer
	f_start = 0
		
	#load universe
	#-------------
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		U_timestep = U.trajectory.dt
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes

		U.trajectory.rewind()
		#sanity check
		if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#create list of index of frames to process
		if args.t_end != -1:
			f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_end < 0:
				print "Error: the starting time specified is before the beginning of the xtc."
				sys.exit(1)
		else:
			f_end = nb_frames_xtc - 1
		if args.t_start != -1:
			f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_start > f_end:
				print "Error: the starting time specified is after the end of the xtc."
				sys.exit(1)
		if (f_end - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)

	#check 'bilayer' selection (only useful to center the z coordinates)
	#-------------------------------------------------------------------
	bilayer = U.selectAtoms("name PO4")

	#create charged particles selections
	#-----------------------------------
	charge_pres = False
	for q in charges_groups.keys():
		charges_groups[q]["sele"] = U.selectAtoms(charges_groups[q]["sele_string"])
		if charges_groups[q]["sele"].numberOfAtoms() == 0:
			print " ->warning: charge selection string '" + str(charges_groups[q]["sele_string"]) + "' returned 0 atoms."
		else:
			charge_pres = True
	if not charge_pres:
		print "Error: no charged particles found, try using --charges to supply correct charges definition."
		sys.exit(1)
		
	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_potential(f_index, box_dim):
	
	#define bins
	tmp_bins = np.linspace(0,box_dim[2],args.slices+1)
	
	#store their z values and volume
	thicks[f_index] = tmp_bins[1]-tmp_bins[0]
	vols[f_index] = thicks[f_index] * box_dim[0]*box_dim[1]
	bins[f_index,:] = tmp_bins - bilayer.centerOfGeometry()[2]
	
	#calculate charge density in each bin
	for q in charges_groups.keys():
		tmp_q_sele = charges_groups[q]["sele"]
		if tmp_q_sele.numberOfAtoms() > 0:
			tmp_coord = tmp_q_sele.coordinates()
			charge_density[f_index,:] += np.histogram(tmp_coord[:,2], tmp_bins)[0] * charges_groups[q]["value"]
	
	#calculate potential
	
	return
def calculate_stats():
	
	for s in range(0, args.slices):
		#bin z position
		bins_stats["avg"][s] = np.average(bins[:,s])
		bins_stats["std"][s] = np.std(bins[:,s])

		#charge density
		charge_density_stats["avg"][s] = np.average(charge_density[:,s])
		charge_density_stats["std"][s] = np.std(charge_density[:,s])	
		
		#potential 
		potential_stats["avg"][s] = np.average(potential[:,s])
		potential_stats["std"][s] = np.std(potential[:,s])
	
	#calculate bin labels
	#--------------------
	bins_stats["avg"][args.slices] = np.average(bins[:,args.slices])
	bins_stats["std"][args.slices] = np.std(bins[:,args.slices])
	for s in range(0, args.slices):
		bins_labels[s] = (bins_stats["avg"][s+1]-bins_stats["avg"][s])/float(2)
	
	return

#=========================================================================================
# outputs
#=========================================================================================

def density_write_charges():
	
	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/density_profile_charges.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [charge density profile - written by cg_potential v" + str(version_nb) + "]\n")
	output_xvg.write("#  -> nb of slices: " + str(args.slices) + " (Angstrom)\n")
	output_xvg.write("#  -> slices thickness: " + str(round(np.average(thicks),2)) + " (" + str(round(np.std(thicks),2)) + ") (Angstrom)\n")
	output_xvg.write("#  -> slices volume: " + str(round(np.average(vols),2)) + " (" + str(round(np.std(vols),2)) + ") (Angstrom3)\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Charge density profile along z\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"charge density (e.Angstrom-3)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 1\n")
	output_xvg.write("@ s0 legend \"total\"\n")
	
	#data
	for s in range(0,args.slices):
		results = str(bins_labels[s])
		results += str(round(bins_labels[s],2))+ "	" + "{:.6e}".format(charge_density_stats["avg"][s])
		output_xvg.write(results + "\n")	
	output_xvg.close()
		
	return
def density_graph_charges():
			
	#filenames
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/density_profile_charges.svg'
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/density_profile_charges.png'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Charge density profile along z")

	#plot data
	ax = fig.add_subplot(111)
	plt.plot(range(0,args.slices), charge_density_stats["avg"], color = 'k', linewidth = 2)
	#plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, min_density_charges, max_density_charges, linestyles = 'dashed')
	#plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, min_density_charges, max_density_charges, linestyles = 'dashed')
	#plt.vlines(bins_nb, min_density_charges, max_density_charges, linestyles = 'dashdot')
	#plt.hlines(0, 0, 2*bins_nb)
	#fontP.set_size("small")
	#ax.legend(prop=fontP)
	plt.xlabel('z distance to bilayer center [$\AA$]')
	plt.ylabel('average charge density [$e.\AA^{-3}$]')
	
	#save figure
	ax.set_xlim(0, args.slices)
	#ax.set_ylim(min_density_charges, max_density_charges)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	xlabel = ax.get_xticks().tolist()
	#for tick_index in range(0,len(xlabel)):
	#	xlabel[tick_index] -= bins_nb
	ax.set_xticklabels(xlabel)
	ax.xaxis.labelpad = 20
	ax.yaxis.labelpad = 20
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================
set_charges()
load_MDA_universe()

#update global variables
thicks = np.zeros(nb_frames_to_process)
vols = np.zeros(nb_frames_to_process)
bins = np.zeros((nb_frames_to_process, args.slices+1))
bins_labels = np.zeros(args.slices)
potential = np.zeros((nb_frames_to_process, args.slices))
charge_density = np.zeros((nb_frames_to_process, args.slices))

#=========================================================================================
# generate data
#=========================================================================================
print "\nCalculating electrosatic potential..."
#case: structure only
#--------------------
if args.xtcfilename=="no":
	calculate_potential(0,U.trajectory.ts.dimensions)
#case: browse xtc frames
#-----------------------
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
		sys.stdout.flush()
		sys.stdout.write(progress)							
		calculate_potential(f_index, U.trajectory.ts.dimensions)
	print ""

#=========================================================================================
# process data
#=========================================================================================
calculate_stats()

#=========================================================================================
# produce outputs
#=========================================================================================
print "\nWriting outputs..."
density_write_charges()
density_graph_charges()
	
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
