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
parser = argparse.ArgumentParser(prog = 'cluster_density_profile', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
******************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/cluster_density_profile
DOI: 
******************************************************

[ DESCRIPTION ]
 
This script plots the density profile along the z axis of groups of particles present
around transmembrane proteins. Density profiles are broken down by the size of the
TM clusters.

A file containing the charged particles can also be supplied to calculate the density
of charges.
 
If you want to calculate densities in a membrane which do not contain any proteins (or
ignore the presence of proteins) just set '--proteins no'.

density calculation
-------------------
Density is calculated based on a cylinder centered on the protein cluster center of 
geometry. All the particles present in this cylinder are binned into cylinder slices
based on their distance (along z) to the center of the lipid bilayer. The cylinder
radius is controlled via --slices_radius and the slices thickness via --slices_thick.

(a) particles selection
  You can specify the particles for which to plot the density by supplying a file
  via the --particles option. Each line of this file should follow the following
  format (without quotation marks):
   -> 'group,label,colour,MDAnalysis selection string'
 
  where 'group' is used to normalise the densities of several particles. By default
  each particle type is normalised with respect to itself and the following densities
  are used:
   -> peptide,peptide,#262626,protein
   -> CHOL,CHOL,#bdbdbd,resname CHOL and name ROH
   -> POPC,POPC,#41ab5d,resname POPC and name PO4
   -> POPE,POPE,#6a51a3,resname POPE and name PO4
   -> POPS,POPS,#cc4c02,resname POPS and name PO4
   -> water,water,#1d91c0,resname W
   -> Na+,Na+,#7bccc4,name NA+
   -> Cl-,Cl-,#fa9fb5,name CL-
  
  Note that the only acceptable protein selection is the whole "protein" string and
  that it should be labeled as "peptide" and in a group of its own also labeled
  "peptide". If you're interested in particular sub-selections of proteins, see (b)
  below.
   
(b) protein details: residue types
  Density profiles can be broken down further for peptides to show the density
  distribution amongst different residue types. You can group the residues as you
  wish by supplying a file via the --residues options. Each line of this file
  should follow the following format:
  -> 'label,colour,resname1,resname2,...'
  
  By default the following residue types are used:
   -> basic,#253494,ARG,LYS
   -> polar,#a1d99b,SER,THR,ASN,GLN,HIS		
   -> negative,#99d8c9,ASP, GLU
   -> hydrophobic,#993404,VAL,ILE,LEU,MET,PRO,CYS,PHE,TYR,TRP
   -> backbone,#969696,ALA,GLY

  WARNING: you need to have defined a particle type called "peptide" in (a) in order
  to show residue details.
  If you do not want to show residue details, just use: '--residues no'.

(c) charge density
  You can specify which particles to take into account for the calculation of the total
  charge density by supplying a file via the --charges option. Each line of this file
  should follow the format (without quotation marks):
   -> 'group_name,colour,charge_name,charge,MDAnalysis selection string for charge'

  The absolute charge for each group will be plotted on the charge density profile. The
  group colour must be specified for each charge.
  
  By default the charged are defined as follows (the peptide charges correspond to that
  of uncapped transportan):
   -> ions,#52A3CC,Na+,1,resname NA+
   -> ions,#52A3CC,CL-,-1,resname CL-
   -> lipids,#b2182b,-1,phosphate,name PO4
   -> lipids,#b2182b,1,amine_choline,name NH3 or name NC3
   -> peptide,#053061,pos,1,(resnum 1 and resname GLY and name BB) or (resname LYS and name SC2)
   -> peptide,#053061,neg,-1, resnum 27 and resname LEU and name BB
  (note that the MDAnalysis selection string should not contain any commas)

  If you do not want to calculate charge density, just use: '--charges no'

(d) colour definition
  Colours can be specified using single letter code (rgbcmykw), hex code  or the name of
  a colour map (see the matplotlib website for a list of the available colour maps).
  In case a colour map is used, its name must be specified as the only colour.

detection of transmembrane protein clusters
-------------------------------------------
Two clustering algorithms can be used to identify protein clusters.

(a) Connectivity based (relies on networkX module):
  A protein is considered in a cluster if it is within a distance less than --nx_cutoff
  from another protein. This means that a single protein can act as a connector between
  two otherwise disconnected protein clusters.
  This algorithm can be ran using either the minimum distante between proteins (default, 
  --algorithm 'min') or the distance between their center of geometry (--algorithm 'cog').
  The 'min' option scales as the square of the number of proteins and can thus be very
  slow for large systems.

(b) Density based (relies on the sklearn module and its implementation of DBSCAN):
  A protein is considered in a cluster if is surrounded by at least --db_neighbours other
  proteins within a radius of --db_radius.
  This density based approach is usually less suited to the detection of protein
  clusters but as a general rule the more compact the clusters, the smaller --db_radius
  the higher --db_neighbours can be - for details on this algorithm see its online
  documentation.
  This algorithm is selected by setting the --algorithm option to 'density'.

The identified protein clusters are considered to be transmembrane only if the closest
lipid headgroup neighbours to the cluster particles are all within the same leaflet.
In addition to the sizes identified, size groups can be defined - see note 7.


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - np
 - sp
 - networkX (if option --algorithm is set to 'min' or 'cog')
 - sklearn (if option --algorithm is set to 'density')


[ NOTES ]

1. The density is calculated with respect to the z axis, not the bilayer normal. So the
   more your system deforms the noiser the less meaningful the results get.

2. Identification of the bilayer leaflets can be controlled via 3 options:
   (a) beads
    By default, the particles taken into account to define leaflet are:e
    -> name PO4 or name PO3 or name B1A
   
    Note that only lipids which contain one of the beads mentioned in the selection string
    will be taken into account. If you wish to specify your own selection string (e.g. to
    choose different beads or add a bead not in the default list in order to take into
    account a particular lipid specie) you can do so by supplying a file via the --beads
    option. This file should contain a single line that can be passed as the argument
    to MDAnalysis selectAtoms() routine and should not contain any quotation marks, e.g.:
     -> name PO4 or name PO3 or name B1A or name AM1
        
   (b) leaflet finding method
    By default leaflets are identified using the MDAnalysis LeafletFinder routine and the
    the optimum cutoff to identify 2 lipids groups is determined using the optimize_cutoff
    routine.
    This optimisation process can take time in large systems and you can specify your own
    cutoff value to skip this step. For instance to use a 15 Angstrom cutoff value:
     -> '--leaflet 15'
   
    In very large systems (more then ~50,000 phospholipids) LeafletFinder (or rather the
    networkX module that it relies on) can fail. To  avoid this you can choose not to use
    this routine by specifying:
     -> '--leaflet large'
    In this case lipids whose headgroups z value is above the average lipids z value will
    be considered to make up the upper leaflet and those whose headgroups z value is below
    the average will be considered to be in the lower leaflet.
    This means that the bilayer should be as flat as possible in the gro file supplied in
    order to get a meaningful outcome.

   (c) flipflopping lipids
    In case lipids flipflop during the trajectory, a file listing them can be supplied
    with the --flipflops option. Each line of this file should follow the format:
     -> 'resname,resid,starting_leaflet,z_bead'
    where starting_leaflet is either 'upper' or 'lower' - e.g. 'POPC,145,lower,PO4'. The
    z_bead is used to track the position of the lipid.
    If flipflopping lipids are not specified they may add significant noise to results as
    they prevent the correct identification of TM clusters.
    Note that, even when specified, flipflopping lipids will be taken into account when
    calculating densities and charges.   

3. Proteins are detected automatically but you can specify an input file to define your
   own selection with the --proteins option.
   In this case the supplied file should contain on each line a protein selection string
   that can be passed as the argument of the MDAnalysis selectAtoms() routine - for 
   instance 'bynum 1:344'.
   
   If you do not have proteins in your bilayer or do not want to take them into account,
   just use '--proteins no'.

4. The densities are calculated for each TM cluster size identified but can also be
   binned into size groups.
   The size groups are defined by supplying a file with --groups, whose lines all
   follow the format:
    -> 'lower_size,upper_size'

   Size groups definition should follow the following rules:
    -to specify an open ended group use 'max', e.g. '3,max'
    -groups should be ordered by increasing size and their boundaries should not overlap
    -boundaries are inclusive so you can specify one size groups with 'size,size,colour'
    -any cluster size not covered will be labeled as 'other'
 
 
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
--range		40 	: distance spanned on either side of the bilayer center
--types			: definition of residue groups, see 'DESCRIPTION'
--charges		: definition of charged particles, see 'DESCRIPTION' 
--slices_thick	0.5 	: z thickness of the slices (Angstrom)
--slices_radius	30 	: radius of the slices (Angstrom)

Lipids identification  
-----------------------------------------------------
--beads			: leaflet identification technique, see note 2(a)
--flipflops		: input file with flipflopping lipids, see note 2(c)
--leaflets	optimise: leaflet identification technique, see note 2(b)

Protein clusters identification
-----------------------------------------------------
--groups		: cluster groups definition file, see note 4
--proteins		: protein selection file, (optional, see note 3)
--algorithm	min	: 'cog','min' or 'density', see 'DESCRIPTION'
--nx_cutoff 	8	: networkX cutoff distance for protein-protein contact (Angstrom)
--db_radius 	20	: DBSCAN search radius (Angstrom)
--db_neighbours	3	: DBSCAN minimum number of neighbours within a circle of radius --db_radius	

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
parser.add_argument('--range', nargs=1, dest='max_z_dist', default=[40], type=float, help=argparse.SUPPRESS)
parser.add_argument('--particles', nargs=1, dest='particlesfilename', default=['mine'], help=argparse.SUPPRESS)
parser.add_argument('--residues', nargs=1, dest='residuesfilename', default=['mine'], help=argparse.SUPPRESS)
parser.add_argument('--charges', nargs=1, dest='chargesfilename', default=['mine'], help=argparse.SUPPRESS)
parser.add_argument('--slices_thick', nargs=1, dest='slices_thick', default=[0.5], type=float, help=argparse.SUPPRESS)
parser.add_argument('--slices_radius', nargs=1, dest='slices_radius', default=[30], type=float, help=argparse.SUPPRESS)

#lipids identification options
parser.add_argument('--beads', nargs=1, dest='beadsfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--flipflops', nargs=1, dest='selection_file_ff', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--leaflets', nargs=1, dest='cutoff_leaflet', default=['optimise'], help=argparse.SUPPRESS)

#radial and protein clusters options
parser.add_argument('--groups', nargs=1, dest='cluster_groups_file', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('--proteins', nargs=1, dest='selection_file_prot', default=['auto'], help=argparse.SUPPRESS)
parser.add_argument('--algorithm', dest='m_algorithm', choices=['cog','min','density'], default='min', help=argparse.SUPPRESS)
parser.add_argument('--nx_cutoff', nargs=1, dest='cutoff_connect', default=[8], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_radius', nargs=1, dest='dbscan_dist', default=[20], type=float, help=argparse.SUPPRESS)
parser.add_argument('--db_neighbours', nargs=1, dest='dbscan_nb', default=[3], type=int, help=argparse.SUPPRESS)

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
args.max_z_dist = args.max_z_dist[0]
args.particlesfilename = args.particlesfilename[0]
args.residuesfilename = args.residuesfilename[0]
args.chargesfilename = args.chargesfilename[0]
args.slices_thick = args.slices_thick[0]
args.slices_radius = args.slices_radius[0]
#lipids identification options
args.beadsfilename = args.beadsfilename[0]
args.cutoff_leaflet = args.cutoff_leaflet[0]
args.selection_file_ff = args.selection_file_ff[0]
#radial and protein clusters options
args.cluster_groups_file = args.cluster_groups_file[0]
args.selection_file_prot = args.selection_file_prot[0]
args.cutoff_connect = args.cutoff_connect[0]
args.dbscan_dist = args.dbscan_dist[0]
args.dbscan_nb = args.dbscan_nb[0]

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
if args.selection_file_ff != "no" and not os.path.isfile(args.selection_file_ff):
	print "Error: file " + str(args.selection_file_ff) + " not found."
	sys.exit(1)
if args.selection_file_prot != "auto" and args.selection_file_prot != "no" and not os.path.isfile(args.selection_file_prot):
	print "Error: file " + str(args.selection_file_prot) + " not found."
	sys.exit(1)
if args.cluster_groups_file != "no" and not os.path.isfile(args.cluster_groups_file):
	print "Error: file " + str(args.cluster_groups_file) + " not found."
	sys.exit(1)
if args.beadsfilename != "no" and not os.path.isfile(args.beadsfilename):
	print "Error: file " + str(args.beadsfilename) + " not found."
	sys.exit(1)
if args.particlesfilename != "mine" and not os.path.isfile(args.particlesfilename):
	print "Error: file " + str(args.particlesfilename) + " not found."
	sys.exit(1)
if args.residuesfilename != "no" and args.residuesfilename != "mine" and not os.path.isfile(args.residuesfilename):
	print "Error: file " + str(args.residuesfilename) + " not found."
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

if args.m_algorithm != "density":
	if '--db_radius' in sys.argv:
		print "Error: --db_radius option specified but --algorithm option set to '" + str(args.m_algorithm) + "'."
		sys.exit(1)
	elif '--db_neighbours' in sys.argv:
		print "Error: --db_neighbours option specified but --algorithm option set to '" + str(args.m_algorithm) + "'."
		sys.exit(1)
else:
	if '--nx_cutoff' in sys.argv:
		print "Error: --nx_cutoff option specified but --algorithm option set to 'density'."
		sys.exit(1)

if args.selection_file_prot == "no" and args.cluster_groups_file != "no":
	print "Error: --groups option specified but --proteins option set to 'no'."
	sys.exit(1)
if args.selection_file_prot == "no" and args.residuesfilename != "no":
	print "Error: --proteins option set to 'no' but option --residues isn't."
	sys.exit(1)

#=========================================================================================
# process options
#=========================================================================================

global lipids_ff_nb
global bins_nb
global bins_nb_max
global bins_labels
global slice_volume
global protein_pres
lipids_ff_nb = 0
bins_nb = int(np.floor(args.max_z_dist/float(args.slices_thick))) 			#actually it's twice as that (-bins_nb,bins_nb) has to be filled
bins_nb_max = bins_nb
bins_labels = [str((n+0.5)*args.slices_thick) for n in range(-bins_nb,bins_nb)]
slice_volume = 2 * math.pi * args.slices_radius * args.slices_thick
if args.selection_file_prot == "no":
	protein_pres = False
	args.residuesfilename = "no"
else:
	protein_pres = True

#leaflet identification
if args.cutoff_leaflet != "large" and args.cutoff_leaflet != "optimise":
	try:
		args.cutoff_leaflet = float(args.cutoff_leaflet)
	except:
		print "Error: the argument of the --leaflets option should be a number or 'large', see note 2"
		sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "cluster_density_profile_" + args.grofilename[:-4]
	else:
		args.output_folder = "cluster_density_profile_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	
	#sizes
	#-----
	os.mkdir(args.output_folder + "/1_sizes")
	os.mkdir(args.output_folder + "/1_sizes/1_particles")
	os.mkdir(args.output_folder + "/1_sizes/1_particles/png")
	os.mkdir(args.output_folder + "/1_sizes/1_particles/xvg")
	if args.chargesfilename != "no":
		os.mkdir(args.output_folder + "/1_sizes/2_charges")
		os.mkdir(args.output_folder + "/1_sizes/2_charges/png")
		os.mkdir(args.output_folder + "/1_sizes/2_charges/xvg")
	if args.residuesfilename != "no":
		os.mkdir(args.output_folder + "/1_sizes/3_residues")
		os.mkdir(args.output_folder + "/1_sizes/3_residues/png")
		os.mkdir(args.output_folder + "/1_sizes/3_residues/xvg")

	#groups
	#------
	if args.cluster_groups_file != "no":
		os.mkdir(args.output_folder + "/2_groups")
		os.mkdir(args.output_folder + "/2_groups/1_particles")
		os.mkdir(args.output_folder + "/2_groups/1_particles/png")
		os.mkdir(args.output_folder + "/2_groups/1_particles/xvg")
		if args.chargesfilename != "no":
			os.mkdir(args.output_folder + "/2_groups/2_charges")
			os.mkdir(args.output_folder + "/2_groups/2_charges/png")
			os.mkdir(args.output_folder + "/2_groups/2_charges/xvg")
		if args.residuesfilename != "no":
			os.mkdir(args.output_folder + "/2_groups/3_residues")
			os.mkdir(args.output_folder + "/2_groups/3_residues/png")
			os.mkdir(args.output_folder + "/2_groups/3_residues/xvg")
	
	#create log
	#----------
	filename_log = os.getcwd() + '/' + str(args.output_folder) + '/cluster_density_profile.log'
	output_log = open(filename_log, 'w')		
	output_log.write("[cluster_density_profile v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python cluster_density_profile.py"
	for c in sys.argv[1:]:
		tmp_log += " " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	
	#copy input files
	#----------------
	if args.selection_file_ff != "no":
		shutil.copy2(args.selection_file_ff,args.output_folder + "/")	
	if args.selection_file_prot != "no" and args.selection_file_prot != "auto":
		shutil.copy2(args.selection_file_prot,args.output_folder + "/")
	if args.cluster_groups_file != "no":
		shutil.copy2(args.cluster_groups_file,args.output_folder + "/")
	if args.beadsfilename != "no":
		shutil.copy2(args.beadsfilename,args.output_folder + "/")
	if args.particlesfilename != "mine":
		shutil.copy2(args.particlesfilename,args.output_folder + "/")
	if args.residuesfilename != "no" and args.residuesfilename != "mine":
		shutil.copy2(args.residuesfilename,args.output_folder + "/")
	if args.chargesfilename != "no" and args.chargesfilename != "mine":
		shutil.copy2(args.chargesfilename,args.output_folder + "/")

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def set_lipids_beads():													#DONE

	global leaflet_sele_string

	#set default beads
	leaflet_sele_string = "name PO4 or name PO3 or name B1A"

	#use users input
	if args.beadsfilename != "no":
		with open(args.beadsfilename) as f:
			lines = f.readlines()
		if len(lines) > 1:
			print "Error: the file " + str(args.beadsfilename) + " should conly ontain 1 line (" + str(len(lines)) + " found), see note 2(a)."
			sys.exit(1)
		else:
			if lines[0][-1] == "\n":
				lines[0] = lines[0][:-1]
			leaflet_sele_string = lines[0]

	return
def set_particles():													#DONE

	global particles_def
	global particles_groups
	particles_def = {k: {} for k in ["group","colour","sele","sele_string"]}
	colours_map = 'custom'
	colours_map_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']
	
	#use default particles definition
	#--------------------------------
	if args.particlesfilename == "mine":
		particles_def["labels"] = ["peptide","POPC","POPE","POPS","CHOL","water","Na+","Cl-"]
		
		#peptide
		particles_def["group"]["peptide"] = "peptide"					#very dark grey
		particles_def["colour"]["peptide"] = "#262626"					#very dark grey
		particles_def["sele_string"]["peptide"] = "protein"
	
		#lipids
		particles_def["group"]["CHOL"] = "CHOL"
		particles_def["group"]["POPC"] = "POPC"
		particles_def["group"]["POPE"] = "POPE"
		particles_def["group"]["POPS"] = "POPS"
		particles_def["colour"]["CHOL"] = "#bdbdbd"						#light grey
		particles_def["colour"]["POPC"] = "#41ab5d"						#dark green
		particles_def["colour"]["POPE"] = "#6a51a3"						#dark purple
		particles_def["colour"]["POPS"] = "#cc4c02"						#dark orange
		particles_def["sele_string"]["POPC"] = "resname POPC and name PO4"		
		particles_def["sele_string"]["POPE"] = "resname POPE and name PO4"
		particles_def["sele_string"]["POPS"] = "resname POPS and name PO4"
		particles_def["sele_string"]["CHOL"] = "resname CHOL and name ROH"

		#solvent
		particles_def["group"]["water"] = "water"
		particles_def["colour"]["water"] = "#1d91c0"					#dark cyan
		particles_def["sele_string"]["water"] = "resname W"

		#ions
		particles_def["group"]["Na+"] = "Na+"
		particles_def["group"]["Cl-"] = "CL-"
		particles_def["colour"]["Cl-"] = "#fa9fb5"						#light red
		particles_def["colour"]["Na+"] = "#7bccc4"						#light blue
		particles_def["sele_string"]["Na+"] = "name NA+"
		particles_def["sele_string"]["Cl-"] = "name CL-"
	
	#use user's particles definition
	#-------------------------------	
	else:
		particles_def["labels"] = []
		with open(args.particlesfilename) as f:
			lines = f.readlines()
		for l_index in range(0, len(lines)):
			line = lines[l_index]
			if line[-1] == "\n":
				line = line[:-1]
			l_content = line.split(',')
			if len(l_content) != 4:
				print "Error: wrong format on line " + str(l_index + 1) + " of " + str(args.particlesfilename) + ", see DESCRIPTION in cluster_density_profile --help."
				sys.exit(1)
			tmp_group = l_content[0]
			tmp_label = l_content[1]
			tmp_color = l_content[2]
			tmp_selec = l_content[3]
			if tmp_label in particles_def["labels"]:
				print "Error: the particle " + str(tmp_label) + " is defined mor than once. Check your --particles option."
				sys.exit(1)
			particles_def["labels"].append(tmp_label)
			particles_def["group"][tmp_label] = tmp_group
			particles_def["colour"][tmp_label] = tmp_color
			particles_def["sele_string"][tmp_label] = tmp_selec

		#generate colours if necessary
		if len(np.unique(particles_def["colour"].values())) == 1:
			if np.unique(particles_def["colour"].values())[0] in colours_map_possible:
				colours_map = np.unique(particles_def["colour"].values())[0]
			else:
				print "Error: either the same colour was specified for all particles or the colour map '" + str(np.unique(particles_def["colour"].values())[0]) + "' is not valid."
				sys.exit(1)
		if colours_map != "custom":
			tmp_cmap = cm.get_cmap(colours_map)
			part_nb = len(particles_def["labels"])
			part_colours_value = tmp_cmap(np.linspace(0, 1, part_nb))
			for part_index in range(0, part_nb):
				particles_def["colour"][particles_def["labels"][part_index]] = part_colours_value[part_index]

	#check whether a protein or peptide group has been defined to calculate residues details
	#---------------------------------------------------------------------------------------
	if args.residuesfilename != "no" and "peptide" not in particles_def["labels"]:
		print "Error: no 'peptide' particles defined, residues details cannot be calculated. Use '--residues no' or update --particles. See DESCRIPTION in cluster_density_profile --help."
		sys.exit(1)

	#build particle groups for normalisation
	#---------------------------------------
	particles_groups = {k: [part for part, g in particles_def["group"].items() if g == k] for k in np.unique(particles_def["group"].values())}	
		
	return
def set_residues():														#DONE

	global residues_def
	residues_def = {k: {} for k in ["colour","sele","res_list","sele_string"]}	
	colours_map = 'custom'
	colours_map_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']
	
	#use default residues definitions
	#--------------------------------
	if args.residuesfilename == "mine":
		residues_def["labels"] = ["basic","polar","hydrophobic","backbone"]
		#basic
		residues_def["colour"]["basic"] = '#253494'						#blue
		residues_def["res_list"]["basic"] = ['ARG','LYS']
		#polar
		residues_def["colour"]["polar"] = '#a1d99b'						#greenish
		residues_def["res_list"]["polar"] = ['SER','THR','ASN','GLN','HIS']		
		#hydrophobic
		residues_def["colour"]["hydrophobic"] = '#993404'				#orange brown
		residues_def["res_list"]["hydrophobic"] = ['VAL','ILE','LEU','MET','PHE','PRO','CYS','TYR','TRP']
		#backbone only
		residues_def["colour"]["backbone"] = '#969696'					#light grey
		residues_def["res_list"]["backbone"] = ['ALA','GLY']	
		#acidic
		#residues_def["colour"]["acidic"] = y
		#residues_def["res_list"]["acidic"] = ['ASP','GLU']
	
	#use user's residues definitions
	#-------------------------------
	else:	
		residues_def["labels"] = []
		with open(args.residuesfilename) as f:
			lines = f.readlines()
		for l_index in range(0, len(lines)):
			line = lines[l_index]
			if line[-1] == "\n":
				line = line[:-1]
			l_content = line.split(',')
			if len(l_content) < 3:
				print "Error: wrong format on line " + str(l_index + 1) + " of " + str(args.residuesfilename) + ", see DESCRIPTION in cluster_density_profile --help."
				sys.exit(1)
			tmp_label = l_content[0]
			tmp_color = l_content[1]
			if tmp_label in residues_def["labels"]:
				print "Error: the residue " + str(tmp_label) + " is defined mor than once. Check your --particles option."
				sys.exit(1)
			residues_def["labels"].append(tmp_label)
			residues_def["colour"][tmp_label] = tmp_colour
			residues_def["res_list"][tmp_label] = [str(l_content[r]) for r in range(2,len(l_content))]
			
		#generate colours from colour map if necessary
		if len(np.unique(residues_def["colour"].values())) == 1:
			if np.unique(residues_def["colour"].values())[0] in colours_map_possible:
				colours_map = np.unique(residues_def["colour"].values())[0]
			else:
				print "Error: either the same colour was specified for all residues or the colour map '" + str(np.unique(residues_colours.values())[0]) + "' is not valid."
				sys.exit(1)
		if colours_map != "custom":
			tmp_cmap = cm.get_cmap(colours_map)
			residues_nb = len(residues_def["labels"])
			residues_colours_value = tmp_cmap(np.linspace(0, 1, residues_nb))
			for r_index in range(0, residues_nb):
				residues_colours[residues_def["labels"][t_index]] = residues_colours_value[r_index]

	#create selection strings for each type
	#--------------------------------------
	for res in residues_def["labels"]:
		residues_def["sele_string"][res] = "resname " + str(residues_def["res_list"][res][0])	
		for r in residues_def["res_list"][res][1:]:
			residues_def["sele_string"][res] += " or resname " + str(r)

	return
def set_charges():														#DONE
	
	global charges_groups
	global charges_colours
	charges_groups = {}
	charges_colours = {}
	charges_colours["total"] = '#262626'								#very dark grey
	colours_map = 'custom'
	colours_map_possible = ['Spectral', 'summer', 'coolwarm', 'pink_r', 'Set1', 'Set2', 'Set3', 'brg_r', 'Dark2', 'hot', 'PuOr_r', 'afmhot_r', 'terrain_r', 'PuBuGn_r', 'RdPu', 'gist_ncar_r', 'gist_yarg_r', 'Dark2_r', 'YlGnBu', 'RdYlBu', 'hot_r', 'gist_rainbow_r', 'gist_stern', 'gnuplot_r', 'cool_r', 'cool', 'gray', 'copper_r', 'Greens_r', 'GnBu', 'gist_ncar', 'spring_r', 'gist_rainbow', 'RdYlBu_r', 'gist_heat_r', 'OrRd_r', 'CMRmap', 'bone', 'gist_stern_r', 'RdYlGn', 'Pastel2_r', 'spring', 'terrain', 'YlOrRd_r', 'Set2_r', 'winter_r', 'PuBu', 'RdGy_r', 'spectral', 'flag_r', 'jet_r', 'RdPu_r', 'Purples_r', 'gist_yarg', 'BuGn', 'Paired_r', 'hsv_r', 'bwr', 'cubehelix', 'YlOrRd', 'Greens', 'PRGn', 'gist_heat', 'spectral_r', 'Paired', 'hsv', 'Oranges_r', 'prism_r', 'Pastel2', 'Pastel1_r', 'Pastel1', 'gray_r', 'PuRd_r', 'Spectral_r', 'gnuplot2_r', 'BuPu', 'YlGnBu_r', 'copper', 'gist_earth_r', 'Set3_r', 'OrRd', 'PuBu_r', 'ocean_r', 'brg', 'gnuplot2', 'jet', 'bone_r', 'gist_earth', 'Oranges', 'RdYlGn_r', 'PiYG', 'CMRmap_r', 'YlGn', 'binary_r', 'gist_gray_r', 'Accent', 'BuPu_r', 'gist_gray', 'flag', 'seismic_r', 'RdBu_r', 'BrBG', 'Reds', 'BuGn_r', 'summer_r', 'GnBu_r', 'BrBG_r', 'Reds_r', 'RdGy', 'PuRd', 'Accent_r', 'Blues', 'Greys', 'autumn', 'cubehelix_r', 'nipy_spectral_r', 'PRGn_r', 'Greys_r', 'pink', 'binary', 'winter', 'gnuplot', 'RdBu', 'prism', 'YlOrBr', 'coolwarm_r', 'rainbow_r', 'rainbow', 'PiYG_r', 'YlGn_r', 'Blues_r', 'YlOrBr_r', 'seismic', 'Purples', 'bwr_r', 'autumn_r', 'ocean', 'Set1_r', 'PuOr', 'PuBuGn', 'nipy_spectral', 'afmhot']
	
	#use default
	#-----------
	if args.chargesfilename == "mine":
		
		#ions
		charges_colours["ions"] = "#52A3CC"								#cyan colour
		charges_groups["ions"] = {}
		charges_groups["ions"]["names"] = ["Na+","Cl-"]
		charges_groups["ions"]["values"] = {}
		charges_groups["ions"]["values"]["Na+"] = 1
		charges_groups["ions"]["values"]["Cl-"] = -1
		charges_groups["ions"]["sele"] = {}
		charges_groups["ions"]["sele_string"] = {}
		charges_groups["ions"]["sele_string"]["Na+"] = "name NA+"
		charges_groups["ions"]["sele_string"]["Cl-"] = "name CL-"
		
		#lipids
		charges_colours["lipids"] = "#b2182b"							#dark red
		charges_groups["lipids"] = {}
		charges_groups["lipids"]["names"] = ["POPS"]
		charges_groups["lipids"]["values"] = {}
		charges_groups["lipids"]["values"]["POPS"] = -1
		charges_groups["lipids"]["sele"] = {}
		charges_groups["lipids"]["sele_string"] = {}
		charges_groups["lipids"]["sele_string"]["POPS"] = "resname POPS and name PO4"
		
		#transportan
		charges_colours["peptide"] = "#053061"							#dark blue
		charges_groups["peptide"] = {}
		charges_groups["peptide"]["names"] = ["pos","neg"]
		charges_groups["peptide"]["values"] = {}
		charges_groups["peptide"]["values"]["pos"] = 1
		charges_groups["peptide"]["values"]["neg"] = -1
		charges_groups["peptide"]["sele"] = {}
		charges_groups["peptide"]["sele_string"] = {}
		charges_groups["peptide"]["sele_string"]["pos"] = "(resnum 1 and resname GLY and name BB) or (resname LYS and name SC2)"
		charges_groups["peptide"]["sele_string"]["neg"] = "resnum 27 and resname LEU and name BB"
	
	#use user supplied
	#-----------------
	else:	
		with open(args.chargesfilename) as f:
			lines = f.readlines()
		for l_index in range(0, len(lines)):
			line = lines[l_index]
			if line[-1] == "\n":
				line = line[:-1]
			l_content = line.split(',')
			if len(l_content) != 5:
				print "Error: wrong format on line " + str(l_index + 1) + " of " + str(args.chargesfilename) + ", see DESCRIPTION in cluster_density_profile --help."
				sys.exit(1)
			tmp_group = l_content[0]
			tmp_colour = l_content[1]
			tmp_name = l_content[2]
			tmp_value = float(l_content[3])
			tmp_sele = l_content[4]
			if tmp_group not in charges_groups.keys():
				charges_groups[tmp_group] = {k: {} for k in ["values","sele","sele_string"]}
				charges_groups[tmp_group]["names"] = []
				charges_colours[tmp_groups] = tmp_colour
			charges_groups[tmp_group]["names"].append(tmp_name)
			charges_groups[tmp_group]["values"][tmp_name] = tmp_value
			charges_groups[tmp_group]["sele_string"][tmp_name] = tmp_sele

		#generate colours from colour map if necessary
		if len(np.unique(charges_colours.values())) == 1:
			if np.unique(charges_colours.values())[0] in colours_map_possible:
				colours_map = np.unique(charges_colours.values())[0]
			else:
				print "Error: either the same colour was specified for all charge groups or the colour map '" + str(np.unique(residues_colours.values())[0]) + "' is not valid."
				sys.exit(1)
		if colours_map != "custom":
			tmp_cmap = cm.get_cmap(colours_map)
			charges_groups_keys = charges_groups.keys()
			charges_groups_nb = len(charges_groups_keys)
			charges_colours_value = tmp_cmap(np.linspace(0, 1, charges_groups_nb))
			for q_index in range(0, charges_groups_nb):
				charges_colours[charges_groups_keys[q_index]] = charges_colours_value[q_index]

	return
def load_MDA_universe():												#DONE
	
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

	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no lipid particles selected, check the --leaflets and/or --beads options."
		sys.exit(1)

	#create selection: particles density
	#-----------------------------------
	water_pres = False
	particles_pres = False
	
	#check for presence of each particle
	for part in particles_def["labels"]:
		particles_def["sele"][part] = U.selectAtoms(particles_def["sele_string"][part])
		if particles_def["sele"][part].numberOfAtoms() == 0:
			print " ->warning: particle selection string '" + str(particles_def["sele_string"][part]) + "' returned 0 atoms."
		else:
			particles_pres = True
			if part == "water":
				water_pres = True
	if not particles_pres:
		print "Error: none of the specified particles was found in the system, check your --particles option."
		sys.exit(1)
				
	#check for the presence of protein
	#---------------------------------
	test_prot = U.selectAtoms("protein")
	if protein_pres:
		if test_prot.numberOfAtoms() == 0:
			print "Error: no proteins found in the system, check your --proteins option (set it to 'no' if there are no proteins)."
			sys.exit(1)
	else:
		if test_prot.numberOfAtoms() > 0:
			print " ->warning: proteins have been detected in the system but the option --proteins is set to 'no'; proteins will be ignored but are likely to influence the density profiles."
		
	#create selection: residues density
	#----------------------------------
	residues_pres = False
	if args.residuesfilename != "no":
		for res in residues_def["labels"]:
			residues_def["sele"][res] = U.selectAtoms(residues_def["sele_string"][res])
			if residues_def["sele"][res].numberOfAtoms() == 0:
				print " ->warning: residues selection string '" + str(residues_def["sele_string"][res]) + "' returned 0 atoms."
			else:
				residues_pres = True		
		if not residues_pres:
			print "Error: none of the specified residues was found in the system, check your --residues option."
			sys.exit(1)
	
	#create charged particles selections
	#-----------------------------------
	charge_pres = False
	if args.chargesfilename != "no":
		for charge_g in charges_groups.keys():
			for q in charges_groups[charge_g]["names"]:
				charges_groups[charge_g]["sele"][q] = U.selectAtoms(charges_groups[charge_g]["sele_string"][q])
			if charges_groups[charge_g]["sele"][q].numberOfAtoms() == 0:
				print " ->warning: charge selection string '" + str(charges_groups[charge_g]["sele_string"][q]) + "' returned 0 atoms."
			else:
				charge_pres = True
		if not charge_pres:
			print "Error: no charged particles found, use '--charges no' or supply correct charges definition."
			sys.exit(1)

	#check whether the arg max is less than half the box size
	#--------------------------------------------------------
	if args.max_z_dist > U.dimensions[2]/float(2):
		print "Warning: the --dist_max option (" + str(args.max_z_dist) + ") is larger than half the box size in the z direction plane (" + str(U.dimensions[2]/float(2)) + "), density profile will be truncated at " + str(U.dimensions[2]/float(2)) + " Angstrom."
		args.max_z_dist = U.dimensions[2]/float(2)
		bins_nb_max = int(floor(args.max_z_dist/float(args.slices_thick)))
		
	return
def identify_ff():														#DONE
	print "\nReading selection file for flipflopping lipids..."
	
	#declare variables
	global lipids_ff_nb
	global lipids_ff_info
	global lipids_ff_resnames
	global lipids_ff_leaflet
	global lipids_ff_u2l_index
	global lipids_ff_l2u_index
	global lipids_sele_ff
	global lipids_sele_ff_bead
	global lipids_sele_ff_bonds
	global lipids_sele_ff_VMD_string
	global leaflet_sele_string
	lipids_ff_nb = 0
	lipids_ff_info = {}
	lipids_ff_resnames = []
	lipids_ff_leaflet = []
	lipids_ff_u2l_index = []
	lipids_ff_l2u_index = []
	lipids_sele_ff = {}
	lipids_sele_ff_bead = {}
	lipids_sele_ff_bonds = {}
	lipids_sele_ff_VMD_string={}
		
	with open(args.selection_file_ff) as f:
		lines = f.readlines()
	lipids_ff_nb = len(lines)
	print " -found " + str(lipids_ff_nb) + " flipflopping lipids"
	leaflet_sele_string = leaflet_sele_string + " and not ("
	for l_index in range(0,lipids_ff_nb):
		line = lines[l_index]
		if line[-1] == "\n":
			line = line[:-1]
		try:
			line_content = line.split(',')
			if len(line_content) != 4:
				print "Error: wrong format for line " + str(l_index+1) + " in " + str(args.selection_file_ff) + ", see note 4 in bilayer_perturbations --help."
				print " ->", line
				sys.exit(1)
			#read current lipid details
			lip_resname = line_content[0]
			lip_resnum = int(line_content[1])
			lip_leaflet = line_content[2]
			lip_bead = line_content[3]
			lipids_ff_info[l_index] = [lip_resname,lip_resnum,lip_leaflet,lip_bead]
						
			#update: starting leaflets
			if lip_leaflet not in lipids_ff_leaflet:
				lipids_ff_leaflet.append(lip_leaflet)

			#update: index in directional lists
			if lip_leaflet == "upper":
				lipids_ff_u2l_index.append(l_index)
			elif lip_leaflet == "lower":
				lipids_ff_l2u_index.append(l_index)
			else:
				print "->unknown starting leaflet '" + str(lip_leaflet) + "'."
				sys.exit(1)
			
			#update: resnames
			if lip_resname not in lipids_ff_resnames:
				lipids_ff_resnames.append(lip_resname)
	
			#update: leaflet selection string
			if l_index==0:
				leaflet_sele_string+="(resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"
			else:
				leaflet_sele_string+=" or (resname " + str(lip_resname) + " and resnum " + str(lip_resnum) + ")"

			#create selections
			lipids_sele_ff[l_index] = U.selectAtoms("resname " + str(lip_resname) + " and resnum " + str(lip_resnum))
			lipids_sele_ff_bead[l_index] = lipids_sele_ff[l_index].selectAtoms("name " + str(lip_bead))
			lipids_sele_ff_VMD_string[l_index]="resname " + str(lipids_ff_info[l_index][0]) + " and resid " + str(lipids_ff_info[l_index][1])
			if lipids_sele_ff[l_index].numberOfAtoms() == 0:
				print "Error:"
				print line
				print "-> no such lipid found."
				sys.exit(1)	
		except:
			print "Error: invalid flipflopping lipid selection string on line " + str(l_index+1) + ": '" + line + "'"
			sys.exit(1)
	leaflet_sele_string+=")"		

	return
def identify_proteins():												#DONE
	print "\nIdentifying proteins..."
	
	#import modules
	if args.m_algorithm == "density":
		global DBSCAN
		from sklearn.cluster import DBSCAN
	else:
		global nx
		import networkx as nx

	#declare variables
	global proteins_nb
	global proteins_sele
	global proteins_sele_string
	global proteins_sele_string_VMD
	global proteins_boundaries
	global proteins_nb_atoms
	global nb_atom_per_protein
	proteins_nb = 0
	proteins_sele = {}
	proteins_sele_string = {}
	proteins_sele_string_VMD = {}
	proteins_boundaries = {}
	
	#check for protein presence
	if U.selectAtoms("protein").numberOfAtoms() == 0:
		print "Error: no protein detected."
		sys.exit(1)
	
	#case: selection file provided
	if args.selection_file_prot != "auto":
		print " -reading protein selection file..."
		with open(args.selection_file_prot) as f:
			lines = f.readlines()
		proteins_nb=len(lines)
		proteins_sele["all"] = MDAnalysis.core.AtomGroup.AtomGroup([])
		for p_index in range(0,proteins_nb):
			line = lines[p_index]
			if line[-1] == "\n":
				line = line[:-1]
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			try:
				print " p[" + str(p_index) + "]=U.selectAtoms(" + line + ")"
				proteins_sele[p_index] = U.selectAtoms(line[1:-2])
				proteins_sele["all"] += proteins_sele[p_index]
				proteins_boundaries[p_index] = [proteins_sele[p_index].indices()[0] + 1, proteins_sele[p_index].indices()[proteins_sele[p_index].numberOfAtoms()]+1]
				proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
				proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			except:
				print "Error:"
				print line
				print "->invalid selection string."
				sys.exit(1)
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
	
	#case: automatic detection
	else:
		#declare local variables
		proteins_ca_nb = {}
		proteins_ca_nmax = 0
		proteins_ca_group = {}
		proteins_boundaries = {}
	
		#retrieve 1st atom info
		proteins_sele["all"] = U.selectAtoms("protein")
		proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
		prec_resnum = proteins_sele["all"][0].resnum
		prec_segid = proteins_sele["all"][0].segid
		prec_atnum = proteins_sele["all"][0].number+1
		prev_atnum = proteins_sele["all"][0].number+1						#atom corresponding to the beginning of the current protein
		#browse following atoms
		for a in proteins_sele["all"][1:]:
			delta_res = a.resnum-prec_resnum
			delta_atm = a.number+1-prec_atnum
			if delta_res < 0 or a.segid != prec_segid or delta_atm > 1:
				proteins_boundaries[proteins_nb] = [prev_atnum,prec_atnum]
				proteins_nb += 1
				prev_atnum = a.number + 1
			prec_resnum = a.resnum
			prec_atnum = a.number + 1
			prec_segid = a.segid		
		#add last protein section
		if prev_atnum < proteins_sele["all"][proteins_nb_atoms-1].number:
			proteins_boundaries[proteins_nb] = [prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
			proteins_nb += 1
		
		#display results
		print " -protein found:", proteins_nb
		print " -protein boundaries (atom numbers): see protein.sele file"
		#create protein selections and save into a txt file
		filename_sele=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
		output_stat = open(filename_sele, 'w')	
		output_stat.write("#This file was generated by the script bilayer_perturbations v" + str(version_nb) +"\n")
		output_stat.write("#The lines below correspond to MDAnalysis section string, e.g. U.selectAtoms(LINE)\n")
		output_stat.write("\n")	
		for p_index in range(0, proteins_nb):
			progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
			sys.stdout.flush()
			sys.stdout.write(progress)
			proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
			proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
			proteins_sele[p_index] = U.selectAtoms(proteins_sele_string[p_index])
			output_stat.write(proteins_sele_string[p_index] + "\n")
		output_stat.close()

	nb_atom_per_protein = proteins_sele[0].numberOfAtoms()
	print ""

	return
def identify_leaflets():												#DONE
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#check the leaflet selection string is valid
	test_beads = U.selectAtoms(leaflet_sele_string)
	if test_beads.numberOfAtoms() == 0:
		print "Error: invalid selection string '" + str(leaflet_sele_string) + "'"
		print "-> no particles selected."
		sys.exit(1)

	#use LeafletFinder:
	if args.cutoff_leaflet != 'large':
		if args.cutoff_leaflet == 'optimise':
			print " -optimising cutoff..."
			cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
		else:
			L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, args.cutoff_leaflet)
	
		if np.shape(L.groups())[0]<2:
			print "Error: imposssible to identify 2 leaflets."
			sys.exit(1)
		if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
			leaflet_sele["upper"] = L.group(0)
			leaflet_sele["lower"] = L.group(1)
		else:
			leaflet_sele["upper"] = L.group(1)
			leaflet_sele["lower"] = L.group(0)
		leaflet_sele["both"] = leaflet_sele["lower"] + leaflet_sele["upper"]
		if np.shape(L.groups())[0] == 2:
			print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		else:
			other_lipids=0
			for g in range(2, np.shape(L.groups())[0]):
				other_lipids += L.group(g).numberOfResidues()
			print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
	#use cog:
	else:
		leaflet_sele["both"] = U.selectAtoms(leaflet_sele_string)
		tmp_lipids_avg_z = leaflet_sele["both"].centerOfGeometry()[2]
		leaflet_sele["upper"] = leaflet_sele["both"].selectAtoms("prop z > " + str(tmp_lipids_avg_z))
		leaflet_sele["lower"] = leaflet_sele["both"].selectAtoms("prop z < " + str(tmp_lipids_avg_z))
		print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
		
	return
def initialise_groups():												#DONE

	global groups_labels
	global groups_number
	global groups_boundaries
	global groups_sizes_dict
	groups_labels = {}
	groups_boundaries = {}
	groups_sizes_dict = {}
	
	#read file
	print "\nReading cluster groups definition file..."
	with open(args.cluster_groups_file) as f:
		lines = f.readlines()
	groups_number = len(lines)
	for g_index in range(0,groups_number):
		line = lines[g_index]
		if line[-1] == "\n":
			line = line[:-1]
		l_content = line.split(',')
		#check format
		if len(l_content) != 2:
			print "Error: the format of line " + str(g_index+1) + " should be 'min,max' (see cluster_density_profile --help, note 4)."
			print "->", line
			sys.exit(1)
		tmp_beg = int(l_content[0])
		tmp_end = l_content[1]
		if tmp_end == "max":
			tmp_end = 100000										#put a stupidly big size to cap the open ended group
		else:
			tmp_end = int(tmp_end)
		groups_boundaries[g_index] = [tmp_beg,tmp_end]
		
	#display results
	print " -found " + str(groups_number) + " cluster groups:"
	for g_index in range(0,groups_number):
		if groups_boundaries[g_index][1] == 100000:
			print "   g" + str(g_index) + " = " + str(groups_boundaries[g_index][0]) + "+"
		else:
			print "   g" + str(g_index) + " = " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])

	#check for boundaries overlapping
	prev_beg = groups_boundaries[0][0]
	prev_end = groups_boundaries[0][1]
	if prev_end < prev_beg:
		print "Error: the max size is smaller than the min size for specified cluster groups " + str(groups_boundaries[0]) + "."
		sys.exit(1)
	for g_index in range(1,groups_number):
		if groups_boundaries[g_index][1] < groups_boundaries[g_index][0]:
			print "Error: the max size is smaller than the min size for group " + str(g_index) + "(" + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + ")."
			sys.exit(1)
		if groups_boundaries[g_index][0] <= prev_end:
			print "Error: specified cluster groups " + str(prev_beg) + "-" + str(prev_end) + " and " + str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1]) + " overlap or are not in increasing order (boundaries are inclusive, see note 7 in bilayer_perturbations --help)."
			sys.exit(1)
		prev_beg = groups_boundaries[g_index][0]
		prev_end = groups_boundaries[g_index][1]
	
	#create equivalency table between groups and sizes
	for g_index in range(0,groups_number):
		bb = groups_boundaries[g_index]
		tmp_beg = bb[0]
		tmp_end = bb[1]
		for tmp_size in range(tmp_beg, tmp_end+1):
			groups_sizes_dict[tmp_size] = g_index
	for tmp_size in list(set(range(1,max(groups_sizes_dict.keys()))) - set(groups_sizes_dict.keys())): 		#this handles potentially unaccounted for sizes up to the maximum specified by the user
		groups_sizes_dict[tmp_size] = groups_number
	if max(groups_sizes_dict.keys())!= 100000:															    #this handles potentially unaccounted for sizes above the maximum specified by the user (in case it's not an open group)
		for tmp_size in range(max(groups_sizes_dict.keys())+1,100001):
			groups_sizes_dict[tmp_size] = groups_number
			
	#create label for each group
	for g_index in range(0, groups_number):
		if groups_boundaries[g_index][1] == 100000:
			groups_labels[g_index] = str(groups_boundaries[g_index][0]) + "+"
		elif groups_boundaries[g_index][0] == groups_boundaries[g_index][1]:
			groups_labels[g_index] = str(groups_boundaries[g_index][0])
		else:
			groups_labels[g_index] = str(groups_boundaries[g_index][0]) + "-" + str(groups_boundaries[g_index][1])
	groups_labels[groups_number] = "other"
		
	return

#=========================================================================================
# data structures
#=========================================================================================

def struct_particles():
	global z_upper
	global z_lower
	global sizes_sampled
	global sizes_coverage
	global sizes_nb_clusters
	z_upper = 0
	z_lower = 0
	sizes_sampled = []
	sizes_coverage = {}
	sizes_nb_clusters = {"all sizes": 0}
	
	#sizes
	global density_particles_sizes_nb
	global density_particles_sizes_pc
	sizes_coverage["particles"] = {k: {"all sizes": {part: 0 for part in particles_def["labels"]}} for k in ["nb","avg","std"]}
	density_particles_sizes_nb = {"all sizes": {part: np.zeros(2*bins_nb) for part in particles_def["labels"]}}
	density_particles_sizes_pc = {"all sizes": {part: np.zeros(2*bins_nb) for part in particles_def["labels"]}}
	
	#groups
	if args.cluster_groups_file != "no":
		global groups_sampled
		global groups_coverage
		global groups_nb_clusters
		global density_particles_groups_nb
		global density_particles_groups_pc
		groups_sampled = []
		groups_coverage = {}
		groups_nb_clusters = {g_index: 0 for g_index in range(0,groups_number+1)}
		groups_coverage["particles"] = {k: {g_index: {part: 0 for part in particles_def["labels"]} for g_index in range(0, groups_number+1)} for k in ["nb","avg","std"]}
		density_particles_groups_nb = {g_index: {part: np.zeros(2*bins_nb) for part in particles_def["labels"]} for g_index in range(0, groups_number+1)}
		density_particles_groups_pc = {g_index: {part: np.zeros(2*bins_nb) for part in particles_def["labels"]} for g_index in range(0, groups_number+1)}		

	return
def struct_residues():

	#sizes
	global density_residues_sizes_nb
	global density_residues_sizes_pc
	sizes_coverage["residues"] = {k: {"all sizes": {res: 0 for res in residues_def["labels"]}} for k in ["nb","avg","std"]}
	density_residues_sizes_nb = {"all sizes": {res: np.zeros(2*bins_nb) for res in residues_def["labels"]}}
	density_residues_sizes_pc = {"all sizes": {res: np.zeros(2*bins_nb) for res in residues_def["labels"]}}
	
	#groups
	if args.cluster_groups_file != "no":
		global density_residues_groups_nb
		global density_residues_groups_pc
		groups_coverage["residues"] = {k: {g_index: {res: 0 for res in residues_def["labels"]} for g_index in range(0, groups_number+1)} for k in ["nb","avg","std"]}
		density_residues_groups_nb = {g_index: {res: np.zeros(2*bins_nb) for res in residues_def["labels"]} for g_index in range(0, groups_number+1)}
		density_residues_groups_pc = {g_index: {res: np.zeros(2*bins_nb) for res in residues_def["labels"]} for g_index in range(0, groups_number+1)}		

	return
def struct_charges():
	
	#sizes
	global density_charges_sizes
	sizes_coverage["charges"] = {k: {"all sizes": {q: 0 for q in charges_groups.keys() + ["total"]}} for k in ["nb","avg","std"]}
	density_charges_sizes = {"all sizes": {q: np.zeros(2*bins_nb) for q in charges_groups.keys() + ["total"]}}

	#groups
	if args.cluster_groups_file != "no":
		global density_charges_groups
		groups_coverage["charges"] = {k: {g_index: {q: 0 for q in charges_groups.keys() + ["total"]} for g_index in range(0, groups_number+1)} for k in ["nb","avg","std"]}
		density_charges_groups = {g_index: {q: np.zeros(2*bins_nb) for q in charges_groups.keys() + ["total"]} for g_index in range(0, groups_number+1)}

	return
	
#=========================================================================================
# core functions
#=========================================================================================

def get_z_coords():														#DONE
	
	global z_upper
	global z_lower	
	global z_middle_instant

	tmp_zu = leaflet_sele["upper"].centerOfGeometry()[2]
	tmp_zl = leaflet_sele["lower"].centerOfGeometry()[2]
	z_middle_instant = tmp_zl + (tmp_zu - tmp_zl)/float(2)
	z_upper += tmp_zu - z_middle_instant
	z_lower += tmp_zl - z_middle_instant
		
	return
def get_distances(box_dim):												#DONE
		
	#method: use minimum distance between proteins
	#---------------------------------------------
	if args.m_algorithm == "min":
		#pre-process: get protein coordinates
		tmp_proteins_coords = np.zeros((proteins_nb, nb_atom_per_protein, 3))
		for p_index in range(0, proteins_nb):
			tmp_proteins_coords[p_index,:] = proteins_sele[p_index].coordinates()

		#store min distance between each proteins
		dist_matrix = 100000 * np.ones((proteins_nb,proteins_nb))
		for n in range(proteins_nb,1,-1):
			dist_matrix[proteins_nb-n,proteins_nb-n+1:proteins_nb] = map(lambda pp: np.min(MDAnalysis.analysis.distances.distance_array(np.float32(tmp_proteins_coords[proteins_nb-n,:]), np.float32(tmp_proteins_coords[pp,:]), box_dim)), range(proteins_nb-n+1,proteins_nb))
			dist_matrix[proteins_nb-n+1:proteins_nb,proteins_nb-n] = dist_matrix[proteins_nb-n,proteins_nb-n+1:proteins_nb]
											
	#method: use distance between cog
	#--------------------------------
	else:
		tmp_proteins_cogs = np.asarray(map(lambda p_index: calculate_cog(proteins_sele[p_index], box_dim), range(0,proteins_nb)))
		dist_matrix = MDAnalysis.analysis.distances.distance_array(np.float32(tmp_proteins_cogs), np.float32(tmp_proteins_cogs), box_dim)

	return dist_matrix
def detect_clusters_connectivity(dist, box_dim):						#DONE
	
	#use networkx algorithm
	connected=(dist<args.cutoff_connect)
	network=nx.Graph(connected)
	groups=nx.connected_components(network)
	
	return groups
def detect_clusters_density(dist, box_dim):								#DONE
	
	#run DBSCAN algorithm
	dbscan_output = DBSCAN(eps = args.dbscan_dist, metric = 'precomputed', min_samples = args.dbscan_nb).fit(dist)

	#build 'groups' structure i.e. a list whose element are all the clusters identified
	groups = []
	for c_lab in np.unique(dbscan_output.labels_):
		tmp_pos = np.argwhere(dbscan_output.labels_ == c_lab)
		if c_lab == -1:
			groups += map(lambda p:p[0] , tmp_pos)
		else:
			groups.append(map(lambda p:p[0] , tmp_pos))

	return groups
def calculate_cog(tmp_coords, box_dim):									#DONE
	
	#this method allows to take pbc into account when calculcating the center of geometry 
	#see: http://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
	
	cog_coord = np.zeros(3)
	tmp_nb_atoms = np.shape(tmp_coords)[0]
	
	for n in range(0,3):
		tet = tmp_coords[:,n] * 2 * math.pi / float(box_dim[n])
		xsi = np.cos(tet)
		zet = np.sin(tet)
		tet_avg = math.atan2(-np.average(zet),-np.average(xsi)) + math.pi
		cog_coord[n] = tet_avg * box_dim[n] / float(2*math.pi)
	
	return cog_coord
def calculate_density(box_dim):											#DONE
	
	global sizes_sampled
	global groups_sampled
	
	#retrieve coordinates arrays (pre-processing saves time as MDAnalysis functions are quite slow and we need to make such calls a few times)
	tmp_lip_coords = {l: leaflet_sele[l].coordinates() for l in ["lower","upper"]}
	get_z_coords()
	
	#case: no proteins
	#=================
	if not protein_pres:

		sizes_nb_clusters["all sizes"] += 1
		
		#density profile: particles
		#--------------------------
		for part in particles_def["labels"]:
			if part != "peptide":
				tmp_part_sele = particles_def["sele"][part]
				if tmp_part_sele.numberOfAtoms() > 0:
														
					#center z coordinates on the bilayer center z coordinate
					tmp_coord = tmp_part_sele.coordinates()
					tmp_coord[:,2] -= z_middle_instant
				
					#deal with pbc and center axis on 0
					tmp_coord[:,0] -= (np.floor(2*tmp_coord[:,0]/float(box_dim[0])) + (1-np.sign(tmp_coord[:,0]))/float(2)) * box_dim[0]
					tmp_coord[:,1] -= (np.floor(2*tmp_coord[:,1]/float(box_dim[1])) + (1-np.sign(tmp_coord[:,1]))/float(2)) * box_dim[1]
					tmp_coord[:,2] -= (np.floor(2*tmp_coord[:,2]/float(box_dim[2])) + (1-np.sign(tmp_coord[:,2]))/float(2)) * box_dim[2]
										
					#keep those within the specified radius
					tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
	
					#update particle coverage (Chang algorithm for on the fly average and std dev - actually Knuth simple version as we add element one by one)
					tmp_pc_coord = np.shape(tmp_coord_within)[0]/float(np.shape(tmp_coord)[0])*100
					delta =  tmp_pc_coord - sizes_coverage["particles"]["avg"]["all sizes"][part]
					sizes_coverage["particles"]["nb"]["all sizes"][part] += 1
					sizes_coverage["particles"]["avg"]["all sizes"][part] += delta / sizes_coverage["particles"]["nb"]["all sizes"][part]
					sizes_coverage["particles"]["std"]["all sizes"][part] += delta * (tmp_pc_coord - sizes_coverage["particles"]["avg"]["all sizes"][part])
					
					#add number of particles within each slice					
					tmp_bins_nb = np.zeros(2*bins_nb)
					bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
					bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)
					if len(bin_abs) > 0:
						tmp_bins_nb = np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0]
					density_particles_sizes_nb["all sizes"][part] += tmp_bins_nb

		#density profile: charges
		#------------------------
		if args.chargesfilename != "no":
			for charge_g in charges_groups.keys():
				tmp_bins_nb = np.zeros(2*bins_nb)
				for q in charges_groups[charge_g]["names"]:
					tmp_q_sele = charges_groups[charge_g]["sele"][q]
					if tmp_q_sele.numberOfAtoms() > 0:
					
						#center z coordinates on the bilayer center z coordinate
						tmp_coord = tmp_q_sele.coordinates()
						tmp_coord[:,2] -= z_middle_instant

						#deal with pbc and center axis on 0
						tmp_coord[:,0] -= (np.floor(2*tmp_coord[:,0]/float(box_dim[0])) + (1-np.sign(tmp_coord[:,0]))/float(2)) * box_dim[0]
						tmp_coord[:,1] -= (np.floor(2*tmp_coord[:,1]/float(box_dim[1])) + (1-np.sign(tmp_coord[:,1]))/float(2)) * box_dim[1]
						tmp_coord[:,2] -= (np.floor(2*tmp_coord[:,2]/float(box_dim[2])) + (1-np.sign(tmp_coord[:,2]))/float(2)) * box_dim[2]
						
						#keep those within the specified radius
						tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
			
						#update particle coverage (Chang algorithm for on the fly average and std dev - actually Knuth simple version as we add element one by one)
						tmp_pc_coord = np.shape(tmp_coord_within)[0]/float(np.shape(tmp_coord)[0])*100
						delta =  tmp_pc_coord - sizes_coverage["charges"]["avg"]["all sizes"][charge_g]
						sizes_coverage["charges"]["nb"]["all sizes"][charge_g] += 1
						sizes_coverage["charges"]["avg"]["all sizes"][charge_g] += delta / sizes_coverage["charges"]["nb"]["all sizes"][charge_g]
						sizes_coverage["charges"]["std"]["all sizes"][charge_g] += delta * (tmp_pc_coord - sizes_coverage["charges"]["avg"]["all sizes"][charge_g])

						#add number of particles within each slice
						bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
						bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)
						if len(bin_abs) > 0:
							tmp_bins_nb += np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0] * charges_groups[charge_g]["values"][q]
				density_charges_sizes["all sizes"][charge_g] += tmp_bins_nb
				density_charges_sizes["all sizes"]["total"] += tmp_bins_nb
	
	#case: proteins
	#==============
	else:
		#identify clusters
		#-----------------
		if args.m_algorithm != "density":
			clusters = detect_clusters_connectivity(get_distances(box_dim), box_dim)
		else:
			clusters = detect_clusters_density(get_distances(box_dim), box_dim)
	
		#process them
		#------------
		for cluster in clusters:		
			#create selection for current cluster and only process it if it's TM (find closest PO4 particles for each particles of clusters, if all are in the same leaflet then it's surfacic [NB: this is done at the CLUSTER level (the same criteria at the protein level would probably fail)])
			c_sele = MDAnalysis.core.AtomGroup.AtomGroup([])
			for p_index in cluster:
				c_sele += proteins_sele[p_index]
			tmp_c_sele_coordinates = c_sele.coordinates()
			dist_min_lower = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["lower"], box_dim), axis = 1)
			dist_min_upper = np.min(MDAnalysis.analysis.distances.distance_array(tmp_c_sele_coordinates, tmp_lip_coords["upper"], box_dim), axis = 1)
			dist = dist_min_upper - dist_min_lower
			if np.size(dist[dist>0]) != np.size(dist) and np.size(dist[dist>0]) !=0:
				
				c_size = np.size(cluster)
				#store new cluster size and add entry if necessary
				if c_size not in sizes_sampled:
					sizes_nb_clusters[c_size] = 0
					sizes_sampled.append(c_size)
					sizes_sampled = sorted(sizes_sampled)
					#particles
					sizes_coverage["particles"]["nb"][c_size] = {part: 0 for part in particles_def["labels"]}
					sizes_coverage["particles"]["avg"][c_size] = {part: 0 for part in particles_def["labels"]}
					sizes_coverage["particles"]["std"][c_size] = {part: 0 for part in particles_def["labels"]}
					density_particles_sizes_nb[c_size] = {part: np.zeros(2*bins_nb) for part in particles_def["labels"]}
					density_particles_sizes_pc[c_size] = {part: np.zeros(2*bins_nb) for part in particles_def["labels"]}
					#residues
					if args.residuesfilename != "no":
						sizes_coverage["residues"]["nb"][c_size] = {res: 0 for res in residues_def["labels"]}
						sizes_coverage["residues"]["avg"][c_size] = {res: 0 for res in residues_def["labels"]}
						sizes_coverage["residues"]["std"][c_size] = {res: 0 for res in residues_def["labels"]}
						density_residues_sizes_nb[c_size] = {res: np.zeros(2*bins_nb) for res in residues_def["labels"]}
						density_residues_sizes_pc[c_size] = {res: np.zeros(2*bins_nb) for res in residues_def["labels"]}
					#charges
					if args.chargesfilename != "no":
						sizes_coverage["charges"]["nb"][c_size] = {q: 0 for q in charges_groups.keys() + ["total"]}
						sizes_coverage["charges"]["avg"][c_size] = {q: 0 for q in charges_groups.keys() + ["total"]}
						sizes_coverage["charges"]["std"][c_size] = {q: 0 for q in charges_groups.keys() + ["total"]}
						density_charges_sizes[c_size] = {q: np.zeros(2*bins_nb) for q in charges_groups.keys() + ["total"]}
	
				if args.cluster_groups_file != "no":
					g_index = groups_sizes_dict[c_size]
					if g_index not in groups_sampled:
						groups_sampled.append(g_index)
						groups_sampled = sorted(groups_sampled)
				
				#add to number of cluster processed for this size
				sizes_nb_clusters[c_size] += 1
				sizes_nb_clusters["all sizes"] += 1
				if args.cluster_groups_file != "no":
					groups_nb_clusters[g_index] += 1
	
				#get z coord of cluster center of geometry
				cluster_cog = calculate_cog(tmp_c_sele_coordinates, box_dim)
										
				#density profile: particles
				#--------------------------
				for part in particles_def["labels"]:
					if part == "peptide":
						tmp_part_sele = c_sele
					else:
						tmp_part_sele = particles_def["sele"][part]
					if tmp_part_sele.numberOfAtoms() > 0:
						
						#center cluster (x,y) coordinates on its cog (x,y) coordinates
						tmp_coord = tmp_part_sele.coordinates()
						tmp_coord[:,0] -= cluster_cog[0]
						tmp_coord[:,1] -= cluster_cog[1]
										
						#center cluster z coordinates on the bilayer center z coordinate
						tmp_coord[:,2] -= z_middle_instant
					
						#deal with pbc and center axis on 0
						# -> we test if we are further away from 0 than half the box dim, if so we translate by a box dim.
						# -> the only trick is that for negative coord floor returns a number offset by 1 compared to what we'd like (-1 = within range, -2 - outside), so we need to add +1 to this for negative coord (hence the bit with sign)
						tmp_coord[:,0] -= (np.floor(2*tmp_coord[:,0]/float(box_dim[0])) + (1-np.sign(tmp_coord[:,0]))/float(2)) * box_dim[0]
						tmp_coord[:,1] -= (np.floor(2*tmp_coord[:,1]/float(box_dim[1])) + (1-np.sign(tmp_coord[:,1]))/float(2)) * box_dim[1]
						tmp_coord[:,2] -= (np.floor(2*tmp_coord[:,2]/float(box_dim[2])) + (1-np.sign(tmp_coord[:,2]))/float(2)) * box_dim[2]

						#keep those within the specified radius
						tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
	
						#update particle coverage (Chang algorithm for on the fly average and std dev - actually Knuth simple version as we add element one by one)
						tmp_pc_coord = np.shape(tmp_coord_within)[0]/float(np.shape(tmp_coord)[0])*100
						for c in [c_size, "all sizes"]:
							delta =  tmp_pc_coord - sizes_coverage["particles"]["avg"][c][part]
							sizes_coverage["particles"]["nb"][c][part] += 1
							sizes_coverage["particles"]["avg"][c][part] += delta / sizes_coverage["particles"]["nb"][c][part]
							sizes_coverage["particles"]["std"][c][part] += delta * (tmp_pc_coord - sizes_coverage["particles"]["avg"][c][part])
						if args.cluster_groups_file != "no":
							delta =  tmp_pc_coord - groups_coverage["particles"]["avg"][g_index][part]
							groups_coverage["particles"]["nb"][g_index][part] += 1
							groups_coverage["particles"]["avg"][g_index][part] += delta / groups_coverage["particles"]["nb"][g_index][part]
							groups_coverage["particles"]["std"][g_index][part] += delta * (tmp_pc_coord - groups_coverage["particles"]["avg"][g_index][part])
						
						#add number of particles within each slice					
						tmp_bins_nb = np.zeros(2*bins_nb)
						bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
						bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)
						if len(bin_abs) > 0:				
							tmp_bins_nb = np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0]
						density_particles_sizes_nb[c_size][part] += tmp_bins_nb
						density_particles_sizes_nb["all sizes"][part] += tmp_bins_nb
						if args.cluster_groups_file != "no":
							density_particles_groups_nb[g_index][part] += tmp_bins_nb
	
				#density profile: residues
				#-------------------------
				if args.residuesfilename != "no":
					for res in residues_def["labels"]:
						#select particles of given residues
						tmp_res_sele = c_sele.selectAtoms(residues_def["sele_string"][res])
						if tmp_res_sele.numberOfAtoms() > 0:
							
							#center cluster (x,y) coordinates on its cog (x,y) coordinates
							tmp_coord = tmp_res_sele.coordinates()
							tmp_coord[:,0] -= cluster_cog[0]
							tmp_coord[:,1] -= cluster_cog[1]
											
							#center cluster z coordinates on the bilayer center z coordinate
							tmp_coord[:,2] -= z_middle_instant
						
							#deal with pbc and center axis on 0
							tmp_coord[:,0] -= (np.floor(2*tmp_coord[:,0]/float(box_dim[0])) + (1-np.sign(tmp_coord[:,0]))/float(2)) * box_dim[0]
							tmp_coord[:,1] -= (np.floor(2*tmp_coord[:,1]/float(box_dim[1])) + (1-np.sign(tmp_coord[:,1]))/float(2)) * box_dim[1]
							tmp_coord[:,2] -= (np.floor(2*tmp_coord[:,2]/float(box_dim[2])) + (1-np.sign(tmp_coord[:,2]))/float(2)) * box_dim[2]
												
							#keep those within the specified radius
							tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
		
							#update particle coverage (Chang algorithm for on the fly average and std dev - actually Knuth simple version as we add element one by one)
							tmp_pc_coord = np.shape(tmp_coord_within)[0]/float(np.shape(tmp_coord)[0])*100
							for c in [c_size, "all sizes"]:
								delta =  tmp_pc_coord - sizes_coverage["residues"]["avg"][c][res]
								sizes_coverage["residues"]["nb"][c][res] += 1
								sizes_coverage["residues"]["avg"][c][res] += delta / sizes_coverage["residues"]["nb"][c][res]
								sizes_coverage["residues"]["std"][c][res] += delta * (tmp_pc_coord - sizes_coverage["residues"]["avg"][c][res])
							if args.cluster_groups_file != "no":
								delta =  tmp_pc_coord - groups_coverage["residues"]["avg"][g_index][res]
								groups_coverage["residues"]["nb"][g_index][res] += 1
								groups_coverage["residues"]["avg"][g_index][res] += delta / groups_coverage["residues"]["nb"][g_index][res]
								groups_coverage["residues"]["std"][g_index][res] += delta * (tmp_pc_coord - groups_coverage["residues"]["avg"][g_index][res])
							
							#add number of particles within each slice					
							tmp_bins_nb = np.zeros(2*bins_nb)
							bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
							bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)
							if len(bin_abs) > 0:
								tmp_bins_nb = np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0]
							density_residues_sizes_nb[c_size][res] += tmp_bins_nb
							density_residues_sizes_nb["all sizes"][res] += tmp_bins_nb
							if args.cluster_groups_file != "no":
								density_residues_groups_nb[g_index][res] += tmp_bins_nb
				
				#density profile: charges
				#------------------------
				if args.chargesfilename != "no":
					for charge_g in charges_groups.keys():
						tmp_bins_nb = np.zeros(2*bins_nb)
						for q in charges_groups[charge_g]["names"]:
							if charge_g == "peptide":
								tmp_q_sele = c_sele.selectAtoms(charges_groups[charge_g]["sele_string"][q])
							else:
								tmp_q_sele = charges_groups[charge_g]["sele"][q]
							if tmp_q_sele.numberOfAtoms() > 0:
													
								#center cluster coordinates on the (x,y) coordinates of its cog
								tmp_coord = tmp_q_sele.coordinates()
								tmp_coord[:,0] -= cluster_cog[0]
								tmp_coord[:,1] -= cluster_cog[1]
							
								#center cluster z coordinates on the bilayer center z coordinate
								tmp_coord[:,2] -= z_middle_instant
		
								#deal with pbc and center axis on 0
								tmp_coord[:,0] -= (np.floor(2*tmp_coord[:,0]/float(box_dim[0])) + (1-np.sign(tmp_coord[:,0]))/float(2)) * box_dim[0]
								tmp_coord[:,1] -= (np.floor(2*tmp_coord[:,1]/float(box_dim[1])) + (1-np.sign(tmp_coord[:,1]))/float(2)) * box_dim[1]
								tmp_coord[:,2] -= (np.floor(2*tmp_coord[:,2]/float(box_dim[2])) + (1-np.sign(tmp_coord[:,2]))/float(2)) * box_dim[2]
								
								#keep those within the specified radius
								tmp_coord_within = tmp_coord[tmp_coord[:,0]**2 + tmp_coord[:,1]**2 < args.slices_radius**2]
					
								#update particle coverage (Chang algorithm for on the fly average and std dev - actually Knuth simple version as we add element one by one)
								tmp_pc_coord = np.shape(tmp_coord_within)[0]/float(np.shape(tmp_coord)[0])*100
								for c in [c_size, "all sizes"]:
									delta =  tmp_pc_coord - sizes_coverage["charges"]["avg"][c][charge_g]
									sizes_coverage["charges"]["nb"][c][charge_g] += 1
									sizes_coverage["charges"]["avg"][c][charge_g] += delta / sizes_coverage["charges"]["nb"][c][charge_g]
									sizes_coverage["charges"]["std"][c][charge_g] += delta * (tmp_pc_coord - sizes_coverage["charges"]["avg"][c][charge_g])
								if args.cluster_groups_file != "no":
									delta =  tmp_pc_coord - groups_coverage["charges"]["avg"][g_index][charge_g]
									groups_coverage["charges"]["nb"][g_index][charge_g] += 1
									groups_coverage["charges"]["avg"][g_index][charge_g] += delta / groups_coverage["charges"]["nb"][g_index][charge_g]
									groups_coverage["charges"]["std"][g_index][charge_g] += delta * (tmp_pc_coord - groups_coverage["charges"]["avg"][g_index][charge_g])
		
								#add number of particles within each slice
								bin_rel = np.floor(tmp_coord_within[:,2]/float(args.slices_thick)).astype(int)
								bin_abs = bin_rel[abs(bin_rel) < bins_nb_max] + int(bins_nb)
								if len(bin_abs) > 0:
									tmp_bins_nb += np.histogram(bin_abs, np.arange(2*bins_nb + 1))[0] * charges_groups[charge_g]["values"][q]
						density_charges_sizes[c_size][charge_g] += tmp_bins_nb 
						density_charges_sizes[c_size]["total"] += tmp_bins_nb
						density_charges_sizes["all sizes"][charge_g] += tmp_bins_nb
						density_charges_sizes["all sizes"]["total"] += tmp_bins_nb
						if args.cluster_groups_file != "no":
							density_charges_groups[g_index][charge_g] += tmp_bins_nb
							density_charges_groups[g_index]["total"] += tmp_bins_nb
	
	return
def calculate_stats():													#DONE
	
	global z_upper
	global z_lower
	global max_density_particles_pc
	max_density_particles_pc  = float("-inf")
	if args.residuesfilename != "no":
		global max_density_residues_pc
		max_density_residues_pc  = float("-inf")
	if args.chargesfilename != "no":
		global max_density_charges
		global min_density_charges
		max_density_charges = float("-inf")
		min_density_charges = float("inf")
	
	#coords
	#======
	z_upper /= float(nb_frames_to_process)
	z_lower /= float(nb_frames_to_process)
	
	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:
		
		#calculate normalisation factor for each group for each size
		tmp_normalisation = {}
		for part_g in particles_groups.keys():
			tmp_normalisation[part_g] = 0
			for part in particles_groups[part_g]:
				tmp_normalisation[part_g] += np.sum(density_particles_sizes_nb[c_size][part])
		
		#density profile: particles
		#--------------------------
		for part in particles_def["labels"]:
			#% of particles covered
			sizes_coverage["particles"]["avg"][c_size][part] = round(sizes_coverage["particles"]["avg"][c_size][part],1)
			if sizes_coverage["particles"]["nb"][c_size][part] > 1:
				sizes_coverage["particles"]["std"][c_size][part] = round(math.sqrt(sizes_coverage["particles"]["std"][c_size][part] / float(sizes_coverage["particles"]["nb"][c_size][part])),1)
			else:
				sizes_coverage["particles"]["std"][c_size][part] = "nan"
		
			#relative density
			if tmp_normalisation[particles_def["group"][part]] > 0:
				density_particles_sizes_pc[c_size][part] = density_particles_sizes_nb[c_size][part] / float(tmp_normalisation[particles_def["group"][part]])
			
			#update scale
			if part != "Na+" and part != "Cl-":
				max_density_particles_pc = max(max_density_particles_pc, max(density_particles_sizes_pc[c_size][part]))
				if (part == "peptide" or part == "water") and args.residuesfilename != "no":
					max_density_residues_pc = max(max_density_residues_pc, max(density_particles_sizes_pc[c_size][part]))
		
		#density profile: residues
		#-------------------------
		if args.residuesfilename != "no":
			for res in residues_def["labels"]:
				#% of particles covered
				sizes_coverage["residues"]["avg"][c_size][res] = round(sizes_coverage["residues"]["avg"][c_size][res],1)
				if sizes_coverage["residues"]["nb"][c_size][res] > 1:
					sizes_coverage["residues"]["std"][c_size][res] = round(math.sqrt(sizes_coverage["residues"]["std"][c_size][res] / float(sizes_coverage["residues"]["nb"][c_size][res])),1)
				else:
					sizes_coverage["residues"]["std"][c_size][res] = "nan"
				
				#relative density
				density_residues_sizes_pc[c_size][res] = density_residues_sizes_nb[c_size][res] / float(tmp_normalisation["peptide"])

		#density profile: charges
		#------------------------
		if args.chargesfilename != "no":
			tmp_normalisation = sizes_nb_clusters[c_size] * slice_volume
			for charge_g in charges_groups.keys() + ["total"]:
				#% of particles covered
				if charge_g != "total":
					sizes_coverage["charges"]["avg"][c_size][charge_g] = round(sizes_coverage["charges"]["avg"][c_size][charge_g],1)
					if sizes_coverage["charges"]["nb"][c_size][charge_g] > 1:
						sizes_coverage["charges"]["std"][c_size][charge_g] = round(math.sqrt(sizes_coverage["charges"]["std"][c_size][charge_g] / float(sizes_coverage["charges"]["nb"][c_size][charge_g])),1)
					else:
						sizes_coverage["charges"]["std"][c_size][charge_g] = "nan"
				
				#absolute density
				density_charges_sizes[c_size][charge_g] /= float(tmp_normalisation)
				
				#update scale
				max_density_charges = max(max_density_charges, max(density_charges_sizes[c_size][charge_g]))
				min_density_charges = min(min_density_charges, min(density_charges_sizes[c_size][charge_g]))
	
	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:
			
			tmp_normalisation = groups_nb_clusters[g_index] * slice_volume
			
			#density profile: particles
			#--------------------------
			for part in particles_def["labels"]:
				#% of particles covered
				groups_coverage["particles"]["avg"][g_index][part] = round(groups_coverage["particles"]["avg"][g_index][part],1)
				if groups_coverage["particles"]["nb"][g_index][part] > 1:
					groups_coverage["particles"]["std"][g_index][part] = round(math.sqrt(groups_coverage["particles"]["std"][g_index][part] / float(groups_coverage["particles"]["nb"][g_index][part])),1)
				else:
					groups_coverage["particles"]["std"][g_index][part] = "nan"
			
				#asbolute density
				density_particles_groups_nb[g_index][part] /= float(tmp_normalisation)

				#relative density
				if np.sum(density_particles_groups_nb[g_index][part]) > 0:
					density_particles_groups_pc[g_index][part] = density_particles_groups_nb[g_index][part] / float(np.sum(density_particles_groups_nb[g_index][part]))
				
				#update scale
				max_density_particles_pc = max(max_density_particles_pc, max(density_particles_groups_pc[g_index][part]))
				if (part == "peptide" or part == "water") and args.residuesfilename != "no":
					max_density_residues_pc = max(max_density_residues_pc, max(density_particles_groups_pc[g_index][part]))
			
			#density profile: residues
			#-------------------------
			if args.residuesfilename != "no":
				for res in residues_def["labels"]:
					#% of particles covered
					groups_coverage["residues"]["avg"][g_index][res] = round(groups_coverage["residues"]["avg"][g_index][res],1)
					if groups_coverage["residues"]["nb"][g_index][res] > 1:
						groups_coverage["residues"]["std"][g_index][res] = round(math.sqrt(groups_coverage["residues"]["std"][g_index][res] / float(groups_coverage["residues"]["nb"][g_index][res])),1)
					else:
						groups_coverage["residues"]["std"][g_index][res] = "nan"
					
					#absolute density
					density_residues_groups_nb[g_index][res] /= float(tmp_normalisation)

					#relative density
					density_residues_groups_pc[g_index][res] = density_residues_groups_nb[g_index][res] / float(np.sum(density_particles_groups_nb[g_index]["peptide"]))
		
			#density profile: charges
			#------------------------
			if args.chargesfilename != "no":
				for charge_g in charges_groups.keys() + ["total"]:
					#% of particles covered
					if charge_g != "total":
						groups_coverage["charges"]["avg"][g_index][charge_g] = round(groups_coverage["charges"]["avg"][g_index][charge_g],1)
						if groups_coverage["charges"]["nb"][g_index][charge_g] > 1:
							groups_coverage["charges"]["std"][g_index][charge_g] = round(math.sqrt(groups_coverage["charges"]["std"][g_index][charge_g] / float(groups_coverage["charges"]["nb"][g_index][charge_g])),1)
						else:
							groups_coverage["charges"]["std"][g_index][charge_g] = "nan"
					
					#absolute density
					density_charges_groups[g_index][charge_g] /= float(tmp_normalisation)
					
					#update scale
					max_density_charges = max(max_density_charges, max(density_charges_groups[g_index][charge_g]))
					min_density_charges = min(min_density_charges, min(density_charges_groups[g_index][charge_g]))
	
	return

#=========================================================================================
# outputs: TO DO
#=========================================================================================

def density_write_particles():											#DONE

	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:		
		#open files
		if c_size == "all sizes":
			tmp_file = '1_1_density_profile_particles_all_sizes'
		else:
			tmp_file = '1_1_density_profile_particles_' + str(c_size)
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_particles/xvg/' + str(tmp_file) + '.xvg'
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_particles/xvg/' + str(tmp_file) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_xvg = open(filename_xvg, 'w')
		
		#general header
		output_txt.write("@[relative particles frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_file) + ".xvg.\n")
		output_xvg.write("# [relative particles frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
		output_xvg.write("# cluster detection method:\n")
		if protein_pres:
			output_xvg.write("#  -> nb of proteins: " + str(proteins_nb) + "\n")
			if args.m_algorithm == "min":
				output_xvg.write("#  -> connectivity based (min distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm == "cog":
				output_xvg.write("#  -> connectivity based (cog distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_xvg.write("#  -> density based (DBSCAN)\n")
				output_xvg.write("#  -> search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		output_xvg.write("# volume properties:\n")
		output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
		output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
		output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
		output_xvg.write("# selection statistics (% of total particles within volume [if protein: % of cluster]):\n")
		for part in particles_def["labels"]:
			output_xvg.write("#  -> " + str(part) + ": " + str(sizes_coverage["particles"]["avg"][c_size][part]) + "% (" + str(sizes_coverage["particles"]["std"][c_size][part]) + ")\n")
		if protein_pres:
			output_xvg.write("# nb of clusters which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(sizes_nb_clusters[c_size]) + "\n")
		else:
			output_xvg.write("# nb of frames which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
		
		#xvg metadata
		if c_size == "all sizes":
			output_xvg.write("@ title \"Relative particles frequency along z around TM clusters\"\n")
		else:
			output_xvg.write("@ title \"Relative particles frequency along z around TM clusters of size " + str(c_size) + "\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"relative particles frequency (Angstrom-3)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(len(particles_def["labels"])) + "\n")
		for part_index in range(0,len(particles_def["labels"])):
			part = particles_def["labels"][part_index]
			output_xvg.write("@ s" + str(part_index) + " legend \"" + str(part) + "\"\n")
			output_txt.write(str(tmp_file) + "," + str(part_index+1) + "," + str(part) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"][part])) + "\n")
		output_txt.close()
		
		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n])
			for part_index in range(0,len(particles_def["labels"])):
				part = particles_def["labels"][part_index]
				results += "	" + "{:.6e}".format(density_particles_sizes_pc[c_size][part][n])
			output_xvg.write(results + "\n")	
		output_xvg.close()
				
	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:
			#open files
			tmp_file = '2_1_density_profile_particles_' + str(groups_labels[g_index])
			filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/1_particles/xvg/' + str(tmp_file) + '.xvg'
			filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/1_particles/xvg/' + str(tmp_file) + '.txt'
			output_txt = open(filename_txt, 'w')
			output_xvg = open(filename_xvg, 'w')
			
			#general header
			output_txt.write("@[relative particles frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
			output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_file) + ".xvg.\n")
			output_xvg.write("# [relative particles frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
			output_xvg.write("# cluster detection method:\n")
			output_xvg.write("#  -> nb of proteins: " + str(proteins_nb) + "\n")
			if args.m_algorithm == "min":
				output_xvg.write("#  -> connectivity based (min distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm == "cog":
				output_xvg.write("#  -> connectivity based (cog distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_xvg.write("#  -> density based (DBSCAN)\n")
				output_xvg.write("#  -> search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
			output_xvg.write("# volume properties:\n")
			output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
			output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
			output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
			output_xvg.write("# selection statistics (% of total particles within volume [if protein: % of cluster]):\n")
			for part in particles_def["labels"]:
				output_xvg.write("#  -> " + str(part) + ": " + str(groups_coverage["particles"]["avg"][g_index][part]) + "% (" + str(groups_coverage["particles"]["std"][g_index][part]) + ")\n")
			output_xvg.write("# nb of clusters which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(groups_nb_clusters[g_index]) + "\n")
			
			#xvg metadata
			output_xvg.write("@ title \"Relative particles frequency along z around TM clusters of sizes " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"relative particles frequency (Angstrom-3)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			output_xvg.write("@ legend length " + str(len(particles_def["labels"])) + "\n")
			for part_index in range(0,len(particles_def["labels"])):
				part = particles_def["labels"][part_index]
				output_xvg.write("@ s" + str(part_index) + " legend \"" + str(part) + "\"\n")
				output_txt.write(str(tmp_file) + "," + str(part_index+1) + "," + str(part) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"][part])) + "\n")
			output_txt.close()
			
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n])
				for part_index in range(0,len(particles_def["labels"])):
					part = particles_def["labels"][part_index]
					results += "	" + "{:.6e}".format(density_particles_groups_pc[g_index][part][n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	return
def density_write_residues():											#DONE

	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:		
		
		#open files
		if c_size == "all sizes":
			tmp_file = '1_3_density_profile_residues_all_sizes'
		else:
			tmp_file = '1_3_density_profile_residues_' + str(c_size)
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/3_residues/xvg/' + str(tmp_file) + '.xvg'
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/3_residues/xvg/' + str(tmp_file) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_xvg = open(filename_xvg, 'w')
		
		#general header
		output_txt.write("@[relative residues frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_file) + ".xvg.\n")
		output_xvg.write("# [relative residues frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
		output_xvg.write("# cluster detection method:\n")
		output_xvg.write("#  -> nb of proteins: " + str(proteins_nb) + "\n")
		if args.m_algorithm == "min":
			output_xvg.write("#  -> connectivity based (min distances)\n")
			output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		elif args.m_algorithm == "cog":
			output_xvg.write("#  -> connectivity based (cog distances)\n")
			output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
		else:
			output_xvg.write("#  -> density based (DBSCAN)\n")
			output_xvg.write("#  -> search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		output_xvg.write("# volume properties:\n")
		output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
		output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
		output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
		output_xvg.write("# selection statistics (% of total particles within volume [if protein: % of cluster]):\n")
		for res in residues_def["labels"]:
			output_xvg.write("#  -> " + str(res) + ": " + str(sizes_coverage["residues"]["avg"][c_size][res]) + "% (" + str(sizes_coverage["residues"]["std"][c_size][res]) + ")\n")
		if protein_pres:
			output_xvg.write("# nb of clusters which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(sizes_nb_clusters[c_size]) + "\n")
		else:
			output_xvg.write("# nb of frames which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
		
		#xvg metadata
		if c_size == "all sizes":
			output_xvg.write("@ title \"Relative residues frequency along z around TM clusters\"\n")
		else:
			output_xvg.write("@ title \"Relative residues frequency along z around TM clusters of size " + str(c_size) + "\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"relative residues frequency (Angstrom-3)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		if water_pres:
			output_xvg.write("@ legend length " + str(len(residues_def["labels"]) + 2) + "\n")
		else:
			output_xvg.write("@ legend length " + str(len(residues_def["labels"]) + 1) + "\n")
		for res_index in range(0,len(residues_def["labels"])):
			res = residues_def["labels"][res_index]
			output_xvg.write("@ s" + str(res_index) + " legend \"" + str(res) + "\"\n")
			output_txt.write(str(tmp_file) + "," + str(res_index+1) + "," + str(res) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(residues_def["colour"][res])) + "\n")
		output_xvg.write("@ s" + str(len(residues_def["labels"])) + " legend \"peptide\"\n")
		output_txt.write(str(tmp_file) + "," + str(len(residues_def["labels"]) + 1) + ",peptide," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"]["peptide"])) + "\n")
		if water_pres:
			output_xvg.write("@ s" + str(len(residues_def["labels"]) + 1) + " legend \"water\"\n")
			output_txt.write(str(tmp_file) + "," + str(len(residues_def["labels"]) + 2) + ",water," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"]["water"])) + "\n")		
		output_txt.close()
		
		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n])
			for res_index in range(0,len(residues_def["labels"])):
				res = residues_def["labels"][res_index]
				results += "	" + "{:.6e}".format(density_residues_sizes_pc[c_size][res][n])
			results += "	" + "{:.6e}".format(density_particles_sizes_pc[c_size]["peptide"][n])
			if water_pres:
				results += "	" + "{:.6e}".format(density_particles_sizes_pc[c_size]["water"][n])
			output_xvg.write(results + "\n")	
		output_xvg.close()

	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:

			#open files
			tmp_file = '2_3_density_profile_residues_' + str(groups_labels[g_index])
			filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/3_residues/xvg/' + str(tmp_file) + '.xvg'
			filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/3_residues/xvg/' + str(tmp_file) + '.txt'
			output_txt = open(filename_txt, 'w')
			output_xvg = open(filename_xvg, 'w')
			
			#general header
			output_txt.write("@[relative residues frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
			output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in " + str(tmp_file) + ".xvg.\n")
			output_xvg.write("# [relative residues frequency profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
			output_xvg.write("# cluster detection method:\n")
			output_xvg.write("#  -> nb of proteins: " + str(proteins_nb) + "\n")
			if args.m_algorithm == "min":
				output_xvg.write("#  -> connectivity based (min distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm == "cog":
				output_xvg.write("#  -> connectivity based (cog distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_xvg.write("#  -> density based (DBSCAN)\n")
				output_xvg.write("#  -> search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
			output_xvg.write("# volume properties:\n")
			output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
			output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
			output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
			output_xvg.write("# selection statistics (% of total particles within volume [if protein: % of cluster]):\n")
			for res in residues_def["labels"]:
				output_xvg.write("#  -> " + str(res) + ": " + str(groups_coverage["residues"]["avg"][g_index][res]) + "% (" + str(groups_coverage["residues"]["std"][g_index][res]) + ")\n")
			output_xvg.write("# nb of clusters which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(groups_nb_clusters[g_index]) + "\n")
			
			#xvg metadata
			output_xvg.write("@ title \"Relative residues frequency along z around TM clusters of sizes " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"relative residues frequency (Angstrom-3)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			if water_pres:
				output_xvg.write("@ legend length " + str(len(residues_def["labels"]) + 2) + "\n")
			else:
				output_xvg.write("@ legend length " + str(len(residues_def["labels"]) + 1) + "\n")
			for res_index in range(0,len(residues_def["labels"])):
				res = residues_def["labels"][res_index]
				output_xvg.write("@ s" + str(res_index) + " legend \"" + str(res) + "\"\n")
				output_txt.write(str(tmp_file) + "," + str(res_index+1) + "," + str(res) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(residues_def["colour"][res])) + "\n")
			output_xvg.write("@ s" + str(len(residues_def["labels"])) + " legend \"peptide\"\n")
			output_txt.write(str(tmp_file) + "," + str(len(residues_def["labels"]) + 1) + ",peptide," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"]["peptide"])) + "\n")
			if water_pres:
				output_xvg.write("@ s" + str(len(residues_def["labels"]) + 1) + " legend \"water\"\n")
				output_txt.write(str(tmp_file) + "," + str(len(residues_def["labels"]) + 2) + ",water," + mcolors.rgb2hex(mcolorconv.to_rgb(particles_def["colour"]["water"])) + "\n")		
			output_txt.close()
			
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n])
				for res_index in range(0,len(residues_def["labels"])):
					res = residues_def["labels"][res_index]
					results += "	" + "{:.6e}".format(density_residues_groups_pc[g_index][res][n])
				results += "	" + "{:.6e}".format(density_particles_groups_pc[g_index]["peptide"][n])
				if water_pres:
					results += "	" + "{:.6e}".format(density_particles_groups_pc[g_index]["water"][n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	return
def density_write_charges():											#DONE
	
	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:

		#open files
		if c_size == "all sizes":
			tmp_file = '1_2_density_profile_charges_all_sizes'
		else:
			tmp_file = '1_2_density_profile_charges_' + str(c_size)
		filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/2_charges/xvg/' + str(tmp_file) + '.xvg'
		filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/2_charges/xvg/' + str(tmp_file) + '.txt'
		output_txt = open(filename_txt, 'w')
		output_xvg = open(filename_xvg, 'w')
		
		#general header
		output_txt.write("@[charge density profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
		output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_thickness_species.xvg.\n")
		output_xvg.write("# [charge density profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
		output_xvg.write("# cluster detection method:\n")
		if protein_pres:
			output_xvg.write("#  -> nb of proteins: " + str(proteins_nb) + "\n")
			if args.m_algorithm == "min":
				output_xvg.write("#  -> connectivity based (min distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm == "cog":
				output_xvg.write("#  -> connectivity based (cog distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_xvg.write("#  -> density based (DBSCAN)\n")
				output_xvg.write("#  -> search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
		output_xvg.write("# volume properties:\n")
		output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
		output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
		output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
		output_xvg.write("# selection statistics (% of total particles within volume):\n")
		for charge_g in charges_groups.keys():
			output_xvg.write("#  -> " + str(charge_g) + ": " + str(sizes_coverage["charges"]["avg"][c_size][charge_g]) + "% (" + str(sizes_coverage["charges"]["std"][c_size][charge_g]) + ")\n")
		if protein_pres:
			output_xvg.write("# nb of clusters which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(sizes_nb_clusters[c_size]) + "\n")
		else:
			output_xvg.write("# nb of frames which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
		
		#xvg metadata
		if c_size == "all sizes":
			output_xvg.write("@ title \"Charge density profile along z around TM clusters\"\n")
		else:
			output_xvg.write("@ title \"Charge density profile along z around TM clusters of size " + str(c_size) + "\"\n")
		output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
		output_xvg.write("@ yaxis label \"charge density (e.Angstrom-3)\"\n")
		output_xvg.write("@ autoscale ONREAD xaxes\n")
		output_xvg.write("@ TYPE XY\n")
		output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
		output_xvg.write("@ legend on\n")
		output_xvg.write("@ legend box on\n")
		output_xvg.write("@ legend loctype view\n")
		output_xvg.write("@ legend 0.98, 0.8\n")
		output_xvg.write("@ legend length " + str(len(charges_groups.keys()) + 1) + "\n")
		charges_groups_keys = charges_groups.keys()
		for q_index in range(0,len(charges_groups_keys)):
			charge_group = charges_groups_keys[q_index]
			output_xvg.write("@ s" + str(q_index) + " legend \"" + str(charge_group) + "\"\n")
			output_txt.write(str(tmp_file) + "," + str(q_index) + "," + str(charge_group) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(charges_colours[charge_group])) + "\n")
		output_xvg.write("@ s" + str(len(charges_groups.keys())) + " legend \"total\"\n")
		output_txt.write(str(tmp_file) + "," + str(len(charges_groups.keys())+1) + ",total," + mcolors.rgb2hex(mcolorconv.to_rgb(charges_colours["total"])) + "\n")
		output_txt.close()
		
		#data
		for n in range(0,2*bins_nb):
			results = str(bins_labels[n])
			for q_index in range(0,len(charges_groups_keys)):
				charge_g = charges_groups_keys[q_index]
				results += "	" + "{:.6e}".format(density_charges_sizes[c_size][charge_g][n])
			results += "	" + "{:.6e}".format(density_charges_sizes[c_size]["total"][n])
			output_xvg.write(results + "\n")	
		output_xvg.close()
		
	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:
			#open files
			tmp_file = '2_2_density_profile_charges_' + str(groups_labels[g_index])
			filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_charges/xvg/' + str(tmp_file) + '.xvg'
			filename_txt = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_charges/xvg/' + str(tmp_file) + '.txt'
			output_txt = open(filename_txt, 'w')
			output_xvg = open(filename_xvg, 'w')
			
			#general header
			output_txt.write("@[charge density profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
			output_txt.write("@Use this file as the argument of the -c option of the script 'xvg_animate' in order to make a time lapse movie of the data in 1_2_thickness_species.xvg.\n")
			output_xvg.write("# [charge density profile - written by cluster_density_profile v" + str(version_nb) + "]\n")
			output_xvg.write("# cluster detection method:\n")
			output_xvg.write("#  -> nb of proteins: " + str(proteins_nb) + "\n")
			if args.m_algorithm == "min":
				output_xvg.write("#  -> connectivity based (min distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			elif args.m_algorithm == "cog":
				output_xvg.write("#  -> connectivity based (cog distances)\n")
				output_xvg.write("#  -> contact cutoff = " + str(args.cutoff_connect) + " Angstrom\n")
			else:
				output_xvg.write("#  -> density based (DBSCAN)\n")
				output_xvg.write("#  -> search radius = " + str(args.dbscan_dist) + " Angstrom, nb of neighbours = " + str(args.dbscan_nb) + "\n")
			output_xvg.write("# volume properties:\n")
			output_xvg.write("#  -> cylinder radius: " + str(args.slices_radius) + " (Angstrom)\n")
			output_xvg.write("#  -> slices thickness: " + str(args.slices_thick) + " (Angstrom)\n")
			output_xvg.write("#  -> slices volume: " + str(round(slice_volume,2)) + " (Angstrom3)\n")
			output_xvg.write("# selection statistics (% of total particles within volume):\n")
			for charge_g in charges_groups.keys():
				output_xvg.write("#  -> " + str(charge_g) + ": " + str(groups_coverage["charges"]["avg"][g_index][charge_g]) + "% (" + str(groups_coverage["charges"]["std"][g_index][charge_g]) + ")\n")
			output_xvg.write("# nb of clusters which contributed to this profile:\n")
			output_xvg.write("# -> weight = " + str(groups_nb_clusters[g_index]) + "\n")
			
			#xvg metadata
			if c_size == "all sizes":
				output_xvg.write("@ title \"Charge density profile along z around TM clusters\"\n")
			else:
				output_xvg.write("@ title \"Charge density profile along z around TM clusters of sizes " + str(groups_labels[g_index]) + "\"\n")
			output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
			output_xvg.write("@ yaxis label \"charge density (e.Angstrom-3)\"\n")
			output_xvg.write("@ autoscale ONREAD xaxes\n")
			output_xvg.write("@ TYPE XY\n")
			output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
			output_xvg.write("@ legend on\n")
			output_xvg.write("@ legend box on\n")
			output_xvg.write("@ legend loctype view\n")
			output_xvg.write("@ legend 0.98, 0.8\n")
			output_xvg.write("@ legend length " + str(len(charges_groups.keys()) + 1) + "\n")
			charges_groups_keys = charges_groups.keys()
			for q_index in range(0,len(charges_groups_keys)):
				charge_group = charges_groups_keys[q_index]
				output_xvg.write("@ s" + str(q_index) + " legend \"" + str(charge_group) + "\"\n")
				output_txt.write(str(tmp_file) + "," + str(q_index) + "," + str(charge_group) + "," + mcolors.rgb2hex(mcolorconv.to_rgb(charges_colours[charge_group])) + "\n")
			output_xvg.write("@ s" + str(len(charges_groups.keys())) + " legend \"total\"\n")
			output_txt.write(str(tmp_file) + "," + str(len(charges_groups.keys())+1) + ",total," + mcolors.rgb2hex(mcolorconv.to_rgb(charges_colours["total"])) + "\n")
			output_txt.close()
			
			#data
			for n in range(0,2*bins_nb):
				results = str(bins_labels[n])
				for q_index in range(0,len(charges_groups_keys)):
					charge_g = charges_groups_keys[q_index]
					results += "	" + "{:.6e}".format(density_charges_groups[g_index][charge_g][n])
				results += "	" + "{:.6e}".format(density_charges_groups[g_index]["total"][n])
				output_xvg.write(results + "\n")	
			output_xvg.close()

	return

def density_graph_particles():											#DONE
			
	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:
		#filenames
		if c_size == "all sizes":
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_particles/png/1_1_density_profile_particles_all_sizes.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_particles/1_1_density_profile_particles_all_sizes.svg'
		else:
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_particles/png/1_1_density_profile_particles_' + str(c_size) +'.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/1_particles/1_1_density_profile_particles_' + str(c_size) +'.svg'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		if c_size == "all sizes":
			fig.suptitle("Relative particles frequency along z around TM clusters")
		else:
			fig.suptitle("Relative particles frequency along z around TM clusters of size " + str(c_size))

		#plot data
		ax = fig.add_subplot(111)
		for part in particles_def["labels"]:
			plt.plot(range(0,2*bins_nb), density_particles_sizes_pc[c_size][part], color = particles_def["colour"][part], label = str(part))
		plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, 0, max_density_particles_pc, linestyles = 'dashed')
		plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, 0, max_density_particles_pc, linestyles = 'dashed')
		plt.vlines(bins_nb, 0, max_density_particles_pc, linestyles = 'dashdot')
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z distance to bilayer center [$\AA$]')
		plt.ylabel('relative particles frequency [$\AA^{-3}$]')
		
		#save figure
		ax.set_xlim(0, 2*bins_nb)
		ax.set_ylim(0, max_density_particles_pc)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
		ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
		xlabel = ax.get_xticks().tolist()
		for tick_index in range(0,len(xlabel)):
			xlabel[tick_index] -= bins_nb
		ax.set_xticklabels(xlabel)
		ax.xaxis.labelpad = 20
		ax.yaxis.labelpad = 20
		plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
		plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
		plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:
			#filenames
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/1_particles/png/2_1_density_profile_particles_' + str(groups_labels[g_index]) +'.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/1_particles/2_1_density_profile_particles_' + str(groups_labels[g_index]) +'.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			fig.suptitle("Relative particles frequency along z around TM clusters of sizes " + str(groups_labels[g_index]))
	
			#plot data
			ax = fig.add_subplot(111)
			for part in particles_def["labels"]:
				plt.plot(range(0,2*bins_nb), density_particles_groups_pc[g_index][part], color = particles_def["colour"][part], label = str(part))
			plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, 0, max_density_particles_pc, linestyles = 'dashed')
			plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, 0, max_density_particles_pc, linestyles = 'dashed')
			plt.vlines(bins_nb, 0, max_density_particles_pc, linestyles = 'dashdot')
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z distance to bilayer center [$\AA$]')
			plt.ylabel('relative particles frequency [$\AA^{-3}$]')
			
			#save figure
			ax.set_xlim(0, 2*bins_nb)
			ax.set_ylim(0, max_density_particles_pc)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
			ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
			xlabel = ax.get_xticks().tolist()
			for tick_index in range(0,len(xlabel)):
				xlabel[tick_index] -= bins_nb
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
def density_graph_residues():											#DONE
			
	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:		

		#filenames
		if c_size == "all sizes":
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/3_residues/png/1_3_density_profile_residues_all_sizes.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/3_residues/1_3_density_profile_residues_all_sizes.svg'		
		else:
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/3_residues/png/1_3_density_profile_residues_' + str(c_size) +'.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/3_residues/1_3_density_profile_residues_' + str(c_size) +'.svg'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		if c_size == "all sizes":
			fig.suptitle("Relative residues frequency along z around TM clusters")
		else:
			fig.suptitle("Relative residues frequency along z around TM clusters of size " + str(c_size))

		#plot data
		ax = fig.add_subplot(111)
		for res in residues_def["labels"]:
			plt.plot(range(0,2*bins_nb), density_residues_sizes_pc[c_size][res], color = residues_def["colour"][res], label = str(res))
		plt.plot(range(0,2*bins_nb), density_particles_sizes_pc[c_size]["peptide"], color = particles_def["colour"]["peptide"], label = "peptide")
		if water_pres:
			plt.plot(range(0,2*bins_nb), density_particles_sizes_pc[c_size]["water"], color = particles_def["colour"]["water"], label = "water")
		plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, 0, max_density_residues_pc, linestyles = 'dashed')
		plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, 0, max_density_residues_pc, linestyles = 'dashed')
		plt.vlines(bins_nb, 0, max_density_residues_pc, linestyles = 'dashdot')
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z distance to bilayer center [$\AA$]')
		plt.ylabel('relative residues frequency [$\AA^{-3}$]')
		
		#save figure
		ax.set_xlim(0, 2*bins_nb)
		ax.set_ylim(0, max_density_residues_pc)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
		ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
		xlabel = ax.get_xticks().tolist()
		for tick_index in range(0,len(xlabel)):
			xlabel[tick_index] -= bins_nb
		ax.set_xticklabels(xlabel)
		ax.xaxis.labelpad = 20
		ax.yaxis.labelpad = 20
		plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
		plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
		plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:
			#filenames
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/3_residues/png/1_3_density_profile_residues_' + str(groups_labels[g_index]) +'.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/3_residues/1_3_density_profile_residues_' + str(groups_labels[g_index]) +'.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			if c_size == "all sizes":
				fig.suptitle("Relative residues frequency along z around TM clusters")
			else:
				fig.suptitle("Relative residues frequency along z around TM clusters of sizes " + str(groups_labels[g_index]))
	
			#plot data
			ax = fig.add_subplot(111)
			for res in residues_def["labels"]:
				plt.plot(range(0,2*bins_nb), density_residues_groups_pc[g_index][res], color = residues_def["colour"][res], label = str(res))
			plt.plot(range(0,2*bins_nb), density_particles_groups_pc[g_index]["peptide"], color = particles_def["colour"]["peptide"], label = "peptide")
			if water_pres:
				plt.plot(range(0,2*bins_nb), density_particles_groups_pc[g_index]["water"], color = particles_def["colour"]["water"], label = "water")
			plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, 0, max_density_residues_pc, linestyles = 'dashed')
			plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, 0, max_density_residues_pc, linestyles = 'dashed')
			plt.vlines(bins_nb, 0, max_density_residues_pc, linestyles = 'dashdot')
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z distance to bilayer center [$\AA$]')
			plt.ylabel('relative residues particles frequency [$\AA^{-3}$]')
			
			#save figure
			ax.set_xlim(0, 2*bins_nb)
			ax.set_ylim(0, max_density_residues_pc)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
			ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
			xlabel = ax.get_xticks().tolist()
			for tick_index in range(0,len(xlabel)):
				xlabel[tick_index] -= bins_nb
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
def density_graph_charges():											#DONE
			
	#sizes
	#=====
	for c_size in sizes_sampled + ["all sizes"]:
		#filenames
		if c_size == "all sizes":
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/2_charges/png/1_2_density_profile_charges_all_sizes.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/2_charges/1_2_density_profile_charges_all_sizes.svg'
		else:
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/2_charges/png/1_2_density_profile_charges_' + str(c_size) +'.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/1_sizes/2_charges/1_2_density_profile_charges_' + str(c_size) +'.svg'

		#create figure
		fig = plt.figure(figsize=(8, 6.2))
		if c_size == "all sizes":
			fig.suptitle("Charge density profile in TM clusters")
		else:
			fig.suptitle("Charge density profile TM clusters of size " + str(c_size))

		#plot data
		ax = fig.add_subplot(111)
		for charge_g in charges_groups.keys() + ["total"]:
			if charge_g == "total":
				plt.plot(range(0,2*bins_nb), density_charges_sizes[c_size][charge_g], color = charges_colours[charge_g], label = str(charge_g), linewidth = 2)
			else:
				plt.plot(range(0,2*bins_nb), density_charges_sizes[c_size][charge_g], color = charges_colours[charge_g], label = str(charge_g))
		plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, min_density_charges, max_density_charges, linestyles = 'dashed')
		plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, min_density_charges, max_density_charges, linestyles = 'dashed')
		plt.vlines(bins_nb, min_density_charges, max_density_charges, linestyles = 'dashdot')
		plt.hlines(0, 0, 2*bins_nb)
		fontP.set_size("small")
		ax.legend(prop=fontP)
		plt.xlabel('z distance to bilayer center [$\AA$]')
		plt.ylabel('average charge density [$e.\AA^{-3}$]')
		
		#save figure
		ax.set_xlim(0, 2*bins_nb)
		ax.set_ylim(min_density_charges, max_density_charges)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
		ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
		xlabel = ax.get_xticks().tolist()
		for tick_index in range(0,len(xlabel)):
			xlabel[tick_index] -= bins_nb
		ax.set_xticklabels(xlabel)
		ax.xaxis.labelpad = 20
		ax.yaxis.labelpad = 20
		plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
		plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
		plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
		fig.savefig(filename_png)
		fig.savefig(filename_svg)
		plt.close()

	#groups
	#======
	if args.cluster_groups_file != "no":
		for g_index in groups_sampled:
			#filenames
			filename_png = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_charges/png/2_2_density_profile_charges_' + str(groups_labels[g_index]) +'.png'
			filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/2_groups/2_charges/2_2_density_profile_charges_' + str(groups_labels[g_index]) +'.svg'
	
			#create figure
			fig = plt.figure(figsize=(8, 6.2))
			if g_index == "all groups":
				fig.suptitle("Charge density profile in TM clusters")
			else:
				fig.suptitle("Charge density profile in TM clusters of sizes " + str(groups_labels[g_index]))
	
			#plot data
			ax = fig.add_subplot(111)
			for charge_g in charges_groups.keys() + ["total"]:
				if charge_g == "total":
					plt.plot(range(0,2*bins_nb), density_charges_groups[g_index][charge_g], color = charges_colours[charge_g], label = str(charge_g), linewidth = 2)
				else:
					plt.plot(range(0,2*bins_nb), density_charges_groups[g_index][charge_g], color = charges_colours[charge_g], label = str(charge_g))
			plt.vlines(np.floor(z_upper/float(args.slices_thick)) + bins_nb, min_density_charges, max_density_charges, linestyles = 'dashed')
			plt.vlines(np.floor(z_lower/float(args.slices_thick)) + bins_nb, min_density_charges, max_density_charges, linestyles = 'dashed')
			plt.vlines(bins_nb, min_density_charges, max_density_charges, linestyles = 'dashdot')
			plt.hlines(0, 0, 2*bins_nb)
			fontP.set_size("small")
			ax.legend(prop=fontP)
			plt.xlabel('z distance to bilayer center [$\AA$]')
			plt.ylabel('average charge density [$e.\AA^{-3}$]')
			
			#save figure
			ax.set_xlim(0, 2*bins_nb)
			ax.set_ylim(min_density_charges, max_density_charges)
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
			ax.xaxis.set_ticks_position('bottom')
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
			ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
			xlabel = ax.get_xticks().tolist()
			for tick_index in range(0,len(xlabel)):
				xlabel[tick_index] -= bins_nb
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
#process inputs
#=========================================================================================

#data loading
set_lipids_beads()
set_particles()
if args.residuesfilename != "no":
	set_residues()
if args.chargesfilename != "no":
	set_charges()
load_MDA_universe()
if args.selection_file_ff != "no":
	identify_ff()
if protein_pres:
	identify_proteins()
identify_leaflets()
if args.cluster_groups_file != "no":
	initialise_groups()

#create data structures
struct_particles()
if args.residuesfilename != "no":
	struct_residues()
if args.chargesfilename != "no":
	struct_charges()

#=========================================================================================
# generate data
#=========================================================================================
print "\nCalculating density profiles..."

#case: structure only
#--------------------
if args.xtcfilename=="no":
	calculate_density(U.trajectory.ts.dimensions)

#case: browse xtc frames
#-----------------------
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
		sys.stdout.flush()
		sys.stdout.write(progress)							
		calculate_density(U.trajectory.ts.dimensions)
	print ""

#=========================================================================================
# process data
#=========================================================================================
calculate_stats()

#=========================================================================================
# produce outputs
#=========================================================================================

print "\nWriting outputs..."

#case: proteins
#--------------
if protein_pres:
	if np.size(sizes_sampled) > 0:
		density_write_particles()
		density_graph_particles()
		if args.residuesfilename != "no":
			density_write_residues()
			density_graph_residues()
		if args.chargesfilename != "no":
			density_write_charges()
			density_graph_charges()
	else:
		print "Warning: no TM cluster identified!"

#case: no proteins
#-----------------
else:
	density_write_particles()
	density_graph_particles()
	if args.chargesfilename != "no":
		density_write_charges()
		density_graph_charges()
	
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
