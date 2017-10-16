# Configuration file for OOFTDA

#--------------------------------------------------------------
# Source the parameters
#--------------------------------------------------------------
source config/parameters.sh

#------------------------------------
# TYPE OF COMPUTATION
# QBTBP   = 0
# NFO2    = 1
# PM      = 2
# PM_TEST = 3
# COMPMAP = 4
#------------------------------------
COMPTYPE=$COMPMAP

#------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#------------------------------------
MODEL=$QBCP 

#------------------------------------
# COORDINATE SYSTEM
# EM = 0; SEM = 1; SE = 2
#------------------------------------
CS=$EM  

#------------------------------------
# Normalization
# FALSE = 0; TRUE = 1
#------------------------------------
ISNORM=$TRUE 

#------------------------------------
# DEFAULT LIBRATION POINT FOR EM & SEM SYSTEM
#------------------------------------
LI_EM=1
LI_SEM=2

#------------------------------------
# PM STYLE
# GRAPH = 0; NORMFORM = 1; MIXED = 2
#------------------------------------
PMS=$GRAPH

#------------------------------------
# STORAGE
# FALSE = 0; TRUE = 1
#------------------------------------
STORAGE=$FALSE

#------------------------------------
# Order to test
#------------------------------------
ORDERS=(20)
NORDERS=${#ORDERS[@]} 

#------------------------------------
# OFS Order to test
#------------------------------------
ORDERS_OFS=(30)
NORDERS_OFS=${#ORDERS_OFS[@]} 

#------------------------------------
# PMAP TYPE
# PMAP = 1; TMAP = 2; EMAP = 3.
#------------------------------------
PMAP_TYPE=$EMAP

#------------------------------------
# PMAP METHODS
# DUAL_INT = 1; DUAL_INT_NO_RESET = 2; 
# DUAL_INT_STEPPED = 3; SINGLE_INT = 4.
#------------------------------------
PMAP_method=$SINGLE_INT

#------------------------------------
# PMAP SETTINGS
#------------------------------------
PMAP_append=0;
PMAP_isPlot=0;
PMAP_isPar=1;

#------------------------------------
# OpenMP settings
#------------------------------------
NUM_THREADS=2

#------------------------------------
# PMAP PARAMETERS
#------------------------------------
PMAP_TF=10000	     #final integration time. needs to be big for PMAP
PMAP_isQBCP=$MODEL   #isQBCP
PMAP_order=20        #ORDER
PMAP_ofs_order=30    #OFS ORDER
PMAP_max_events=1000 #Maximum number of events allowed (warning: all directions of crossing are considered!)
PMAP_t0=0.0          #Initial time (-1 for inner computation in C routine)
PMAP_dHv=0.01        #Energy
PMAP_gsize=10        #Number of steps on the grid
PMAP_gmin=0.0        #left boundary
PMAP_gmax=2.0        #right boundary

#------------------------------------
# Numerical constants for all computation
#------------------------------------
OFS_ORDER=30
OFTS_ORDER=20
REDUCED_NV=4




