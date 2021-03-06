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
CS=$SEM  

#------------------------------------
# Normalization
# FALSE = 0; TRUE = 1
#------------------------------------
ISNORM=$TRUE 

#------------------------------------
# DEFAULT LIBRATION POINT FOR EM & SEM SYSTEM
#------------------------------------
LI_EM=2
LI_SEM=1

#------------------------------------
# PM STYLE
# GRAPH = 0; NORMFORM = 1; MIXED = 2
#------------------------------------
PMS=$GRAPH

#------------------------------------
# TYPE OF MANIFOLD
# MAN_CENTER    = 0
# MAN_CENTER_S  = 1
# MAN_CENTER_U  = 2
# MAN_CENTER_US = 3
#------------------------------------
MANTYPE_EM=$MAN_CENTER
MANTYPE_SEM=$MAN_CENTER


#------------------------------------
# STORAGE
# FALSE = 0; TRUE = 1
#------------------------------------
STORAGE=$FALSE

#------------------------------------
# PMAP TYPE
# PMAP = 1; TMAP = 2; EMAP = 3; IMAP=4
#------------------------------------
PMAP_TYPE=$IMAP

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
NUM_THREADS=50

#------------------------------------
# Numerical constants for all computation
#------------------------------------
OFS_ORDER=30
OFTS_ORDER=30
REDUCED_NV=4

#------------------------------------
# PMAP PARAMETERS
#------------------------------------
PMAP_order=$OFTS_ORDER       #ORDER
PMAP_ofs_order=$OFS_ORDER    #OFS ORDER

PMAP_TF=10000	     #final integration time. needs to be big for PMAP
PMAP_isQBCP=$MODEL   #isQBCP
PMAP_max_events=100  #Maximum number of events allowed (warning: all directions of crossing are considered!)
PMAP_t0=0.0          #Initial time (-1 for inner computation in C routine)
PMAP_dHv=0.1         #Energy
PMAP_gsize=10000     #Number of steps on the grid
PMAP_gmin=-0.5       #left boundary
PMAP_gmax=0.5        #right boundary
PMAP_PROJ_FREQ=-1    #projection frequency on manifold. If -1, no projection


#------------------------------------
# Order to test
#------------------------------------
ORDERS=(5 10 15 20 25 30)
#ORDERS=(3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30)
NORDERS=${#ORDERS[@]} 

#------------------------------------
# OFS Order to test
#------------------------------------
ORDERS_OFS=($OFS_ORDER)
NORDERS_OFS=${#ORDERS_OFS[@]} 



