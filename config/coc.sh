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
# PMAP    = 4
# COC     = 5
#------------------------------------
COMPTYPE=$COC

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
STORAGE=$TRUE

#------------------------------------
# Numerical constants for all computation
#------------------------------------
OFS_ORDER=0
OFTS_ORDER=8
REDUCED_NV=4






