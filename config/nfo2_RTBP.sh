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
COMPTYPE=$NFO2

#------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2; ERTBP = 3
#------------------------------------
MODEL=$CRTBP

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
LI_EM=2
LI_SEM=2

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
# Numerical constants for all computation
#------------------------------------
OFS_ORDER=30
OFTS_ORDER=8
REDUCED_NV=4



