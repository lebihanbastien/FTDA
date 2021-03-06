# Configuration file for OOFTDA

#--------------------------------------------------------------
# Source the parameters
#--------------------------------------------------------------
source config/parameters.sh

#------------------------------------
# TYPE OF COMPUTATION (COMPTYPE)
# QBTBP   = 0
# NFO2    = 1
# PM      = 2
# PM_TEST = 3
# COMPMAP = 4
# COC     = 5
# TRAJ    = 6
#------------------------------------
COMPTYPE=$PM_TEST

#------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2
#------------------------------------
MODEL=$CRTBP

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
LI_EM=1
LI_SEM=2

#------------------------------------
# TYPE OF MANIFOLD
# MAN_CENTER    = 0
# MAN_CENTER_S  = 1
# MAN_CENTER_U  = 2
# MAN_CENTER_US = 3
#------------------------------------
MANTYPE_EM=$MAN_CENTER_U
MANTYPE_SEM=$MAN_CENTER_U

#------------------------------------
# PM STYLE
# GRAPH = 0; NORMFORM = 1; MIXED = 2
#------------------------------------
PMS=$MIXED

#------------------------------------
# STORAGE
# FALSE = 0; TRUE = 1
#------------------------------------
STORAGE=$FALSE

#------------------------------------
# Point to test
#------------------------------------
IC=(0.1 0.1 0.2 0.0 1e-6)


#------------------------------------
# Order to to test
#------------------------------------
ORDERS=(5 10 12 15)
NORDERS=${#ORDERS[@]} 

#------------------------------------
# Numerical constants for all computation
#------------------------------------
OFTS_ORDER=15





