#!/bin/sh
# Parameters file for OOFTDA

#------------------------------------
# TYPE OF COMPUTATION (COMPTYPE)
# QBTBP   = 0
# NFO2    = 1
# PM      = 2
# PM_TEST = 3
# COMPMAP = 4
# COC     = 5
# DYNEQ   = 6
#------------------------------------
export QBTBP=0
export NFO2=1
export PM=2
export PM_TEST=3
export COMPMAP=4
export COC=5
export DYNEQ=6
export TRAJ=7
export FT_TEST=8

# List of possible computations for YAD dialog
export LIST_TYPE_OF_COMPUTATIONS=$(echo "QBTBP|NFO2|PM|PM_TEST|DYNEQ|FT_TEST")

read -r -d '' LABEL_TYPE_OF_COMPUTATION << EOM
- QBTBP = solving the Quasi-Bicircular Three-Body Problem.
- NFO2 = computing the Normal Form of Order 2 of the Hamiltonian about a libration point.
- PM = the Parameterization Method.
- PM_TEST = tests on the Parameterization Method.
- DYNEQ   = Dynamical Equivalents to the Libration Points.
- FT_TEST = test of the Fourier-Taylor Algebra.
EOM
export LABEL_TYPE_OF_COMPUTATION

#- COMPMAP = computing various Maps (Invariance error, Poincaré, stroboscopic...)
#- COC = tests on the various Changes Of Coordinates.

#------------------------------------
# TYPE OF MANIFOLD
# MAN_CENTER    = 0
# MAN_CENTER_S  = 1
# MAN_CENTER_U  = 2
# MAN_CENTER_US = 3
#------------------------------------
export MAN_CENTER=0
export MAN_CENTER_S=1
export MAN_CENTER_U=2
export MAN_CENTER_US=3

# List of possible manifolds for YAD dialog
export LIST_MANIFOLD_TYPE=$(echo "MAN_CENTER|MAN_CENTER_S|MAN_CENTER_U|MAN_CENTER_US")

read -r -d '' LABEL_MANIFOLD_TYPE << EOM
MAN_CENTER = Center.
MAN_CENTER_S = Center-Stable.
MAN_CENTER_U = Center-Unstable. 
MAN_CENTER_US = Center-Hyperbolic.:LBL
EOM
export LABEL_MANIFOLD_TYPE


#------------------------------------
# MODEL
# RTBP = 0; QBCP = 1; BCP = 2; ERTBP = 3
#------------------------------------
export RTBP=0
export QBCP=1
export BCP=2  
export ERTBP=3  

#------------------------------------
# COORDINATE SYSTEM
# EM = 0; SEM = 1; SE = 2
#------------------------------------
export EM=0
export SEM=1
export SE=2  

# List of possible coordinate systems for YAD dialog
export LIST_COORDINATE_SYSTEM=$(echo "EM,SEM")

read -r -d '' LABEL_COORDINATE_SYSTEM << EOM
EM = Earth-Moon, SEM = Sun-Earth:LBL
EOM
export LABEL_COORDINATE_SYSTEM


#------------------------------------
# Normalization
# FALSE = 0; TRUE = 1
#------------------------------------
export FALSE=0
export TRUE=1 

#------------------------------------
# PM STYLE
# GRAPH = 0; NORMAL FORM = 1; MIXED = 2
#------------------------------------
export GRAPH=0
export NORMFORM=1
export MIXED=2

# List of possible PM styles for YAD dialog
export LIST_PM_STYLE=$(echo "GRAPH|NORMFORM|MIXED")

read -r -d '' LABEL_PM_STYLE<< EOM
GRAPH = Graph style.
NORMFORM = Normal Form Style.
MIXED = Mixed style.:LBL
EOM
export LABEL_PM_STYLE

#------------------------------------
# STORAGE
# FALSE = 0; TRUE = 1
#------------------------------------
read -r -d '' LABEL_BOOL<< EOM
FALSE = 0, TRUE = 1:LBL
EOM
export LABEL_BOOL

#------------------------------------
# PMAP TYPE
# PMAP = 1; TMAP = 2; EMAP = 3; IMAP = 4
#------------------------------------
export PMAP=1
export TMAP=2
export EMAP=3
export IMAP=4
export HMAP=5
export IMAPPLANAR=6

#------------------------------------
# PMAP METHODS
# DUAL_INT = 1; DUAL_INT_NO_RESET = 2; 
# DUAL_INT_STEPPED = 3; SINGLE_INT = 4.
#------------------------------------
export DUAL_INT=1
export DUAL_INT_NO_RESET=2 
export DUAL_INT_STEPPED=3
export SINGLE_INT=4

#------------------------------------
# LIBRATION POINT
# EML1=1
# EML2=2
# SEL1=3
# SEL2=4
#------------------------------------
export EML1=1
export EML2=2
export SEL1=3
export SEL2=4

# List of possible lib points for YAD dialog
export LIST_LIBRATION_POINTS=$(echo "EML1|EML2|SEL1|SEL2")

