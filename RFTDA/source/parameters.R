# Config files for parameters of the QBCP/BCP Sun-Earth-Moon system.
# Object Oriented.
# BLB 2016
#
#---------------------------------------------------------

#---------------------------------------------------------
# EML1
#---------------------------------------------------------
EML1 = structure(list(), class = "LIBRATION_POINT")

# From CRTBP
EML1$MU    = +1.215058162343360e-02 #mass ratio [-]
EML1$LPRIM = +384400 # distance between primaries [km]
EML1$TPRIM = +2.360584684800000e+06 # period of rotation of the primaries [s]
EML1$GAMMA = +0.150934272990064 # distance Li-smallest primary [-]
EML1$C1    = -5.5448979798087 # C1 constant [-]
# From QBCP
EML1$TQBCP = 6.791193871907917e+00 #[-]


#---------------------------------------------------------
# EML2
#---------------------------------------------------------
EML2 = structure(list(), class = "LIBRATION_POINT")

# From CRTBP
EML2$MU    = +1.215058162343360e-02 #mass ratio [-]
EML2$LPRIM = +384400 # distance between primaries [km]
EML2$TPRIM = +2.360584684800000e+06 # period of rotation of the primaries [s]
EML2$GAMMA = +0.16783273173707 # distance Li-smallest primary [-]
EML2$C1    = -6.8859163415402 # C1 constant [-]
# From QBCP
EML2$TQBCP = 6.791193871907917e+00 #[-]


#---------------------------------------------------------
# SEML1
#---------------------------------------------------------
SEML1 = structure(list(), class = "LIBRATION_POINT")

# From CRTBP
SEML1$MU    = +3.040277404912506e-06 #mass ratio [-]
SEML1$LPRIM = +149.60e6 # distance between primaries [km]
SEML1$TPRIM = +3.155814950400000e+07 # period of rotation of the primaries [s]
SEML1$GAMMA = +0.0100108175327472 # distance Li-smallest primary [-]
SEML1$C1    = -98.8916378512805 # C1 constant [-]
# From QBCP
SEML1$TQBCP = 5.080085647283307e-01 #[-]


#---------------------------------------------------------
# SEML2
#---------------------------------------------------------
SEML2 = structure(list(), class = "LIBRATION_POINT")

# From CRTBP
SEML2$MU    = +3.040277404912506e-06 #mass ratio [-]
SEML2$LPRIM = +149.60e6 # distance between primaries [km]
SEML2$TPRIM = +3.155814950400000e+07 # period of rotation of the primaries [s]
SEML2$GAMMA = +0.0100780785917628 # distance Li-smallest primary [-]
SEML2$C1    = -100.224961432011 # C1 constant [-]
# From QBCP
SEML2$TQBCP = 5.080085647283307e-01 #[-]


#---------------------------------------------------------
# Routine to select the libration point
#---------------------------------------------------------
WHICH_LIB_POINT <- function(Li, FWRK)
{
  if (FWRK == "EM")
  {
    if (Li == "L1")
      LIB = EML1
    if (Li == "L2")
      LIB = EML2
  }
  if (FWRK == "SEM")
  {
    if (Li == "L1")
      LIB = SEML1
    if (Li == "L2")
      LIB = SEML2
  }
  return(LIB)
}
