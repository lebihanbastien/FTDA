# script file for OOFTDA
# RECALL: make it executable with "$ chmod +x script.sh"

#--------------------------------------------------------------
# Source the configuration file passed as argument
#--------------------------------------------------------------
source $1

#--------------------------------------------------------------
# Title
#--------------------------------------------------------------
echo
echo "#------------------------------------------#"
echo "#          OOFTDA configuration            #"
echo "#------------------------------------------#"
echo "Source file: " "$1"
#--------------------------------------------------------------
# Additional tests
#--------------------------------------------------------------

#----------------------------------------
# OFS_ORDER = Order of the Fourier series
#----------------------------------------
#Check that the variable $MODEL exist
if [ -z ${MODEL+x} ]; then
	#If not, default values to QBCP
	echo 'WARNING: the variable MODEL is not set.'
	echo 'QBCP is used by default.'
	MODEL=$QBCP
else
	if [ "$MODEL" == "$RTBP" ]; then
		OFS_ORDER=0
	else  
		OFS_ORDER=30
	fi
fi


#----------------------------------------
# MANTYPE    = Type of manifold
# REDUCED_NV = Number of reduced variables
#----------------------------------------
#Check that the variable $MANTYPE_EM exist
if [ -z ${MANTYPE_EM+x} ]; then
	#If not, default values to center manifold
	echo 'WARNING: the variable MANTYPE_EM is not set.'
	echo 'By default, MANTYPE_EM is set to MAN_CENTER and REDUCED_NV to 4.'
	MANTYPE=$MAN_CENTER
	REDUCED_NV=4
else
	#If yes, set the REDUCED_NV in correspondance
	case $MANTYPE_EM in
		$MAN_CENTER)   
			        REDUCED_NV=4
				#Check that PMS exists
				if [ -z ${PMS+x} ]; then
					#If not, send a warning and use GRAPH by default
					echo 'WARNING: the variable PMS is not set.'
					echo 'GRAPH is used by default.'
					PMS=$GRAPH
				fi	
			;;
		$MAN_CENTER_S)  
			        REDUCED_NV=5
			        PMS=$MIXED
			;;
		$MAN_CENTER_U)  
				REDUCED_NV=5
			        PMS=$MIXED
			;;
		$MAN_CENTER_US) REDUCED_NV=6
			        PMS=$MIXED
			;;
	esac
fi

#Check that the variable $MANTYPE_SEM exist
if [ -z ${MANTYPE_SEM+x} ]; then
	#If not, default values to center manifold
	echo 'WARNING: the variable MANTYPE_SEM is not set.'
	echo 'By default, MANTYPE_SEM is set to MAN_CENTER and REDUCED_NV to 4.'
	MANTYPE=$MAN_CENTER
	REDUCED_NV=4
else
	#If yes, set the REDUCED_NV in correspondance
	case $MANTYPE_SEM in
		$MAN_CENTER)   
			        REDUCED_NV=4
				#Check that PMS exists
				if [ -z ${PMS+x} ]; then
					#If not, send a warning and use GRAPH by default
					echo 'WARNING: the variable PMS is not set.'
					echo 'GRAPH is used by default.'
					PMS=$GRAPH
				fi	
			;;
		$MAN_CENTER_S)  
			        REDUCED_NV=5
			        PMS=$MIXED
			;;
		$MAN_CENTER_U)  
				REDUCED_NV=5
			        PMS=$MIXED
			;;
		$MAN_CENTER_US) REDUCED_NV=6
			        PMS=$MIXED
			;;
	esac
fi


#--------------------------------------------------------------
# Display current set of parameters
#--------------------------------------------------------------
echo "#------------------------------------------#"
echo "# PARAMETERS: "
echo "#------------------------------------------#"
# MANIFOLD TYPE (EM)
#If yes, set the REDUCED_NV in correspondance
case $MANTYPE_EM in
	$MAN_CENTER)    echo 'MANTYPE_EM   = MAN_CENTER'
			;;
	$MAN_CENTER_S)  echo 'MANTYPE_EM   = MAN_CENTER_S'
			;;
	$MAN_CENTER_U)  echo 'MANTYPE_EM   = MAN_CENTER_U'
			;;
	$MAN_CENTER_US) echo 'MANTYPE_EM   = MAN_CENTER_US'
			;;
	*) 	        echo "MANTYPE_EM   = " $MANTYPE_EM". Unknown type." 
esac

# MANIFOLD TYPE (SEM)
#If yes, set the REDUCED_NV in correspondance
case $MANTYPE_SEM in
	$MAN_CENTER)    echo 'MANTYPE_SEM  = MAN_CENTER'
			;;
	$MAN_CENTER_S)  echo 'MANTYPE_SEM  = MAN_CENTER_S'
			;;
	$MAN_CENTER_U)  echo 'MANTYPE_SEM  = MAN_CENTER_U'
			;;
	$MAN_CENTER_US) echo 'MANTYPE_SEM  = MAN_CENTER_US'
			;;
	*) 	        echo "MANTYPE_SEM  = " $MANTYPE_SEM". Unknown type." 
esac
echo

# TYPE OF COMPUTATION
case $COMPTYPE in
	$QBTBP)   echo 'COMPTYPE = QBTBP'
	;;
	$NFO2) 	  echo 'COMPTYPE = NFO2'
	;;
	$PM)      echo 'COMPTYPE = PM'
	;;
	$PM_TEST) echo 'COMPTYPE = PM_TEST'
	;;
	$COMPMAP) echo 'COMPTYPE = COMPMAP'
	;;
	$COC)     echo 'COMPTYPE = COC'
	;;
	$TRAJ)    echo 'COMPTYPE = TRAJ'
	;;
	*) 	  echo "COMPTYPE = " $COMPTYPE ". Unknown type."
esac

# MODEL
case $MODEL in
	$RTBP)  echo 'MODEL    = RTBP'
	;;
	$QBCP) 	echo 'MODEL    = QBCP'
	;;
	$BCP)   echo 'MODEL    = BCP'
	;;
	$ERTBP) echo 'MODEL    = ERTBP'
	;;
	*)      echo "MODEL    = " $MODEL ". Unknown type."
esac


# COORDINATE SYSTEM
case $CS in
	$EM)  echo 'CS       = EM'
	;;
	$SEM) echo 'CS       = SEM'
	;;
	$SE)  echo 'CS       = SE'
	;;
	*)    echo "CS       = " $CS ". Unknown type."
esac


# NORMALIZATION
echo "ISNORM   =" $ISNORM


# DEFAULT LIBRATION POINT
echo "LI_EM    =" $LI_EM
echo "LI_SEM   =" $LI_SEM


# PM STYLE
case $PMS in
	$GRAPH)    echo 'PMS      = GRAPH'
	;;
	$NORMFORM) echo 'PMS      = NORMFORM'
	;;
	$MIXED)    echo 'PMS      = MIXED'
	;;
	*)         echo "PMS      = " $PMS ". Unknown type."
esac


# STORAGE
echo "STORAGE  =" $STORAGE
echo

#--------------------------------------------------------------
# Display initial conditions if necessary: PM TESTING
#--------------------------------------------------------------
if [ $COMPTYPE == $PM_TEST ]; then

	#Check the number of initial conditions match REDUCED_NV
	if [ ${#IC[@]} != $REDUCED_NV ]; then
		echo 'WARNING: the number of IC does not match REDUCED_NV. Stop.'
		exit 1
	fi

	echo "#------------------------------------------#"
	echo "# INITIAL CONDITIONS FOR PM TESTING:"
	echo "#------------------------------------------#"
	echo "si = [" ${IC[*]} "]"
	echo

	echo "#------------------------------------------#"
	echo "# ORDERS FOR PM TESTING:"
	echo "#------------------------------------------#"
	echo "NORDERS = " $NORDERS
	echo "ORDERS  = [" ${ORDERS[*]} "]"
	echo
	
fi


#--------------------------------------------------------------
# Display initial conditions if necessary: PMAP
#--------------------------------------------------------------
if [ $COMPTYPE == $COMPMAP ]; then

	echo "#------------------------------------------#"
	echo "# ORDERS FOR PMAP:"
	echo "#------------------------------------------#"
	echo "NORDERS = " $NORDERS
	echo "ORDERS  = [" ${ORDERS[*]} "]"
	echo

	echo "#------------------------------------------#"
	echo "# OFS ORDERS FOR PMAP:"
	echo "#------------------------------------------#"
	echo "NORDERS_OFS =  " $NORDERS_OFS
	echo "ORDERS_OFS  = [" ${ORDERS_OFS[*]} "]"
	echo

	echo "#------------------------------------------#"
	echo "# PMAP TYPE:"
	echo "#------------------------------------------#"
	case $PMAP_TYPE in
		$PMAP) echo 'PMAP_TYPE = PMAP'
		;;
		$TMAP) echo 'PMAP_TYPE = TMAP'
		;;
		$EMAP) echo 'PMAP_TYPE = EMAP'
		;;
		$IMAP) echo 'PMAP_TYPE = IMAP'
		;;
		$HMAP) echo 'PMAP_TYPE = HMAP'
		;;
		*)     echo "PMAP_TYPE = " $PMAP_TYPE ". Unknown type."
	esac
	echo

	echo "#------------------------------------------#"
	echo "# PMAP PARAMETERS:"
	echo "# Note: t0 = -1 for inner computation "
	echo "# in C routine."
	echo "#------------------------------------------#"
	echo "tf         =" $PMAP_TF
	echo "t0         =" $PMAP_t0
	echo "isQBCP     =" $PMAP_isQBCP
	echo "order      =" $PMAP_order
	echo "ofs_order  =" $PMAP_ofs_order
	echo "max_events =" $PMAP_max_events
	echo "dHv        =" $PMAP_dHv
	echo "gsize      =" $PMAP_gsize
	echo "gmin       =" $PMAP_gmin
	echo "gmax       =" $PMAP_gmax
	echo

	echo "#------------------------------------------#"
	echo "# PMAP SETTINGS:"
	echo "#------------------------------------------#"
	echo "PMAP_append =" $PMAP_append
	echo "PMAP_isPlot =" $PMAP_isPlot
	echo "PMAP_isPar  =" $PMAP_isPar
	echo

	echo "#------------------------------------------#"
	echo "# PMAP METHODS:"
	echo "#------------------------------------------#"
	case $PMAP_method in
		$DUAL_INT) 	    echo 'PMAP_method = DUAL_INT'
		;;
		$DUAL_INT_NO_RESET) echo 'PMAP_method = DUAL_INT_NO_RESET'
		;;
		$DUAL_INT_STEPPED)  echo 'PMAP_method = DUAL_INT_STEPPED'
		;;
		$SINGLE_INT)        echo 'PMAP_method = SINGLE_INT'
		;;
		*)                  echo "PMAP_method =" $PMAP_method ". Unknown type."
	esac
	echo

	echo "#------------------------------------------#"
	echo "# OpenMP SETTINGS:"
	echo "#------------------------------------------#"
	echo "NUM_THREADS =" $NUM_THREADS
	echo

	echo "#------------------------------------------#"
	echo "# Frequency of projection on the MAN:"
	echo "#------------------------------------------#"
	echo "PMAP_PROJ_FREQ =" $PMAP_PROJ_FREQ
	echo	
fi



echo "#------------------------------------------#"
echo "# GLOBAL CONSTANTS:"
echo "#------------------------------------------#"
echo "OFS_ORDER  =" $OFS_ORDER
echo "OFTS_ORDER =" $OFTS_ORDER
echo "REDUCED_NV =" $REDUCED_NV
echo

#--------------------------------------------------------------
# Go on with the implementation?
#--------------------------------------------------------------
echo -e "Do you want to go on with the computation (y/n)? \c "
read  ans
# bash check the answer
if [ "$ans" == "y" ]; then
	
	#-------------------------------
	#Build the array of coefficients
	#-------------------------------
	case $COMPTYPE in
	
	$QBTBP) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE)
		;;

	$NFO2) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE)
		;;

	$PM) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE)
		;;

	$COC) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE)
		;;

	$PM_TEST) 
		COEFFS=($OFTS_ORDER $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE ${IC[*]} $NORDERS ${ORDERS[*]})
		;;

	

	$COMPMAP) 
		#Building the PMAP vector
		PMAP=($PMAP_PROJ_FREQ $PMAP_TYPE $PMAP_TF $PMAP_isQBCP $PMAP_order $PMAP_ofs_order)
 	        PMAP=(${PMAP[*]} $PMAP_max_events $PMAP_t0 $PMAP_dHv $PMAP_gsize $PMAP_gmin $PMAP_gmax)
		PMAP=(${PMAP[*]} $PMAP_append $PMAP_isPlot $PMAP_isPar $PMAP_method)
		#Number of PMAP parameters 
		NPMAP=${#PMAP[@]} 
		#Coeffs
		COEFFS=($OFTS_ORDER  $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM)
		COEFFS=(${COEFFS[*]} $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE)
		COEFFS=(${COEFFS[*]} $NORDERS ${ORDERS[*]} )
		COEFFS=(${COEFFS[*]} $NORDERS_OFS ${ORDERS_OFS[*]} )
		COEFFS=(${COEFFS[*]} $NPMAP ${PMAP[*]} $NUM_THREADS ) 
		;;
	$TRAJ) 
		#Building the PMAP vector
		PMAP=($PMAP_TYPE $PMAP_TF $PMAP_FACTOR $PMAP_isQBCP)
 	        PMAP=(${PMAP[*]} $PMAP_max_events $PMAP_t0 $PMAP_dHv)
		#Number of PMAP parameters 
		NPMAP=${#PMAP[@]} 
		#Coeffs
		COEFFS=($OFTS_ORDER  $OFS_ORDER $REDUCED_NV $COMPTYPE $MODEL $CS $ISNORM)
		COEFFS=(${COEFFS[*]} $LI_EM $LI_SEM $PMS $MANTYPE_EM $MANTYPE_SEM $STORAGE)
		COEFFS=(${COEFFS[*]} $NPMAP ${PMAP[*]} $NUM_THREADS ${IC[*]} )	   
		;;
	esac

	#-------------------------------
	#Call software
	#-------------------------------
	#Check the NOHUP condition
	if [ "$ISNOHUP" == "1" ]; then
			#If true, check that an output file has been set
			if [ -z ${OUT+x} ]; then
				echo 'WARNING: OUT variable has not been set.'
				echo 'OUT = default.out by default.'
				OUT='default.out'
			fi
			#Call software
			nohup bin/Release/OOFTDA ${COEFFS[*]} > $OUT &
	else
			bin/Release/OOFTDA ${COEFFS[*]}
	fi

	
	
else  
	echo "Stop. No computation."
fi 
