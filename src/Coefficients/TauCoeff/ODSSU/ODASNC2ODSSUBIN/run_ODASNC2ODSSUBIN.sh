#!/bin/sh

# $Id$

usage()
{
  echo " Usage: run_ODASNC2ODSSUBIN.sh -l|b[fh]"
  echo
  echo "   Convert any netCDF format ODAS SSU TauCoeff files in the current"
  echo "   directory to Binary format."
  echo
  echo "    l       Produce little-endian binary output files"
  echo
  echo "    b       Produce big-endian binary output files"
  echo
  echo "    f       Force overwrite of output file if it already exists."
  echo "            Default behaviour is to skip conversion if the output"
  echo "            file is already present."
  echo
  echo "    h       Print this message and exit"
  echo
  echo "   Note the endian-ness option is only setup for the following"
  echo "   compilers that allow for run-time conversion via environment"
  echo "   variables:"
  echo "     - AIX xlf"
  echo "     - Linux gfortran"
  echo "     - Linux ifort (Intel)"
  echo "     - Linux g95"
  echo
}


# Define defaults
SCRIPT_NAME="`basename $0`"
ENDIAN_TYPE="NONE"
OVERWRITE="NO"
REMOVE="rm -f"


# Parse command line options
while getopts :hlbf OPTVAL; do

  # If option argument looks like another option exit the loop
  case ${OPTARG} in
    -*) break;;
  esac

  # Parse the valid options here
  case ${OPTVAL} in
    l)  ENDIAN_TYPE="little";;
    b)  ENDIAN_TYPE="big";;
    f)  OVERWRITE="YES";;
    h)  usage
        exit 0;;
    :|\?) OPTVAL=${OPTARG}
          break;;
  esac

done

# Remove the options processed
shift `expr ${OPTIND} - 1`

# Now output invalidities based on OPTVAL
# Need to do this as getopts does not handle
# the situations where an option is passed
# as an argument to another option.
case ${OPTVAL} in

  # If OPTVAL contains nothing, then all options
  # have been successfully parsed
  \?) if [ $# -ne 0 ]; then
        ( echo "${SCRIPT_NAME}: Invalid argument(s) $*" ; echo ; usage ) | more
        exit 2
      fi;;

  # Invalid option
  ?) ( echo "${SCRIPT_NAME}: Invalid option '-${OPTARG}'" ; usage ) | more
     exit 2;;

esac


# Check endian type for run time options
# ...Save the current envar options
XLFRTEOPTS_SAVE=${XLFRTEOPTS}
GFORTRAN_CONVERT_UNIT_SAVE=${GFORTRAN_CONVERT_UNIT}
F_UFMTENDIAN_SAVE=${F_UFMTENDIAN}
G95_ENDIAN_SAVE=${G95_ENDIAN}
# ...Set the non-switchable envars
case ${ENDIAN_TYPE} in
  "little") export XLFRTEOPTS="ufmt_littleendian=-100";;
  "big") export XLFRTEOPTS="";;
  *) ( echo "${SCRIPT_NAME}: Must specify and endian type, -l or -b" ; echo ; usage ) | more
     exit 2;;
esac
# ...Switchable envars
export GFORTRAN_CONVERT_UNIT="${ENDIAN_TYPE}_endian"
export F_UFMTENDIAN="${ENDIAN_TYPE}"
export G95_ENDIAN="${ENDIAN_TYPE}"


# Assign processing parameters
EXE_FILE="ODASNC2ODSSUBIN"
LOG_FILE="${EXE_FILE}.${ENDIAN_TYPE}_endian.log"



# Process netCDF CloudCoeff files
SSU_SENSOR_ID_LIST="ssu_tirosn ssu_n06 ssu_n07 ssu_n08 ssu_n09 ssu_n11 ssu_n14"
for SENSOR_ID in ${SSU_SENSOR_ID_LIST}; do

  # Create filenames
  INPUT_FILE="${SENSOR_ID}.input.txt"
  BIN_FILE="${SENSOR_ID}.TauCoeff.bin"
  SIGNAL_FILE="${BIN_FILE}.signal"

  # Check to see if the output file exists
  if [ -f ${BIN_FILE} -a ${OVERWRITE} = "NO" ]; then
    echo " Output file ${BIN_FILE} already exists. Skipping to next sensor..."
    continue
  fi

  # Delete signal file if it exists
  if [ -f ${SIGNAL_FILE} ]; then
    ${REMOVE} ${SIGNAL_FILE}
    if [ $? -ne 0 ]; then
      echo "   ERROR deleting signal file ${SIGNAL_FILE}. Skipping to next sensor..."
      continue
    fi
  fi
  
  # Convert the file
  echo " Converting ${SENSOR_ID} netCDF files to ${ENDIAN_TYPE} endian binary file ${BIN_FILE}..."
  ${EXE_FILE} < ${INPUT_FILE} >> ${LOG_FILE}

  # Check that signal file was created
  if [ -f ${SIGNAL_FILE} ]; then
    ${REMOVE} ${SIGNAL_FILE} > /dev/null
  else
    echo "   ERROR creating output file ${BIN_FILE}."
  fi

done


# Restore the run-time option environment variables
export XLFRTEOPTS="${XLFRTEOPTS_SAVE}"
export GFORTRAN_CONVERT_UNIT="${GFORTRAN_CONVERT_UNIT_SAVE}"
export F_UFMTENDIAN="${F_UFMTENDIAN_SAVE}"
export G95_ENDIAN="${G95_ENDIAN_SAVE}"
