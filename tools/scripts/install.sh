#!/bin/bash

usage(){
   echo ""
   echo "Description:"
   echo "                Build, Test and Install MQS apps"
   echo ""
   echo "Usage:"
   echo "                install.sh [ <option> ] ... ]"
   echo ""
   echo "Options:"
   echo "-j <n>          Specify the number of nodes to be used for compilation."
   echo "-c <cxx>        Specify the CXX compiler to be used (either g++ or clang++). Default is clang++."
   echo "-b <builddir>   Specify the build directory (if VScode detected set by default)."
   echo "-i <installdir> Specify the install directory (if VScode detected set by default)."
   echo "-t              Enable test"
   echo "-d              Disable update of VCS"
   echo "-w              Enable WEB server for testing docs (use default port 8080)"
   echo ""
   echo "-h              Prints this help information"
   echo ""
   exit 1
}

#########################################################
## parse parameters
##########################################################
while getopts "hj:s:c:b:dti:w" option ; do
   case $option in
       h ) usage ;;
       j ) NP=$OPTARG ;;
       c ) CXX=$OPTARG ;;
       s ) SRCDIR=$OPTARG ;;
       b ) BUILDDIR=$OPTARG ;;
       i ) INSTALLDIR=$OPTARG ;;
       d ) DEBUG=1 ;;
       t ) TESTON=1 ;;
       w ) WEBON=1 ;;
       ? ) usage ;;
   esac
done
# shift to have the good number of other args
shift $((OPTIND - 1))


: ${DEBUG:=0}
: ${TESTON:=0}
: ${WEBON:=0}
: ${NP:=0}
: ${NP_MAX:=$(/usr/bin/getconf _NPROCESSORS_ONLN)}
: ${CXX:=clang++-9}
: ${CC:=clang-9}
: ${FC:=gfortran}

if [ ${NP} -eq 0 ]; then
   NJOBS=1
fi
CURDIR=$PWD

if [ "${TERM_PROGRAM}" = "vscode" ]; then
   echo "**** VScode detected: overwrite main directories definition ****"
   SRCDIR=/workspaces/mqs
   BUILDDIR=/workspaces/mqs/build
   INSTALLDIR=/workspaces/mqs/install
fi

# check if SRCDIR  exists
if [ -z ${SRCDIR} ]; then
   echo "SRCDIR is undefined"
   exit 1
fi

if [ ! -d ${SRCDIR} ]; then
   echo "SRCDIR=${SRCDIR} is not defined"
   exit 1
fi

# check if NJOBS is 
if [ ${NJOBS} -ge ${NP_MAX} ]; then
   echo "NJOBS=${NJOBS} is invalid - should be less than or equal to ${NP_MAX}"
   exit 1
fi

# check if CXX exists

# actually build
mkdir -p ${BUILDDIR}
cd  ${BUILDDIR}

cmake ${SRCDIR} \
 -DCMAKE_CXX_COMPILER:FILEPATH=/usr/bin/${CXX} \
 -DCMAKE_C_COMPILER:FILEPATH=/usr/bin/${CC} \
 -DCMAKE_fortran_COMPILER:FILEPATH=/usr/bin/${FC} \
 -DCMAKE_VERBOSE_MAKEFILE=ON \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_INSTALL_PREFIX=${INSTALLDIR} > cmake.log 2>&1

make -j${NJOBS} > make.log 2>&1
make install > install.log 2>&1

if [ "${TESTON}" = "1" ]; then
   mkdir check
   cd check
   mpirun --bind-to core --mca btl vader,self -np ${NJOBS} ${INSTALLDIR}/bin/feelpp_p_myapp 
fi

if [ "${WEBON}" = "1" ]; then
   # get VERSION from main mqs CMakeLists.txt
   VERSION="0.1.0"
   echo "************************************************************"
   echo "Starting WEB server on port 8080"
   sudo npm i -g live-server
   if [ "${TERM_PROGRAM}" = "vscode" ]; then 
      echo "To access the docs:"
      echo "Tape F1 and Forward port 8080"
   fi
   nohup live-server --wait=1000 ${INSTALLDIR}/share/doc/mqs > output.log &
   
   # get PID
   WEBPID=$(ps -ef|grep /usr/local/bin/live-server | head -1 |  tr -s ' ' | cut -d' ' -f 2)
   
   echo "To view the docs; run Firefox in private mode and load localhost:8080/mqs/${VERSION}"
   echo "If you make change to the docs, do not forget to run make docs in ${BUILDDIR}"
   echo "To stop the server: kill $WEBPID"
   echo "************************************************************"
fi

cd ${CURDIR} 
