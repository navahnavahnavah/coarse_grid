#!/bin/tcsh
#
# richardd, 7 Oct 11
#
# converted to tcsh, 17 Oct 11
#
# example PBS script for berserker
#
# put this script in the same directory as the program
# you want to run.
#
# set the name of the job
#PBS -N coarse_cons
#
# set the output and error files
#PBS -o /home/navah/coarse_grid/mOutG.txt
#PBS -e /home/navah/coarse_grid/mErrG.txt
#PBS -m abe -M navah@uchicago.edu
# set the number of nodes to use, and number of processors
# to use per node


##PBS -l nodes=compute-0-0:ppn=12+compute-0-1:ppn=12+compute-0-2:ppn=12+compute-0-3:ppn=12+compute-0-4:ppn=12+compute-1-1:ppn=12+compute-1-2:ppn=12+compute-1-3:ppn=12+compute-1-4:ppn=12+compute-1-5:ppn=12+compute-1-6:ppn=12+compute-1-7:ppn=12+compute-1-8:ppn=12+compute-1-9:ppn=12+compute-1-0:ppn=12


##PBS -l nodes=compute-0-1:ppn=10+compute-0-2:ppn=10+compute-0-3:ppn=10+compute-0-4:ppn=10+compute-0-0:ppn=10
#PBS -l nodes=compute-1-1:ppn=10+compute-1-2:ppn=10+compute-1-3:ppn=10+compute-1-4:ppn=10+compute-1-0:ppn=10
##PBS -l nodes=compute-1-5:ppn=10+compute-1-6:ppn=10+compute-1-7:ppn=10+compute-1-8:ppn=10+compute-1-9:ppn=10

# or, if using only one node, you can do it this way too
##PBS -l ncpus=5


# NEW STUFF MAY 2016

# in this example, I'm using the intel compilers and mvapich2
#
# bring in the module settings

sleep 10s 
source /etc/profile.d/modules.csh
module load intel/intel-12
module load mpi/mvapich2/intel

  
 

# model parameters go here i guess
set PARAM_O='100'
set PARAM_O_RHS='-50'
set PARAM_W='300'
set PARAM_W_RHS='25'
set PARAM_H='200'
set PARAM_TSW='2750'
set PARAM_DIC='232.5'
set PARAM_SCOPE='5'

set PARAM_RESTART='1'
set PARAM_CRASHSTEP='2'
set PARAM_TRACE='0'


set PARAM_CH='600'

set PARAM_PAQ = '1e-12' 

set PARAM_CH_RHS = '600'

## FRACTURE PARAMS
set PARAM_F_DX = '-1.6'
set PARAM_F_K = '1e-5'
#set PARAM_F_FREQ = '0.0005'
set PARAM_F_FREQ = '20'
set PARAM_F_POR = '6e-4' 
 



set PARAM_PATH='/home/navah/coarse_grid/output/coarse_cons/'
set PARAM_ISO_PATH='/home/navah/coarse_grid/output/coarse_cons/'

 
echo $PARAM_PATH
# set PROGNAME to the name of your program
set PROGNAME=massacr
 
# figure out which mpiexec to use
set LAUNCH=/usr/mpi/intel/openmpi-1.4.3-qlc/bin/mpirun
 
# working directory
set WORKDIR=${HOME}
set WORKDIR=/home/navah/coarse_grid
 
set NCPU=`wc -l < $PBS_NODEFILE`
set NNODES=`uniq $PBS_NODEFILE | wc -l`
# set this to zero to turn OFF debugging, 1 to turn it on
set PERMDIR=${HOME}
set SERVPERMDIR=${PBS_O_HOST}:${PERMDIR}

set DEBUG=1
if ( $DEBUG ) then 
	echo ------------------------------------------------------
	echo ' This job is allocated on '${NCPU}' cpu(s)'
	echo 'Job is running on node(s): '
	cat $PBS_NODEFILE
	echo ------------------------------------------------------
	echo PBS: qsub is running on $PBS_O_HOST
	echo PBS: originating queue is $PBS_O_QUEUE
	echo PBS: executing queue is $PBS_QUEUE
	echo PBS: working directory is $PBS_O_WORKDIR
	echo PBS: execution mode is $PBS_ENVIRONMENT
	echo PBS: job identifier is $PBS_JOBID
	echo PBS: job name is $PBS_JOBNAME
	echo PBS: node file is $PBS_NODEFILE
	echo PBS: number of nodes is $NNODES
	echo PBS: current home directory is $PBS_O_HOME
	echo PBS: PATH = $PBS_O_PATH
	echo ------------------------------------------------------
	echo workdir is $WORKDIR
	echo permdir is $PERMDIR
	echo servpermdir is $SERVPERMDIR
	echo ------------------------------------------------------
	echo 'Job is running on node(s): '
	cat $PBS_NODEFILE
	echo ------------------------------------------------------
	echo ${LAUNCH} -n {$NCPU} -f ${PBS_NODEFILE} ${WORKDIR}/${PROGNAME}
	echo ' '
	echo ' '
endif

echo $NCPU
echo $PBS_NODEFILE
echo $PARAM_PATH

# run the program

# set SCRDIR=/state/partition1/navah/${PBS_JOBID}/
# mkdir -p $SCRDIR




#cd $SCRDIR
cd ${WORKDIR}

find $PARAM_PATH -name 'z_sec*.txt' -exec rm -f {} \;
find ${PARAM_PATH}'ch_s/' -name '*.txt' -exec rm -f {} \;
find ${PARAM_PATH}'ch_a/' -name '*.txt' -exec rm -f {} \;
find ${PARAM_PATH}'ch_b/' -name '*.txt' -exec rm -f {} \;
find ${PARAM_PATH}'ch_d/' -name '*.txt' -exec rm -f {} \;
# find $PARAM_PATH + 'ch_a/' -name '*.txt' -exec rm -f {} \;
# find $PARAM_PATH + 'ch_b/' -name '*.txt' -exec rm -f {} \;
##wait
${LAUNCH} -n {$NCPU} -hostfile ${PBS_NODEFILE} --nooversubscribe ${WORKDIR}/${PROGNAME} ${PARAM_RESTART} ${PARAM_PATH} ${PARAM_PATH} ${PARAM_CRASHSTEP} ${PARAM_O} ${PARAM_W} ${PARAM_W_RHS} ${PARAM_H} ${PARAM_O_RHS} ${PARAM_TSW} ${PARAM_DIC} ${PARAM_SCOPE} ${PARAM_TRACE} ${PARAM_CH} ${PARAM_PAQ} ${PARAM_CH_RHS} ${PARAM_F_DX} ${PARAM_F_K} ${PARAM_F_FREQ} ${PARAM_F_POR} ${PARAM_ISO_PATH}
# ${LAUNCH} -n {$NCPU} -hostfile ${PBS_NODEFILE} --nooversubscribe ${WORKDIR}/${PROGNAME} ${PARAM_RESTART} ${SCRDIR} ${PARAM_PATH} ${PARAM_CRASHSTEP} ${PARAM_O} ${PARAM_W} ${PARAM_W_RHS} ${PARAM_H} ${PARAM_O_RHS} ${PARAM_TSW} ${PARAM_DIC} ${PARAM_SCOPE} ${PARAM_TRACE} ${PARAM_CH} ${PARAM_PAQ} ${PARAM_CH_RHS} ${PARAM_F_DX} ${PARAM_F_K} ${PARAM_F_FREQ} ${PARAM_F_POR}

# wait
# ssh $PBS_NODEFILE
# cd $SCRDIR
# mv * $PARAM_PATH_ALT
