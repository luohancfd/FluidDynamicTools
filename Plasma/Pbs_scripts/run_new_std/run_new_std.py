#!/usr/bin/env python3
# %%
import os
import argparse
import shlex
import time
import socket
from subprocess import PIPE, Popen

# %%
PBS_TEMPLATE = '''#PBS -N {jobname:s}
#PBS -l walltime={walltime:d}:00:00
#PBS -l nodes={node:d}:ppn=20
#PBS -q {queue:s}
#PBS -m bea

set -e
module load impi
module load comsol/5.4
module load matlab/R2018a

cd $PBS_O_WORKDIR

# ========== setup of comsol ========================
export COMSOL_LIC="{lic:s}"
export RECOVER_DIR=$(pwd)/recovery
export BACKUP_DIR=$(pwd)/backup
export TEMP_DIR=$(pwd)/temp
export NN=${{PBS_NUM_NODES}}
export NP=${{PBS_NUM_PPN}}
echo "Running on the following nodes:"
export hostfile="$PBS_JOBID.hostfile"
cat $PBS_NODEFILE|uniq >  $hostfile
mkdir -p ${{RECOVER_DIR}}
mkdir -p ${{TEMP_DIR}}
mkdir -p ${{BACKUP_DIR}}
export AUX_ARG="{auxargs:s}"

# =========== setup of the case to run ===============
export CASENAME="{casename:s}"  # name of the case
export INPUTFILE="${{CASENAME}}.mph"
export TL="{timeLength:d}"      # new study length
export TD="{timeDiff:d}"        # difference between output
export NEW_STD="{newStd:d}"

export OLD_STD="{oldStd:d}"      # OLDSTD
export OLD_OUTPUT="${{CASENAME}}_std${{OLD_STD}}.mph"
export OLD_LOG="${{CASENAME}}_std${{OLD_STD}}.mph.log"
export OLD_STATUS="${{CASENAME}}_std${{OLD_STD}}.mph.status"
export OLD_RECOVER="${{CASENAME}}_std${{OLD_STD}}.mph.recovery"

export CREATE_NEW_CASE="0"

if [ -f ${{OLD_OUTPUT}} ]; then
    if [ ! -z ${{OLD_STATUS}} ]; then
        if [ -f ${{OLD_STATUS}} ]; then
            if grep -q "Done" ${{OLD_STATUS}}; then
                CREATE_NEW_CASE="1"
            fi
        fi
    fi
fi

# =========== Select which case to run =====================
if [ ${{CREATE_NEW_CASE}} -eq "0" ]; then
    if [ -f ${{OLD_RECOVER}} ]; then
        echo "Recover old case ${{OLD_OUTPUT}}" >> log.txt
        export AUX_ARG="${{AUX_ARG}} -recover"
        export OLD_LOG="${{CASENAME}}_std${{OLD_STD}}_$(date +%s).mph.log"
    fi
    export COMSOL_ARG="-mpi intel -nn ${{NN}} -np ${{NP}} ${{COMSOL_LIC}} -f ${{hostfile}} -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -autosave on"
    comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{OLD_OUTPUT}} -study "std${{OLD_STD}}" -batchlog ${{OLD_LOG}} ${{AUX_ARG}}
elif [ ${{OLD_STD}} != ${{NEW_STD}} ]; then
    matlab -nodesktop -nosplash -r "setupCase('${{CASENAME}}',${{TL}},${{TD}},0,${{OLD_STD}},${{NEW_STD}});exit" > ${{CASENAME}}_std${{NEW_STD}}_matlab.output

    export NEW_STD_TEST=$(cat stdList.txt)
    export NEW_OUTPUT="${{CASENAME}}_std${{NEW_STD}}.mph"

    if [ -z ${{NEW_STD_TEST}} ]; then
        echo "stdList.txt is empty" >> log.txt
        exit 1
    fi

    if [ ${{NEW_STD_TEST}} -lt "0" ] || [ ${{NEW_STD}} != ${{NEW_STD_TEST}} ]; then
        echo "something wrong" >> log.txt
        echo "something wrong"
        exit 1
    fi

    if [ ! -f ${{NEW_OUTPUT}} ]; then
        echo "${{NEW_OUTPUT}} is not generated" >> log.txt
        exit 1
    else
        cp ${{INPUTFILE}} "backup/${{CASENAME}}_$(date +%s).mph"
        mv ${{NEW_OUTPUT}} ${{INPUTFILE}}
    fi

    export NEW_LOG="${{CASENAME}}_std${{NEW_STD}}.mph.log"
    export NEW_STATUS="${{CASENAME}}_std${{NEW_STD}}.mph.status"
    export NEW_RECOVER="${{CASENAME}}_std${{NEW_STD}}.mph.recovery"

    echo "Run new case ${{NEW_OUTPUT}}" >> log.txt
    export COMSOL_ARG="-mpi intel -nn ${{NN}} -np ${{NP}} ${{COMSOL_LIC}} -f ${{hostfile}} -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -autosave on"
    comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{NEW_OUTPUT}} -study "std${{NEW_STD}}" -batchlog ${{NEW_LOG}} ${{AUX_ARG}}
fi
'''

SLURM_TEMPLATE = '''#!/bin/bash
#SBATCH -J {jobname:s}
#SBATCH -t {walltime:d}:00:00
#SBATCH -N {node:d}
#SBATCH --ntasks-per-node=20
#SBATCH -A {queue:s}
#SBATCH --mail-type=ALL

set -e
set -x
module load impi
module load comsol/5.4
module load matlab/R2018a

cd ${{SLURM_SUBMIT_DIR}}

# ========== setup of comsol ========================
export COMSOL_LIC="{lic:s}"
export RECOVER_DIR=$(pwd)/recovery
export BACKUP_DIR=$(pwd)/backup
export TEMP_DIR=$(pwd)/temp
export NN=${{SLURM_JOB_NUM_NODES}}
export NP=20
echo "Running on the following nodes:"
export hostfile="${{SLURM_JOB_ID}}.$(hostname).hostfile"
scontrol show hostnames $SLURM_JOB_NODELIST >  $hostfile
mkdir -p ${{RECOVER_DIR}}
mkdir -p ${{TEMP_DIR}}
mkdir -p ${{BACKUP_DIR}}
export AUX_ARG="{auxargs:s}"

# =========== setup of the case to run ===============
export CASENAME="{casename:s}"  # name of the case
export INPUTFILE="${{CASENAME}}.mph"
export TL="{timeLength:d}"      # new study length
export TD="{timeDiff:d}"        # difference between output
export NEW_STD="{newStd:d}"

export OLD_STD="{oldStd:d}"      # OLDSTD
export OLD_OUTPUT="${{CASENAME}}_std${{OLD_STD}}.mph"
export OLD_LOG="${{CASENAME}}_std${{OLD_STD}}.mph.log"
export OLD_STATUS="${{CASENAME}}_std${{OLD_STD}}.mph.status"
export OLD_RECOVER="${{CASENAME}}_std${{OLD_STD}}.mph.recovery"

export CREATE_NEW_CASE="0"

if [ -f ${{OLD_OUTPUT}} ]; then
    if [ ! -z ${{OLD_STATUS}} ]; then
        if [ -f ${{OLD_STATUS}} ]; then
            if grep -q "Done" ${{OLD_STATUS}}; then
                CREATE_NEW_CASE="1"
            fi
        fi
    fi
fi

# =========== Select which case to run =====================
if [ ${{CREATE_NEW_CASE}} -eq "0" ]; then
    if [ -f ${{OLD_RECOVER}} ]; then
        echo "Recover old case ${{OLD_OUTPUT}}" >> log.txt
        export AUX_ARG="${{AUX_ARG}} -recover"
        export OLD_LOG="${{CASENAME}}_std${{OLD_STD}}_$(date +%s).mph.log"
    fi
    export COMSOL_ARG="-mpi intel -nn ${{NN}} -np ${{NP}} ${{COMSOL_LIC}} -f ${{hostfile}} -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -autosave on"
    comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{OLD_OUTPUT}} -study "std${{OLD_STD}}" -batchlog ${{OLD_LOG}} ${{AUX_ARG}}
elif [ ${{OLD_STD}} != ${{NEW_STD}} ]; then
    matlab -nodesktop -nosplash -r "setupCase('${{CASENAME}}',${{TL}},${{TD}},0,${{OLD_STD}},${{NEW_STD}});exit" > ${{CASENAME}}_std${{NEW_STD}}_matlab.output

    export NEW_STD_TEST=$(cat stdList.txt)
    export NEW_OUTPUT="${{CASENAME}}_std${{NEW_STD}}.mph"

    if [ -z ${{NEW_STD_TEST}} ]; then
        echo "stdList.txt is empty" >> log.txt
        exit 1
    fi

    if [ ${{NEW_STD_TEST}} -lt "0" ] || [ ${{NEW_STD}} != ${{NEW_STD_TEST}} ]; then
        echo "something wrong" >> log.txt
        echo "something wrong"
        exit 1
    fi

    if [ ! -f ${{NEW_OUTPUT}} ]; then
        echo "${{NEW_OUTPUT}} is not generated" >> log.txt
        exit 1
    else
        cp ${{INPUTFILE}} "backup/${{CASENAME}}_$(date +%s).mph"
        mv ${{NEW_OUTPUT}} ${{INPUTFILE}}
    fi

    export NEW_LOG="${{CASENAME}}_std${{NEW_STD}}.mph.log"
    export NEW_STATUS="${{CASENAME}}_std${{NEW_STD}}.mph.status"
    export NEW_RECOVER="${{CASENAME}}_std${{NEW_STD}}.mph.recovery"

    echo "Run new case ${{NEW_OUTPUT}}" >> log.txt
    export COMSOL_ARG="-mpi intel -nn ${{NN}} -np ${{NP}} ${{COMSOL_LIC}} -f ${{hostfile}} -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -autosave on"
    comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{NEW_OUTPUT}} -study "std${{NEW_STD}}" -batchlog ${{NEW_LOG}} ${{AUX_ARG}}
fi
'''

if __name__ == "__main__":
    # %%
    parser = argparse.ArgumentParser(
        description="Continuously generate new COMSOL study")
    parser.add_argument('CASE', help='full path of the case')
    parser.add_argument('STD', type=int, help='ID of the first study to run')
    parser.add_argument('NRUN', type=int, help='Number of runs')
    parser.add_argument('--tag', default='', help='Tag infront of job name')
    parser.add_argument('--depid', type=int, default=0,
                        help='Dependent job of the first job (default: None)')
    parser.add_argument('-n', '--node', type=int, default=1,
                        help='Number of nodes to use (default: %(default)s)')
    parser.add_argument('-q', '--queue', default="standby",
                        help='PBS queue (default: %(default)s)')
    parser.add_argument('-w', '--walltime', type=int, default=4,
                        help='PBS wall time in hours (default: %(default)s)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Dry run with fake job id')
    parser.add_argument('--sleep', type=int, default=0,
                        help='Sleep for number of hours at the end of script (default: %(default)s)')
    parser.add_argument('--tl', type=int, default=500,
                        help='Number of periods (default: %(default)s)')
    parser.add_argument('--td', type=int, default=50,
                        help='Output interval (default: %(default)s)')
    parser.add_argument('--remote', default='slurm',
                        help='Submit the job to PBS or SLURM system (default: %(default)s)')


    hostname = socket.gethostname()
    COMSOL_LIC = None
    parser.add_argument('--lic', default=COMSOL_LIC,
                        help='location of license file (default: %(default)s)')


    ARGS = parser.parse_args()

    if ARGS.remote.lower() == 'slurm':
        TEMPLATE = SLURM_TEMPLATE
        submitter = 'sbatch'
        depArgs = '--dependency=afterany:'

        def getJobID(o):
            w = [i for i in o.split() if o.strip()]
            return int(w[-1])
    else:
        TEMPLATE = PBS_TEMPLATE
        submitter = 'qsub'
        depArgs = '-W depend=afterany:'

        def getJobID(o):
            return int(o[0:8])

    cwd = os.getcwd()
    abspath = os.path.abspath(ARGS.CASE)
    if not os.path.isfile(abspath):
        raise FileNotFoundError('{:s} doesn''t exists'.format(ARGS.CASE))
    folder = os.path.dirname(abspath)
    casename = os.path.basename(abspath).split('.')[0]

    # create pbs files
    os.chdir(folder)
    for i, istd in enumerate(range(ARGS.STD, ARGS.STD + ARGS.NRUN)):
        pbs_file = '{:s}_std{:d}.pbs'.format(casename, istd)
        if ARGS.tag:
            jobname = '{:s}_{:s}_std{:d}'.format(
                ARGS.tag, casename, istd)
        else:
            jobname = '{:s}_std{:d}'.format(casename, istd)
        with open(pbs_file, 'w', encoding='utf8') as f:
            auxargs = '-continue'
            print(istd)
            if istd == 1:
                oldStd = 1
                newStd = 1
            else:
                oldStd = istd - 1
                newStd = istd
            f.write(TEMPLATE.format(jobname=jobname,
                                    walltime=ARGS.walltime,
                                    node=ARGS.node,
                                    queue=ARGS.queue,
                                    casename=casename,
                                    timeLength=ARGS.tl,
                                    timeDiff=ARGS.td,
                                    newStd=newStd,
                                    oldStd=oldStd,
                                    lic="-c "+ARGS.lic if ARGS.lic else "",
                                    auxargs=auxargs))
            if ARGS.sleep > 0:
                f.write("\nsleep {:d}h\n".format(ARGS.sleep))

    # submit job
    if ARGS.depid != 0:
        prev_job_id = ARGS.depid
    else:
        prev_job_id = 0

    for i, istd in enumerate(range(ARGS.STD, ARGS.STD + ARGS.NRUN)):
        pbs_file = '{:s}_std{:d}.pbs'.format(casename, istd)
        if i == 0:
            if ARGS.depid == 0:
                command = "{:s} {:s}".format(submitter, pbs_file)
            else:
                command = "{:s} {:s}{:d} {:s}".format(
                        submitter, depArgs, prev_job_id, pbs_file)
        else:
            command = "{:s} {:s}{:d} {:s}".format(
                        submitter, depArgs, prev_job_id, pbs_file)

        if ARGS.dry_run:
            prev_job_id = istd
            o = "Dry run {:s}".format(pbs_file)
        else:
            c = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE)
            o, e = c.communicate()
            o = o.decode("utf-8")
            prev_job_id = getJobID(o)

        print('{:2d} PBS file: {:s}'.format(istd, pbs_file))
        print('{:2d} Command: {:s}'.format(istd, command))
        print('{:2d} Output:  {:s}'.format(istd, o))
        print("{:2d} JOB_ID:{:d}".format(istd, prev_job_id))

        print("")
        time.sleep(2)

    os.chdir(cwd)
