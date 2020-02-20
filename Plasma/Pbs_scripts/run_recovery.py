#!/usr/bin/env python3
# %%
import os
import readline
import glob
import argparse
import shlex
import stat
import time
import socket
from subprocess import PIPE, Popen

# %%
PBS_TEMPLATE = '''#PBS -N {jobname:s}
#PBS -l walltime={walltime:d}:00:00
#PBS -l nodes={node:d}:ppn=20
#PBS -q {queue:s}
#PBS -m bea

cd $PBS_O_WORKDIR

export CASENAME="{casename:s}"
export INPUTFILE="./${{CASENAME}}.mph"
export STUDY=std{studyid:d}
export STDFILE="${{CASENAME}}_${{STUDY}}"
export OUTPUTFILE="./${{STDFILE}}.mph"
export STATUS_FILE="${{STDFILE}}.mph.status"
export RECOVER_FILE="${{STDFILE}}.mph.recovery"

export NRUN={runid:d}
export PREV_RUN=`echo ${{NRUN}}-1 | bc -l`
export LOGFILE="./${{STDFILE}}_${{NRUN}}.log"

if [ -f $OUTPUTFILE ]; then
    mkdir -p backup
    cp $OUTPUTFILE "backup/${{STDFILE}}_run_${{PREV_RUN}}.mph"
fi

export AUX_ARG="{auxargs:s}"
# Check if previous run is finished
if [ ! -z ${{STATUS_FILE}} ]; then
if [ -f ${{STATUS_FILE}} ]; then
    if grep -q "Done" ${{STATUS_FILE}}; then
        echo "Job was finished, quit" | mail -s "${{CASENAME}}" luo160@purdue.edu
        exit $?
        else
        echo "Recover job"
        export AUX_ARG="${{AUX_ARG}} -recover"
    fi
else
    echo "Normal run"
fi
fi

export ALIVE_TIME=60  # time between saving, every 5 mintues
export NN=${{PBS_NUM_NODES}}
export NP=${{PBS_NUM_PPN}}
export COMSOL_LIC="{lic:s}"
export RECOVER_DIR=$(pwd)/recovery
export TEMP_DIR=$(pwd)/temp


module load impi
module load comsol/5.4
echo "Running on the following nodes:"
export hostfile="$PBS_JOBID.hostfile"
cat $PBS_NODEFILE|uniq >  $hostfile
mkdir -p ${{RECOVER_DIR}}
mkdir -p ${{TEMP_DIR}}

export COMSOL_ARG="-mpi intel -nn $NN -np $NP $COMSOL_LIC -f $hostfile -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -alivetime ${{ALIVE_TIME}} -autosave on"

comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{OUTPUTFILE}} -study ${{STUDY}} -batchlog ${{LOGFILE}} ${{AUX_ARG}}
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
cd ${{SLURM_SUBMIT_DIR}}

export CASENAME="{casename:s}"
export INPUTFILE="./${{CASENAME}}.mph"
export STUDY=std{studyid:d}
export STDFILE="${{CASENAME}}_${{STUDY}}"
export OUTPUTFILE="./${{STDFILE}}.mph"
export STATUS_FILE="${{STDFILE}}.mph.status"
export RECOVER_FILE="${{STDFILE}}.mph.recovery"

export NRUN={runid:d}
export PREV_RUN=`echo ${{NRUN}}-1 | bc -l`
export LOGFILE="./${{STDFILE}}_${{NRUN}}.log"

if [ -f $OUTPUTFILE ]; then
    mkdir -p backup
    cp $OUTPUTFILE "backup/${{STDFILE}}_run_${{PREV_RUN}}.mph"
fi

export AUX_ARG="{auxargs:s}"
# Check if previous run is finished
if [ ! -z ${{STATUS_FILE}} ]; then
if [ -f ${{STATUS_FILE}} ]; then
    if grep -q "Done" ${{STATUS_FILE}}; then
        echo "Job was finished, quit" | mail -s "${{CASENAME}}" luo160@purdue.edu
        exit $?
        else
        echo "Recover job"
        export AUX_ARG="${{AUX_ARG}} -recover"
    fi
else
    echo "Normal run"
fi
fi

export ALIVE_TIME=60  # time between saving, every 5 mintues
export NN=${{SLURM_JOB_NUM_NODES}}
export NP=20
export COMSOL_LIC="{lic:s}"
export RECOVER_DIR=$(pwd)/recovery
export TEMP_DIR=$(pwd)/temp


module load impi
module load comsol/5.4
echo "Running on the following nodes:"
export hostfile="${{SLURM_JOB_ID}}.$(hostname).hostfile"
scontrol show hostnames $SLURM_JOB_NODELIST >  $hostfile
mkdir -p ${{RECOVER_DIR}}
mkdir -p ${{TEMP_DIR}}

export COMSOL_ARG="-mpi intel -nn $NN -np $NP $COMSOL_LIC -f $hostfile -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -alivetime ${{ALIVE_TIME}} -autosave on"

comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{OUTPUTFILE}} -study ${{STUDY}} -batchlog ${{LOGFILE}} ${{AUX_ARG}}
'''

# ====================================================================================
ECN_TEMPLATE = '''#!/bin/bash
set -e
set -x
cd {workdir:s}

export CASENAME="{casename:s}"
export INPUTFILE="./${{CASENAME}}.mph"
export STUDY=std{studyid:d}
export STDFILE="${{CASENAME}}_${{STUDY}}"
export OUTPUTFILE="./${{STDFILE}}.mph"
export STATUS_FILE="${{STDFILE}}.mph.status"
export RECOVER_FILE="${{STDFILE}}.mph.recovery"
export LOGFILE="./${{STDFILE}}.log"

if [ -f $OUTPUTFILE ]; then
    exit 1
fi

export AUX_ARG="{auxargs:s}"
# Check if previous run is finished
if [ ! -z ${{STATUS_FILE}} ]; then
if [ -f ${{STATUS_FILE}} ]; then
    if grep -q "Done" ${{STATUS_FILE}}; then
        echo "Job was finished, quit" | mail -s "${{CASENAME}}" luo160@purdue.edu
        exit $?
        else
        echo "Recover job"
        export AUX_ARG="${{AUX_ARG}} -recover"
    fi
else
    echo "Normal run"
fi
fi

export ALIVE_TIME=60  # time between saving, every 5 mintues
export COMSOL_LIC="{lic:s}"
export RECOVER_DIR=$(pwd)/recovery
export TEMP_DIR=$(pwd)/temp

if [[ $(hostname) == *"ecn"* ]]; then
    module load impi
fi
mkdir -p ${{RECOVER_DIR}}
mkdir -p ${{TEMP_DIR}}

export COMSOL_ARG="-mpi intel $COMSOL_LIC -recoverydir ${{RECOVER_DIR}} -tmpdir ${{TEMP_DIR}} -alivetime ${{ALIVE_TIME}} -autosave on"

nohup bash -c "comsol batch ${{COMSOL_ARG}} -inputfile ${{INPUTFILE}} -outputfile ${{OUTPUTFILE}} -study ${{STUDY}} -batchlog ${{LOGFILE}} ${{AUX_ARG}}" > /dev/null 2>&1 &

echo $! > pid.log
'''

# %%


class tabCompleter(object):
    """
    A tab completer that can either complete from
    the filesystem or from a list.

    Partially taken from:
    http://stackoverflow.com/questions/5637124/tab-completion-in-pythons-raw-input
    https://gist.github.com/iamatypeofwalrus/5637895
    """

    def pathCompleter(self, text, state):
        """
        This is the tab completer for systems paths.
        Only tested on *nix systems
        """

        # line = readline.get_line_buffer().split()

        if '~' in text:
            text = os.path.expanduser('~')

        # autocomplete directories with having a trailing slash
        if os.path.isdir(text):
            text += '/'

        return [x for x in glob.glob(text + '*')][state]

    def createListCompleter(self, ll):
        """
        This is a closure that creates a method that autocompletes from
        the given list.

        Since the autocomplete function can't be given a list to complete from
        a closure is used to create the listCompleter function with a list to complete
        from.
        """
        def listCompleter(text, state):
            line = readline.get_line_buffer()

            if not line:
                return [c + " " for c in ll][state]

            else:
                return [c + " " for c in ll if c.startswith(line)][state]

        self.listCompleter = listCompleter


if __name__ == "__main__":
    # %%
    parser = argparse.ArgumentParser(
        description="Generate multiple comsol run")
    parser.add_argument('CASE', help='full path of the case')
    parser.add_argument('STD', type=int, help='ID of the study to run')
    parser.add_argument('NRUN', type=int, help='Number of runs')
    parser.add_argument('--tag', default='', help='Tag in front of job name')
    parser.add_argument('--id', type=int, default=0,
                        help='Run ID of the first job (default: %(default)s) ')
    parser.add_argument('--depid', type=int, default=0,
                        help='Dependent job of the first job (default: None)')
    parser.add_argument('-n', '--node', type=int, default=2,
                        help='Number of nodes to use (default: %(default)s)')
    parser.add_argument('-q', '--queue', default="standby",
                        help='PBS/SLURM queue (default: %(default)s)')
    parser.add_argument('-w', '--walltime', type=int, default=4,
                        help='PBS wall time in hours (default: %(default)s)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Dry run with fake job id')
    parser.add_argument('--sleep', type=int, default=0,
                        help='Sleep for number of hours at the end of script (default: %(default)s)')
    parser.add_argument('--local', action='store_true',
                        help='run the job on local system with nohup')
    parser.add_argument('--remote', default='slurm',
                        help='Submit the job to PBS or SLURM system (default: %(default)s)')

    hostname = socket.gethostname()
    COMSOL_LIC = None

    parser.add_argument('--lic', default=COMSOL_LIC,
                        help='location of license file (default: %(default)s)')

    ARGS = parser.parse_args()

    cwd = os.getcwd()
    abspath = os.path.abspath(ARGS.CASE)
    if not os.path.isfile(abspath):
        raise FileNotFoundError('{:s} doesn''t exists'.format(ARGS.CASE))
    folder = os.path.dirname(abspath)
    casename = os.path.basename(abspath).split('.')[0]

    # check if this is local run
    if ARGS.local:
        TEMPLATE = ECN_TEMPLATE
        ARGS.NRUN = 1
        submitter = 'bash'
        print("Local run only allow nrun=1")
    else:
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


    # create pbs files
    os.chdir(folder)
    for i, irun in enumerate(range(ARGS.id, ARGS.id + ARGS.NRUN)):
        if ARGS.local:
            pbs_file = '{:s}_std{:d}_{:d}.sh'.format(casename, ARGS.STD, irun)
        else:
            pbs_file = '{:s}_std{:d}_{:d}.pbs'.format(casename, ARGS.STD, irun)
        if ARGS.tag:
            jobname = '{:s}_{:s}_std{:d}_{:d}'.format(
                ARGS.tag, casename, ARGS.STD, irun)
        else:
            jobname = '{:s}_std{:d}_{:d}'.format(casename, ARGS.STD, irun)
        with open(pbs_file, 'w', encoding='utf8') as f:
            auxargs = '-continue'
            if ARGS.local:
                f.write(TEMPLATE.format(workdir=folder,
                                        casename=casename,
                                        studyid=ARGS.STD,
                                        auxargs=auxargs,
                                        lic="-c "+ARGS.lic if ARGS.lic else ""
                                        ))
            else:
                f.write(TEMPLATE.format(jobname=jobname,
                                        node=ARGS.node,
                                        casename=casename,
                                        studyid=ARGS.STD,
                                        runid=irun,
                                        auxargs=auxargs,
                                        queue=ARGS.queue,
                                        lic="-c "+ARGS.lic if ARGS.lic else "",
                                        walltime=ARGS.walltime))
                if ARGS.sleep > 0:
                    f.write("\nsleep {:d}h\n".format(ARGS.sleep))

        # Make it executable
        st = os.stat(pbs_file)
        os.chmod(pbs_file, st.st_mode | stat.S_IEXEC)

    # submit job
    if ARGS.local:
        for i, irun in enumerate(range(ARGS.id, ARGS.id + ARGS.NRUN)):
            pbs_file = '{:s}_std{:d}_{:d}.sh'.format(casename, ARGS.STD, irun)
            command = "./{:s}".format(pbs_file)
            if ARGS.dry_run:
                print(command)
            else:
                c = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE)
                o, e = c.communicate()
                o = o.decode("utf-8")

            print('{:2d} SH file: {:s}'.format(irun, pbs_file))
            print('{:2d} Command: {:s}'.format(irun, command))
    else:
        if ARGS.depid != 0:
            prev_job_id = ARGS.depid
        else:
            prev_job_id = 0

        for i, irun in enumerate(range(ARGS.id, ARGS.id + ARGS.NRUN)):
            pbs_file = '{:s}_std{:d}_{:d}.pbs'.format(casename, ARGS.STD, irun)
            if i == 0:
                if ARGS.depid == 0:
                    command = "{:s} {:s}".format(submitter, pbs_file)
                else:
                    command = "{:s} {:s}{:d} {:s}".format(
                        submitter, depArgs, prev_job_id, pbs_file)
            else:
                command = "{:s} {:s}{:d} {:s}".format(submitter, depArgs,
                                                      prev_job_id, pbs_file)

            if ARGS.dry_run:
                prev_job_id = irun
                o = "Dry run {:s}".format(pbs_file)
            else:
                c = Popen(shlex.split(command), stdout=PIPE, stderr=PIPE)
                o, e = c.communicate()
                o = o.decode("utf-8")
                prev_job_id = getJobID(o)

            print('{:2d} PBS file: {:s}'.format(irun, pbs_file))
            print('{:2d} Command: {:s}'.format(irun, command))
            print('{:2d} Output:  {:s}'.format(irun, o))
            print("{:2d} JOB_ID:{:d}".format(irun, prev_job_id))

            print("")
            time.sleep(2)

    os.chdir(cwd)
