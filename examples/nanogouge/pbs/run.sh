#!/bin/bash
# Same use as the centre's sample (cp /home/app/lammps/run.sh .): edit the lines marked EDIT, then   qsub run.sh
# One stage by hand.  (./go.sh submits the whole chain; this is for a single stage or a test.)
#PBS -N nanogouge
#PBS -q EDIT_QUEUE                       # EDIT: queue name (copy the "#PBS -q" line of the centre's run.sh; delete this line for the default queue)
#PBS -l select=4:ncpus=128:mpiprocs=128  # EDIT: nodes and cores per node (copy from the centre's run.sh); all cores of every node are used
#PBS -l walltime=24:00:00                # EDIT
#PBS -j oe
export STAGE=01                          # EDIT: 01 equilibration, 02 dehydration, 03 consolidation, 04 shear, 05 hold
export SIGMA=5.6                         # EDIT: normal stress in MPa (stages 03-05)
cd "${PBS_O_WORKDIR:-$PWD}"              # run folder: in.*, job.sh, config.env, gouge.data and gouge.settings are here (go.sh / submit_chain.sh make it: RUN_N12)
exec bash ./job.sh
