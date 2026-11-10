#!/bin/bash
# Make the flat, self-contained pack  nanogouge_pack.tar.gz  (no git, no soapff-gui needed on the supercomputer).
set -eu
cd "$(dirname "$0")"
P=nanogouge_pack; rm -rf $P; mkdir $P
cp ../montmorillonite_edge/out_clay_montmorillonite/system.data $P/particle.data
cp ../montmorillonite_edge/out_clay_montmorillonite/system.settings $P/particle.settings
cp build_gouge.py go.sh lammps/in.* pbs/job.sh pbs/status.sh pbs/autoconfig.sh pbs/run.sh pbs/submit_chain.sh pbs/config.env analysis/*.py $P/
cp README.md $P/README.md
tar czf $P.tar.gz $P && rm -rf $P
ls -la $P.tar.gz
