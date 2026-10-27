"""Run the reference MD of build_refs.py on another computer over SSH (soapff_gui.remote, as the GUI's remote runs): the inputs of
refs/<system>/ (made with "build_refs.py none") are copied to <remote_root>/library_water/<system>/, run there by a queue of
`jobs` LAMMPS processes at a time (one MPI process each), and copied back.  The password is taken from the environment variable
SOAPFF_PW (never written to a file; empty or unset: key login).

    python run_remote.py start user@host:port jobs system ...     (with the Python of soapff-gui)
    python run_remote.py status user@host:port
    python run_remote.py fetch user@host:port system ...           ->  refs/<system>/ (log.lammps, traj.lammpstrj, final.data ...)"""
import os
import shlex
import subprocess
import sys
import tempfile

from soapff_gui.remote import Remote

HERE = os.path.dirname(os.path.abspath(__file__))
mode, where = sys.argv[1], sys.argv[2]
user, _, hp = where.rpartition("@")
host, _, port = hp.partition(":")
r = Remote(host, user, int(port or 22), remote_root="soapff_runs/library_water", password=os.environ.get("SOAPFF_PW") or None)
r.check()
try:
    if mode == "start":
        jobs, names = int(sys.argv[3]), sys.argv[4:]
        for n in names:
            r.upload(os.path.join(HERE, "refs", n), n)
        lmp = r.lmp
        queue = ("#!/bin/sh\n# the reference MD of library_water: a queue of LAMMPS runs, `jobs` at a time\n"
                 f"cd \"$(dirname \"$0\")\"\nprintf '%s\\n' {' '.join(shlex.quote(n) for n in names)} | "
                 f"xargs -P {jobs} -I@ sh -c 'cd @ && OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl {lmp} -in in.lammps -log log.lammps "
                 "> screen.txt 2>&1; echo $? > exit_code'\n")
        with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as f:
            f.write(queue)
        subprocess.run([r.scp_cmd, "-P", str(r.port)] + r._opts() + [f.name, f"{r.target}:{r.remote_root}/queue.sh"],
                       env=r._env(), check=True, capture_output=True)
        os.remove(f.name)
        _, out, _ = r.run(f"cd {r.remote_root} || exit 1; nohup sh queue.sh > queue.log 2>&1 < /dev/null & echo $!")
        print("queue started, pid", out.strip())
    elif mode == "status":
        _, out, _ = r.run(f"cd {r.remote_root} && for d in */; do printf '%s ' $d; cat $d/exit_code 2>/dev/null || "
                          "(grep -E '^ +[0-9]+ +[-0-9.]' $d/log.lammps 2>/dev/null | tail -1 | awk '{print \"step\", $1}'); echo; done; uptime")
        print(out)
    elif mode == "fetch":
        for n in sys.argv[3:]:
            got = r.download_all(n, os.path.join(HERE, "refs", n))
            print(n, len(got), "files")
finally:
    r.close()
