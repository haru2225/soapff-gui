"""Running LAMMPS on another computer over SSH (the system ssh / scp; no extra Python package).

The force field is written here (soapff-cui), the three input files are copied to <remote_root>/<name>/ on the remote
computer, LAMMPS is started there with nohup (it keeps running if the connection drops), its screen output is followed
with tail, and the results are copied back.  One SSH connection is shared by all commands (ControlMaster), so a password
is asked once; it is handed to ssh through SSH_ASKPASS and never written to disk (an empty password: key login)."""
import os
import shlex
import stat
import subprocess
import tempfile

INPUTS = ("system.data", "system.settings", "in.lammps")
OPTIONAL_INPUTS = ("water.mol",)        # sent when present (the water molecule of the swelling curve)
RESULTS = ("log.lammps", "screen.txt", "exit_code", "traj.lammpstrj", "final.data", "msd_*.dat", "swelling.dat",
           "swelling_states.lammpstrj")


class RemoteError(RuntimeError):
    pass


def _family(text):
    t = (text or "").lower()
    if "stubs" in t:
        return "none"
    if "open mpi" in t or "openrte" in t or "ompi" in t:
        return "openmpi"
    if "mpich" in t or "hydra" in t:
        return "mpich"
    if "intel" in t:
        return "intel"
    return ""


def mpi_warning(info, mpirun):
    """a Japanese warning when lmp cannot run in parallel with that mpirun ('' when fine or unknown)."""
    lm = info.get("lmp_mpi", "")
    if _family(lm) == "none":
        return ("この lmp は MPI なしで作られています（lmp -h：「" + lm.strip() + "」）。mpirun で何並列にしても、"
                "1 並列の LAMMPS が並列数の分だけ別々に動きます。並列数を 1 にするか、MPI 版の lmp を指定してください")
    run = next((v for k, v in info.items() if k.startswith("mpirun_version[") and k[15:-1] == mpirun), "")
    a, b = _family(lm), _family(run)
    if a and b and a != b:
        return (f"lmp は {lm.strip()} で作られていますが、mpirun（{mpirun}）は別の MPI（{run.strip()}）です。"
                "このままだと 1 並列の LAMMPS が並列数の分だけ別々に動きます。lmp と同じ MPI の mpirun を指定してください")
    return ""


def explain(stderr):
    """a Japanese hint for the common ssh failures."""
    s = stderr or ""
    if "Permission denied" in s:
        return ("ログインできませんでした（認証に失敗）。パスワードを入れるか、「鍵を登録」で公開鍵を相手に登録してください。"
                "\n" + s.strip())
    if "Connection refused" in s:
        return "接続を拒否されました。ホスト名・ポート番号と、相手の SSH サーバーが動いているかを確認してください。\n" + s.strip()
    if "timed out" in s or "Operation timed out" in s:
        return "接続できませんでした（時間切れ）。ネットワーク（VPN、Tailscale など）とアドレスを確認してください。\n" + s.strip()
    if "Could not resolve hostname" in s:
        return "ホスト名が見つかりません。\n" + s.strip()
    if "Host key verification failed" in s or "REMOTE HOST IDENTIFICATION HAS CHANGED" in s:
        return "ホストの鍵が以前と違います（~/.ssh/known_hosts を確認してください）。\n" + s.strip()
    return s.strip()


class Remote:
    """host, user, port; key (optional private key file); remote_root (folder for the runs); lmp, mpirun (paths or
    names on the remote computer); setup (shell commands run first, e.g. 'module load lammps' or
    'source ~/.bashrc').  password: kept in memory only."""

    def __init__(self, host, user="", port=22, key="", remote_root="soapff_runs", lmp="lmp", mpirun="mpirun", setup="",
                 password=None, ssh="ssh", scp="scp"):
        if not host:
            raise RemoteError("ホストが指定されていません")
        self.host, self.user, self.port = host.strip(), user.strip(), int(port or 22)
        self.key, self.remote_root = os.path.expanduser(key.strip()) if key else "", (remote_root or "soapff_runs").strip()
        self.lmp, self.mpirun, self.setup = lmp.strip() or "lmp", mpirun.strip() or "mpirun", setup.strip()
        self.password = password or None
        self.ssh_cmd, self.scp_cmd = ssh, scp
        d = os.path.join(os.path.expanduser("~"), ".ssh")
        # one shared connection per Remote (per job): with one per host, closing it (end of a job, a connection test)
        # cut the connections of the other windows or jobs to the same computer
        Remote._count = getattr(Remote, "_count", 0) + 1
        self.control = os.path.join(d if os.path.isdir(d) else tempfile.gettempdir(), f"sf{os.getpid()}-{Remote._count}-%C")

    @property
    def target(self):
        return f"{self.user}@{self.host}" if self.user else self.host

    def label(self):
        return f"{self.target}:{self.port}"

    # ----------------------------------------------------------------------------------------------- ssh plumbing
    def _opts(self):
        o = ["-o", "ConnectTimeout=15", "-o", "StrictHostKeyChecking=accept-new", "-o", "ServerAliveInterval=30",
             "-o", "ControlMaster=auto", "-o", f"ControlPath={self.control}", "-o", "ControlPersist=30m"]
        o += ["-o", "BatchMode=no" if self.password else "BatchMode=yes"]
        if self.key:
            o += ["-i", self.key]
        return o

    def _env(self):
        env = dict(os.environ)
        if self.password:
            if not getattr(self, "_helper", None):
                fd, self._helper = tempfile.mkstemp(prefix="soapff-askpass-", suffix=".sh")
                with os.fdopen(fd, "w") as f:
                    f.write('#!/bin/sh\nprintf "%s\\n" "$SOAPFF_SSH_SECRET"\n')   # the password stays in the environment
                os.chmod(self._helper, stat.S_IRWXU)
            env.update(SSH_ASKPASS=self._helper, SSH_ASKPASS_REQUIRE="force", DISPLAY=env.get("DISPLAY", ":0"),
                       SOAPFF_SSH_SECRET=self.password)
        return env

    def _ssh_args(self, command, tty=False):
        return [self.ssh_cmd, "-p", str(self.port)] + self._opts() + (["-tt"] if tty else []) + [self.target, command]

    def run(self, command, timeout=60, check=True):
        """run a shell command on the remote computer; returns (returncode, stdout, stderr)."""
        p = subprocess.run(self._ssh_args(command), env=self._env(), stdin=subprocess.DEVNULL, capture_output=True,
                           text=True, timeout=timeout)
        if check and p.returncode != 0:
            raise RemoteError(explain(p.stderr) or f"リモートのコマンドが失敗しました（{p.returncode}）：{command}")
        return p.returncode, p.stdout, p.stderr

    def _shell(self, body):
        """a command line for the remote shell: the setup commands first (non-interactive shells skip ~/.bashrc)."""
        return f"{self.setup}; {body}" if self.setup else body

    # ----------------------------------------------------------------------------------------------- operations
    def check(self):
        """log in and report the remote computer: host name, OS, CPUs, where lmp and mpirun are, which MPI lmp was built
        with ('MPI v...' of lmp -h; 'STUBS' = no MPI) and which MPI the mpirun belongs to.  When an mpirun sits next to
        lmp (same installation, same MPI), it is used instead of the one on PATH (self.mpirun_used)."""
        q = shlex.quote
        body = ("echo host=$(hostname); echo os=$(uname -sm); "
                "echo cpus=$( (nproc || sysctl -n hw.ncpu) 2>/dev/null); "
                # physical cores (Open MPI allows that many processes by default; hyperthreads hardly speed up LAMMPS)
                "C=$(lscpu -p=core,socket 2>/dev/null | grep -v '^#' | sort -u | wc -l | tr -d ' '); "
                "[ \"${C:-0}\" -gt 0 ] 2>/dev/null || C=$(sysctl -n hw.physicalcpu 2>/dev/null); echo cores=$C; "
                f"L=$(command -v {q(self.lmp)}); echo lmp=${{L:-NOTFOUND}}; "
                f"M=$(command -v {q(self.mpirun)}); echo mpirun=${{M:-NOTFOUND}}; "
                'if [ -n "$L" ] && [ -x "$(dirname "$L")/mpirun" ]; then echo mpirun_next=$(dirname "$L")/mpirun; fi; '
                'if [ -n "$L" ]; then echo lmp_mpi=$("$L" -h 2>/dev/null | grep -m1 "^MPI v" | tr -s " "); fi; '
                'for X in "$M" "$(dirname "${L:-/x}")/mpirun"; do [ -x "$X" ] && echo "mpirun_version[$X]=$("$X" --version 2>&1 | '
                'grep -m1 -i -E "open mpi|openrte|hydra|mpich|intel" | tr -s " ")"; done; '
                f"echo root=$(mkdir -p {self.remote_root} && cd {self.remote_root} && pwd)")
        _, out, _ = self.run(self._shell(body), timeout=40)
        info = dict(l.split("=", 1) for l in out.splitlines() if "=" in l)
        self.mpirun_used = self.mpirun
        nxt = info.get("mpirun_next")
        if nxt and "/" not in self.mpirun:        # 'mpirun' (from PATH) given: prefer the one of the lmp installation
            self.mpirun_used = nxt
        info["mpirun_used"] = self.mpirun_used
        info["mpi_warning"] = mpi_warning(info, self.mpirun_used)
        run = next((v for k, v in info.items() if k.startswith("mpirun_version[") and k[15:-1] == self.mpirun_used), "")
        self.openmpi = _family(run or info.get("lmp_mpi", "")) == "openmpi"
        self.cores = int(info["cores"]) if str(info.get("cores", "")).isdigit() else 0
        return info

    def run_dir(self, name):
        return f"{self.remote_root}/{name}"

    def upload(self, local_dir, name):
        """copy the input files of local_dir to <remote_root>/<name>/ (created, old results removed)."""
        rdir = self.run_dir(name)
        q = shlex.quote(rdir)
        self.run(f"mkdir -p {q} && cd {q} && rm -f {' '.join(RESULTS)}")
        files = [os.path.join(local_dir, f) for f in INPUTS]
        files += [os.path.join(local_dir, f) for f in OPTIONAL_INPUTS if os.path.isfile(os.path.join(local_dir, f))]
        missing = [f for f in files if not os.path.isfile(f)]
        if missing:
            raise RemoteError(f"入力ファイルがありません：{missing}")
        p = subprocess.run([self.scp_cmd, "-P", str(self.port)] + self._opts() + files + [f"{self.target}:{rdir}/"],
                           env=self._env(), stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=600)
        if p.returncode != 0:
            raise RemoteError(explain(p.stderr) or "入力ファイルを送れませんでした")
        return rdir

    def start(self, name, nproc=1):
        """start LAMMPS with nohup in the run folder; returns the process id of the job (a shell that waits for LAMMPS
        and writes its exit code to exit_code)."""
        rdir = shlex.quote(self.run_dir(name))
        mpirun = getattr(self, "mpirun_used", None) or self.mpirun
        opts = ""
        if getattr(self, "openmpi", False) and int(nproc) > 1:
            # several jobs at once: Open MPI binds each job to the same first cores ('--bind-to none' leaves the placement
            # to the OS); more processes than physical cores need the hardware threads counted as slots
            opts = " --bind-to none" + (" --use-hwthread-cpus" if getattr(self, "cores", 0) and int(nproc) > self.cores else "")
        lmp = f"{mpirun}{opts} -np {int(nproc)} {self.lmp}" if int(nproc) > 1 else self.lmp
        # HWLOC_COMPONENTS=-gl: Open MPI's hwloc would ask the X server of a desktop machine ('Authorization required')
        job = (f"OMP_NUM_THREADS=1 HWLOC_COMPONENTS=-gl {lmp} -in in.lammps -log log.lammps > screen.txt 2>&1; "
               "echo $? > exit_code")
        # only the nohup command goes to the background ('cd ... && nohup ... &' would background the whole list and keep
        # the ssh output open until LAMMPS ends)
        body = f"cd {rdir} || exit 1; nohup sh -c {shlex.quote(self._shell(job))} > /dev/null 2>&1 < /dev/null & echo $!"
        _, out, _ = self.run(body)
        pid = out.strip().splitlines()[-1].strip() if out.strip() else ""
        if not pid.isdigit():
            raise RemoteError(f"LAMMPS を起動できませんでした：{out}")
        return int(pid)

    def follow(self, name, pid):
        """yields the screen output of the remote LAMMPS line by line until the job ends (portable sh: tail -f in the
        background, polled with kill -0)."""
        rdir = shlex.quote(self.run_dir(name))
        body = (f"cd {rdir} && touch screen.txt && (tail -n +1 -f screen.txt & T=$!; "
                f"while kill -0 {int(pid)} 2>/dev/null; do sleep 2; done; sleep 1; kill $T)")
        self.proc = subprocess.Popen(self._ssh_args(body), env=self._env(), stdin=subprocess.DEVNULL,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
        for line in self.proc.stdout:
            yield line.rstrip("\n")
        self.proc.wait()

    def exit_code(self, name):
        rc, out, _ = self.run(f"cat {shlex.quote(self.run_dir(name))}/exit_code 2>/dev/null", check=False)
        try:
            return int(out.strip())
        except ValueError:
            return None

    def stop(self, pid):
        """stop the remote job: the whole process tree under the job's shell (mpirun, its proxies and the LAMMPS
        processes), deepest first, then the shell."""
        tree = ("k() { for c in $(pgrep -P $1 2>/dev/null); do k $c; done; kill -TERM $1 2>/dev/null; }; "
                f"k {int(pid)}; true")
        self.run(tree, check=False)

    def list_dir(self, path=None):
        """the folders and files in a remote folder (default: the folder of the runs): [(name, is_dir)]."""
        path = path or self.remote_root
        _, out, _ = self.run(f"cd {shlex.quote(path)} && for f in * ; do [ -e \"$f\" ] || continue; "
                             f"if [ -d \"$f\" ]; then echo \"D $f\"; else echo \"F $f\"; fi; done", timeout=40)
        return [(l[2:], l[0] == "D") for l in out.splitlines() if len(l) > 2]

    def fetch(self, remote_path, local_parent):
        """copy a remote folder (or file) with everything in it into local_parent/ (scp -r, like cp -r)."""
        os.makedirs(local_parent, exist_ok=True)
        p = subprocess.run([self.scp_cmd, "-r", "-P", str(self.port)] + self._opts() + [f"{self.target}:{remote_path}",
                           local_parent.rstrip("/") + "/"], env=self._env(), stdin=subprocess.DEVNULL, capture_output=True,
                           text=True, timeout=24 * 3600)
        if p.returncode != 0:
            raise RemoteError(explain(p.stderr) or f"受け取れませんでした：{remote_path}")
        return os.path.join(local_parent, os.path.basename(remote_path.rstrip("/")))

    def download_all(self, name, local_dir):
        """everything in the run folder (cp -r) into local_dir: copied into a temporary folder next to it, then moved in
        (files of the same name are replaced)."""
        import shutil
        parent = os.path.dirname(os.path.abspath(local_dir))
        tmp = tempfile.mkdtemp(prefix=".soapff-fetch-", dir=parent)
        try:
            got = self.fetch(self.run_dir(name), tmp)
            names = []
            for root, dirs, files in os.walk(got):
                rel = os.path.relpath(root, got)
                dest = os.path.join(local_dir, rel) if rel != "." else local_dir
                os.makedirs(dest, exist_ok=True)
                for f in files:
                    shutil.move(os.path.join(root, f), os.path.join(dest, f))
                    names.append(os.path.normpath(os.path.join(rel, f)))
            return sorted(names)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def close(self):
        """end the shared connection (and remove the password helper)."""
        if getattr(self, "_helper", None) and os.path.isfile(self._helper):
            os.remove(self._helper)
            self._helper = None
        subprocess.run([self.ssh_cmd, "-p", str(self.port), "-o", f"ControlPath={self.control}", "-O", "exit", self.target],
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=10)
