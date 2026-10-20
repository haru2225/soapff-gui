"""A double-clickable macOS application for soapff-GUI: ~/Applications/soapff-GUI.app starts this Python with
'-m soapff_gui' (Finder, Launchpad, Dock).  Built with osacompile (an AppleScript applet: a shell-script bundle is
refused by recent macOS).  soapff-gui --install-app [folder]"""
import os
import shutil
import subprocess
import sys


def install(folder=None):
    folder = os.path.expanduser(folder or "~/Applications")
    os.makedirs(folder, exist_ok=True)
    app = os.path.join(folder, "soapff-GUI.app")
    if os.path.isdir(app) and os.path.isfile(os.path.join(app, "Contents", "Info.plist")):
        shutil.rmtree(app)                       # an earlier soapff-GUI.app (made by this function)
    py = sys.executable
    # apps started from Finder get a short PATH: put this Python's bin (lmp, mpirun of the same environment) first
    sh = (f"export PATH={os.path.dirname(py)}:/usr/local/bin:/opt/homebrew/bin:$PATH; "
          f"nohup '{py}' -m soapff_gui > /dev/null 2>&1 &")
    script = 'do shell script "' + sh.replace("\\", "\\\\").replace('"', '\\"') + '"'
    subprocess.run(["osacompile", "-o", app, "-e", script], check=True)
    return app
