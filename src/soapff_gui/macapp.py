"""A double-clickable macOS application for soapff-GUI: ~/Applications/soapff-GUI.app starts this Python with
'-m soapff_gui' (Finder, Launchpad, Dock).  Built with osacompile (an AppleScript applet: a shell-script bundle is
refused by recent macOS); its icon is replaced by data/icon.png (sips, iconutil) and the bundle signed again ad hoc
(codesign -s -), as changing it breaks the signature of osacompile.  soapff-gui --install-app [folder]"""
import os
import shutil
import subprocess
import sys
import tempfile

ICON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "icon.png")


def set_icon(app, png=ICON):
    """data/icon.png -> app/Contents/Resources/applet.icns (all the sizes of an iconset); the applet's asset catalog
    icon (Assets.car, CFBundleIconName), which macOS would use first, is removed."""
    res = os.path.join(app, "Contents", "Resources")
    with tempfile.TemporaryDirectory() as tmp:
        iconset = os.path.join(tmp, "icon.iconset")
        os.makedirs(iconset)
        for n in (16, 32, 128, 256, 512):
            for scale, name in ((1, f"icon_{n}x{n}.png"), (2, f"icon_{n}x{n}@2x.png")):
                subprocess.run(["sips", "-z", str(n * scale), str(n * scale), png, "--out", os.path.join(iconset, name)],
                               check=True, capture_output=True)
        subprocess.run(["iconutil", "-c", "icns", iconset, "-o", os.path.join(res, "applet.icns")], check=True)
    if os.path.exists(os.path.join(res, "Assets.car")):
        os.remove(os.path.join(res, "Assets.car"))
    subprocess.run(["plutil", "-remove", "CFBundleIconName", os.path.join(app, "Contents", "Info.plist")],
                   capture_output=True)


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
    try:
        set_icon(app)
    except (OSError, subprocess.CalledProcessError) as e:       # the app works with the default icon
        print(f"icon not set: {e}")
    subprocess.run(["codesign", "--force", "--deep", "-s", "-", app], check=True, capture_output=True)
    subprocess.run(["touch", app])                               # Finder and the Dock read the new icon
    return app
