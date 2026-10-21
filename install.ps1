# Install soapff-GUI on Windows (PowerShell):  powershell -ExecutionPolicy Bypass -File soapff-gui\install.ps1
# A Python 3.10-3.12 from python.org (it has tkinter) is needed; the py launcher finds it.  LAMMPS is installed separately.
# The environment goes to %USERPROFILE%\.soapff-gui\venv; a shortcut "soapff-GUI" is put on the desktop.
$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Venv = if ($env:SOAPFF_VENV) { $env:SOAPFF_VENV } else { Join-Path $env:USERPROFILE ".soapff-gui\venv" }
$Py = $null
foreach ($v in "3.12", "3.11", "3.10", "3.13") {
    try { & py "-$v" -c "import tkinter, sys; assert tkinter.TkVersion >= 8.6" 2>$null; if ($LASTEXITCODE -eq 0) { $Py = @("py", "-$v"); break } } catch {}
}
if (-not $Py) { Write-Error "No Python 3.10-3.13 with tkinter was found. Install one from https://www.python.org (with the py launcher) and run this again." }
Write-Host "Python: $($Py -join ' ')"
& $Py[0] $Py[1] -m venv $Venv
$VPy = Join-Path $Venv "Scripts\python.exe"
& $VPy -m pip install -q --upgrade pip
& $VPy -m pip install -q -e $Here
$Exe = Join-Path $Venv "Scripts\soapff-gui.exe"
& $VPy -c "import sys; assert sys.version_info >= (3, 11)" 2>$null
if ($LASTEXITCODE -eq 0) { & $VPy -m pip install -q sovapy; if ($LASTEXITCODE -ne 0) { Write-Host "WARNING: pip install sovapy failed (PDF / RDF / S(Q) are not available)" } }
else { Write-Host "WARNING: SOVA (sovapy) needs Python 3.11 or newer (PDF / RDF / S(Q) are not available)" }
Write-Host "NOTE: LAMMPS is not installed here: get the Windows installer from https://packages.lammps.org/windows.html (or use WSL2 and install.sh); OVITO from https://www.ovito.org."
Write-Host "installed: $Exe"
$Lnk = Join-Path ([Environment]::GetFolderPath("Desktop")) "soapff-GUI.lnk"
$Sh = (New-Object -ComObject WScript.Shell).CreateShortcut($Lnk)
$Sh.TargetPath = Join-Path $Venv "Scripts\pythonw.exe"
$Sh.Arguments = "-m soapff_gui"
$Sh.WorkingDirectory = $env:USERPROFILE
$Sh.Save()
Write-Host "start it from the desktop shortcut soapff-GUI, or with: $Exe"
