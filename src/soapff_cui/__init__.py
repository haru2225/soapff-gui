"""soapff-cui: a structure file in, a LAMMPS ClayFF force field out, with every atomic charge taken from SOAP.

  python -m soapff_cui build-library REF.data [REF2.data ...] -o library.npz
      reference structures with known charges -> table of (SOAP vector, charge, type) for every atom
      (clay, water and ions alike)
  python -m soapff_cui check-library library.npz
      developer check: nearest other entry in the library has the same charge
  python -m soapff_cui make STRUCTURE -o OUTDIR [--library library.npz]
      STRUCTURE: any file ASE reads with a cell (extxyz, cif, POSCAR, ...) or a LAMMPS data file (atom_style full).
      Writes OUTDIR/system.data, system.settings, in.lammps, charges.csv, report.json.
  python -m soapff_cui open OUTDIR
      the trajectory (ClayFF type names and charges) in OVITO
  python -m soapff_cui analyze OUTDIR
      -> OUTDIR/analysis/: self-diffusion coefficients from LAMMPS compute msd (water, ions; total, in plane, normal)
      and PDF g(r) / RDF from SOVA (sovapy) (soapff_cui/analysis.py, soapff_cui/sova_structure.py)

Charges: for every atom the SOAP vector (same element) most similar to it in the library gives the charge, unchanged.
Atoms whose best similarity is below the threshold, or whose element is not in the library, get the ClayFF type
charge instead (water flexible SPC, ions formal) and are flagged in charges.csv / report.json.
Force field: ClayFF exactly as in Cygan et al. (2004) Tables 1 and 2 (soapff_cui/clayff.py: charges, LJ, O-H bond,
H-O-H and M-O-H angle) with its flexible SPC water, 0.5 fs time step (writer and typing: soapff_cui/_clayff, copied from mincharge).
"""

# bundled with soapff-gui (https://github.com/haru2225/soapff-cui is the stand-alone version)
__version__ = "0.1.0"
