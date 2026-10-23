"""Element data shared by every step: formal charges, roles, bond cut-offs, masses.

Scope: Si, Al, Mg, O, H, Na, K, Ca (v1) + Cs (adsorbate / interlayer cation) and Cl (electrolyte anion, never bonded).
"""
FORMAL = {"Si": 4, "Al": 3, "Mg": 2, "Ca": 2, "Na": 1, "K": 1, "Cs": 1, "H": 1, "O": -2, "Cl": -1,
          "Ti": 4, "Fe": 3}   # Ti(IV), Fe(III) added for the general charge assignment (no ClayFF type)
ALLOWED = set(FORMAL)

# framework cations are part of the bonded mineral network: a cut leaves them under-coordinated, so they are
# removed or capped with O.  Exchangeable cations (interlayer / cavity) are never capped; they are moved or removed
# to keep the slab neutral and non-polar.
FRAMEWORK = ("Si", "Al", "Mg", "Ti", "Fe")   # Ti, Fe appended (general charges; absent from all earlier data)
EXCHANGEABLE = ("Na", "K", "Ca", "Cs")

# M-O bond cut-offs (A).  Identical to the ones used for the DFT environment dictionary (test4/minerals_dict.py),
# so that environment classes computed here match the classes of the training data.
CUT = {"H": 1.2, "Si": 2.0, "Al": 2.2, "Mg": 2.4, "Ca": 2.8, "Na": 2.9, "K": 3.3, "Cs": 3.6,   # Cs added (Cs-O 3.0-3.5 A in micas)
       "Ti": 2.3, "Fe": 2.3}                                                          # Ti-O 1.9-2.0, Fe(III)-O 1.9-2.1 A

MASS = {"H": 1.00794, "O": 15.9994, "Na": 22.98977, "Mg": 24.305, "Al": 26.981538, "Si": 28.0855, "K": 39.0983, "Ca": 40.078,
        "Cs": 132.90545, "Cl": 35.453, "Ti": 47.867, "Fe": 55.845}

# covalent-ish radii used only for clash detection (A)
RCLASH = {"H": 0.35, "O": 0.70, "Si": 0.85, "Al": 0.85, "Mg": 0.85, "Na": 1.0, "K": 1.3, "Ca": 1.0, "Cs": 1.6, "Cl": 1.4, "Ti": 0.9, "Fe": 0.9}

OH_LENGTH = 0.97          # A, length used for newly placed O-H bonds
