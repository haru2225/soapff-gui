"""ClayFF (Cygan, Liang, Kalinichev, J. Phys. Chem. B 108, 1255 (2004)), Tables 1 and 2, exactly as published, with the
flexible SPC water of ClayFF.  No other terms are added: the M-O-H bending constant is the one of the original ClayFF.

Table 1, nonbond: charge q (e), D0 (kcal/mol), R0 (A) of  E = D0 [ (R0/r)^12 - 2 (R0/r)^6 ]  (the 12-6 form with
epsilon = D0 and sigma = R0 / 2^(1/6), which is what LAMMPS lj/cut takes).
Table 2, bond stretch  E = k1 (r - r0)^2  and angle bend  E = k2 (theta - theta0)^2.
Type names: o* -> ow, h* -> hw (water), the others as in the table.  Na and Cl- are Smith and Dang's, K and Ca Koneshan
et al.'s, Cs Smith and Dang's, Ba Aqvist's (the table's footnotes)."""

# type: (charge, D0 kcal/mol, R0 A)
TABLE1 = {
    "hw": (0.4100, 0.0, 0.0), "ho": (0.4250, 0.0, 0.0),
    "ow": (-0.8200, 0.1554, 3.5532), "oh": (-0.9500, 0.1554, 3.5532), "ob": (-1.0500, 0.1554, 3.5532),
    "obos": (-1.1808, 0.1554, 3.5532), "obts": (-1.1688, 0.1554, 3.5532), "obss": (-1.2996, 0.1554, 3.5532),
    "ohs": (-1.0808, 0.1554, 3.5532),
    "st": (2.1000, 1.8405e-6, 3.7064), "ao": (1.5750, 1.3298e-6, 4.7943), "at": (1.5750, 1.8405e-6, 3.7064),
    "mgo": (1.3600, 9.0298e-7, 5.9090), "mgh": (1.0500, 9.0298e-7, 5.9090),
    "cao": (1.3600, 5.0298e-6, 6.2484), "cah": (1.0500, 5.0298e-6, 6.2428),
    "feo": (1.5750, 9.0298e-6, 5.5070), "lio": (0.5250, 9.0298e-6, 4.7257),
    "Na": (1.0, 0.1301, 2.6378), "K": (1.0, 0.1000, 3.7423), "Cs": (1.0, 0.1000, 4.3002),
    "Ca": (2.0, 0.1000, 3.2237), "Ba": (2.0, 0.0470, 4.2840), "Cl": (-1.0, 0.1001, 4.9388),
}
SIGMA_FACTOR = 2.0 ** (1.0 / 6.0)
Q = {t: v[0] for t, v in TABLE1.items()}
LJ = {t: (v[1], v[2] / SIGMA_FACTOR) for t, v in TABLE1.items()}          # (epsilon kcal/mol, sigma A) for LAMMPS
BOND_OH = (554.1349, 1.0000)                    # o*-h*, oh-ho, ohs-ho (k1 kcal/mol/A^2, r0 A)
ANGLE_HOH = (45.7696, 109.47)                   # h*-o*-h* (k2 kcal/mol/rad^2, theta0 deg)
ANGLE_MOH = (30.0, 109.47)                      # Metal-oh-ho and Metal-ohs-ho
WATER_Q = (Q["ow"], Q["hw"])                    # flexible SPC

# montmorillonite edges: only with the optional library built with the JAEA particle (clay_montmorillonite.data, the edge model
# that Underwood & Bourg, J. Phys. Chem. C 124, 3702 (2020) take from Lammers et al.).  The edge types take the LJ of their parent
# ClayFF type (as the JAEA settings file); their charges are not in Table 1: they are those of the JAEA file, kept from the
# library (type_charges() leaves them out).  The edge model adds a bending of each edge octahedral metal - O - octahedral metal
# angle (the JAEA file, angle type 2).
EDGE_PARENT = {"oed": "ob", "oab": "ob", "omb": "ob", "oaa": "ob", "oma": "ob", "alb": "ao", "ala": "ao", "mgb": "mgo", "mga": "mgo"}
EDGE_HYDROXYL_O = ("oed", "oab", "omb")         # edge O carrying one H (Si-OH, Al-OH, Mg-OH at the edge)
EDGE_METALS = ("alb", "ala", "mgb", "mga")       # octahedral Al, Mg at the edge
ANGLE_MOM_EDGE = (800.0, 100.0)                 # edge metal - O - octahedral metal (k2 kcal/mol/rad^2, theta0 deg; JAEA file)
LJ.update({t: LJ[p] for t, p in EDGE_PARENT.items()})

# charges: only the values of ClayFF Table 1 as published (4 decimals; Cygan et al. 2004), one per type, the ions as in the table.
# The charge of every atom is the charge of its type (the type comes from SOAP).  The net charge is not corrected: Table 1 leaves
# +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al substitution.
CHARGE_SCHEMES = {"clayff2004": "ClayFF (Cygan et al. 2004, Table 1), the published values"}
# water model (O, H charges): the flexible SPC water of ClayFF (Table 1 o*, h*; Table 2 bond and angle), the only one.
WATER_MODELS = {"spc": (TABLE1["ow"][0], TABLE1["hw"][0])}
WATER_LABELS = {"spc": "flexible SPC (ClayFF)"}


def type_charges(scheme="clayff2004", water="spc", types=()):
    """type -> charge: the values of ClayFF Table 1 as published, and the water model (`types`: unused, kept for the callers)"""
    if scheme not in CHARGE_SCHEMES:
        raise ValueError(f"charge scheme must be one of {list(CHARGE_SCHEMES)}")
    qO, qH = WATER_MODELS[water]
    return dict({t: v[0] for t, v in TABLE1.items()}, ow=qO, hw=qH)
