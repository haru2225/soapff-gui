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
ANGLE_MOM = (800.0, 100.0)                      # M-O-M harmonic restraint (k kcal/mol/rad^2, theta0 deg) of the particle edges
                                                # (JAEA file, "analogous to Ho et al."); only with the Underwood & Bourg charges,
                                                # for the triples with an edge metal (_e)
WATER_Q = (Q["ow"], Q["hw"])                    # flexible SPC

# ------------------------------------------------------------------------------------------ surface and edge types
# Bond-valence rule of the edge procedure of Lammers et al. (used by Underwood & Bourg):  q_O = -2 + sum_i (Z_i - q_i) / CN_i  over
# the cations (and H) bonded to the O.  With the ClayFF charges it gives every O type of Table 1 (ob, obos, obts, obss, oh, ohs, and
# the O of talc, brucite, kaolinite ...), so it extends ClayFF to O that lost cations at a surface or an edge:
#   oe_<codes>   an O that lost bonded cations; codes = the cations (and H) it keeps, in the order below, e.g. oe_sh = Si-OH,
#                oe_sa = Si-O-Al(oct), oe_ah = Al(oct)-OH, oe_th = Al(tet)-OH, oe_ta = Al(tet)-O-Al(oct), oe_ggh = Mg2-OH (brucite)
#   <metal>_e<n> an octahedral metal whose n singly coordinated OH2 were removed (Underwood & Bourg); it takes their charge,
#                n (Z - q) / CN, e.g. ao_e1 = 1.575 + 0.2375
# LJ terms: those of the parent type (all O are alike in ClayFF; the metal's own).
BV_CODE = {"s": "st", "t": "at", "a": "ao", "m": "mgo", "g": "mgh", "c": "cao", "k": "cah", "f": "feo", "l": "lio", "h": "ho"}
BV_ZCN = {"st": (4, 4), "at": (3, 4), "ao": (3, 6), "mgo": (2, 6), "mgh": (2, 6), "cao": (2, 6), "cah": (2, 6), "feo": (3, 6),
          "lio": (1, 6), "ho": (1, 1)}
CODE_OF = {t: c for c, t in BV_CODE.items()}
# type names of the JAEA particle file (clay_montmorillonite.data) -> the names above
EDGE_ALIAS = {"oed": "oe_sh", "oaa": "oe_sa", "oab": "oe_ah", "omb": "oe_mh", "oma": "oe_sm", "alb": "ao_e1", "ala": "ao_e1",
              "mgb": "mgo_e1", "mga": "mgo_e1"}


def edge_type(name):
    """("O", [cation types kept]) for oe_<codes>, ("M", base, n) for <base>_e<n>, None otherwise"""
    if name.startswith("oe_") and len(name) > 3 and all(c in BV_CODE for c in name[3:]):
        return ("O", [BV_CODE[c] for c in name[3:]])
    base, sep, n = name.rpartition("_e")
    if sep and base in BV_ZCN and base != "ho" and n.isdigit():
        return ("M", base, int(n))
    return None


# the bonding of the O types of Table 1 (the same charge by the bond-valence rule): an O at a surface that keeps one of these
# sets of cations is that type
BULK_O = {"oe_ss": "ob", "oe_saa": "ob", "oe_sggg": "ob", "oe_sam": "obos", "oe_st": "obts", "oe_taa": "obts",
          "oe_tggg": "obts", "oe_tam": "obss", "oe_aah": "oh", "oe_gggh": "oh", "oe_amh": "ohs"}


def o_type_name(cations):
    """type of an O bonded to these cation (and ho) types: the Table 1 type of that bonding (BULK_O) or oe_<codes>"""
    order = list(BV_CODE.values())
    name = "oe_" + "".join(CODE_OF[t] for t in sorted(cations, key=order.index))
    return BULK_O.get(name, name)


def bv_charge(name, cation_q):
    """charge of an oe_/_e type from the cation charges cation_q (type -> charge) by the bond-valence rule"""
    e = edge_type(name)
    share = lambda t: (BV_ZCN[t][0] - cation_q[t]) / BV_ZCN[t][1]
    if e[0] == "O":
        return -2.0 + sum(share(t) for t in e[1])
    return cation_q[e[1]] + e[2] * share(e[1])


def edge_parent(name):
    e = edge_type(EDGE_ALIAS.get(name, name))
    return None if e is None else ("ob" if e[0] == "O" else e[1])


def is_edge_type(name):
    return edge_parent(name) is not None


def is_edge_metal(name):
    e = edge_type(EDGE_ALIAS.get(name, name))
    return e is not None and e[0] == "M"


def is_edge_hydroxyl(name):
    e = edge_type(EDGE_ALIAS.get(name, name))
    return e is not None and e[0] == "O" and "ho" in e[1]


def register_types(names):
    """add the oe_/_e types among names to Q, LJ and the element table (ClayFF charges); unknown names are left alone"""
    from ._clayff.ff import ELEMENT
    for n in set(map(str, names)):
        par = edge_parent(n)
        if par is None or n in LJ:
            continue
        LJ[n] = LJ[par]
        ELEMENT[n] = ELEMENT[par]
        Q[n] = bv_charge(EDGE_ALIAS.get(n, n), Q)

# charges: only the values of ClayFF Table 1 as published (4 decimals; Cygan et al. 2004), one per type, the ions as in the table.
# The charge of every atom is the charge of its type (the type comes from SOAP).  The net charge is not corrected: Table 1 leaves
# +0.0002 e per octahedral Mg and -0.0002 e per tetrahedral Al substitution.
CHARGE_SCHEMES = {"clayff2004": "ClayFF (Cygan et al. 2004, Table 1), the published values"}
# water models (O, H charges); the O-H bond and H-O-H angle of ClayFF (r0 1.0 A, 109.47 deg, the geometry of SPC and SPC/E):
# SPC flexible (ClayFF), SPC/E rigid (SHAKE on the water bond and angle, shake_line)
WATER_MODELS = {"spc": (-0.82, 0.41), "spce": (-0.8476, 0.4238)}
WATER_LABELS = {"spc": "flexible SPC (ClayFF)", "spce": "rigid SPC/E (SHAKE)"}
RIGID_WATER = ("spce",)


def shake_line(bond_type, angle_type, fix_id="shake1"):
    """the LAMMPS fix that keeps the water molecules rigid (water O-H bond type and H-O-H angle type)"""
    return f"fix {fix_id} all shake 0.0001 20 0 b {bond_type} a {angle_type}"


def type_charges(scheme="clayff2004", water="spc", types=()):
    """type -> charge: the values of ClayFF Table 1 as published, and the water model.  The oe_/_e types among `types` (edges,
    made only by make_slab(allow_bond_cuts=True), not used for now) by the bond-valence rule of Lammers et al. 2017 with them."""
    if scheme not in CHARGE_SCHEMES:
        raise ValueError(f"charge scheme must be one of {list(CHARGE_SCHEMES)}")
    qO, qH = WATER_MODELS[water]
    out = dict({t: v[0] for t, v in TABLE1.items()}, ow=qO, hw=qH)
    for t in set(map(str, types)):
        if is_edge_type(t):
            out[t] = bv_charge(EDGE_ALIAS.get(t, t), out)
    return out
