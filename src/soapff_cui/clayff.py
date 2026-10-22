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
    # particle-edge types of the JAEA montmorillonite particle (clay_montmorillonite.data); not in ClayFF Table 1.  LJ terms of the
    # parent type (EDGE_PARENT).  Charges from the bond-valence rule of the edge procedure of Lammers et al.,
    # q_O = -2 + sum_i (Z_i - q_i)/CN_i over the bonded neighbours (H 1/0.425/1, Si 4/2.1/4, Al 3/1.575/6, Mg 2/1.36/6) with the
    # ClayFF charges: oed Si-OH, oaa Si-O-Al, oab Al-OH, omb Mg-OH, oma Si-O-Mg; the edge metal takes the share of the removed OH2
    # (Al 1.575 + 0.2375, Mg 1.36 + 0.1067).  (The JAEA file itself uses Mg 1.361, O next to it -1.081 / -1.181 and Mg edge 1.4675,
    # which makes that particle exactly neutral; here the ClayFF values are kept and the small residual is removed by make().)
    "oed": (-0.9500, 0.1554, 3.5532), "oab": (-1.1875, 0.1554, 3.5532), "omb": (-1.3183, 0.1554, 3.5532),
    "oaa": (-1.2875, 0.1554, 3.5532), "oma": (-1.4183, 0.1554, 3.5532),
    "alb": (1.8125, 1.3298e-6, 4.7943), "ala": (1.8125, 1.3298e-6, 4.7943),
    "mgb": (1.4667, 9.0298e-7, 5.9090), "mga": (1.4667, 9.0298e-7, 5.9090),
}
EDGE_PARENT = {"oed": "ob", "oab": "ob", "omb": "ob", "oaa": "ob", "oma": "ob", "alb": "ao", "ala": "ao", "mgb": "mgo", "mga": "mgo"}
SIGMA_FACTOR = 2.0 ** (1.0 / 6.0)
Q = {t: v[0] for t, v in TABLE1.items()}
LJ = {t: (v[1], v[2] / SIGMA_FACTOR) for t, v in TABLE1.items()}          # (epsilon kcal/mol, sigma A) for LAMMPS
BOND_OH = (554.1349, 1.0000)                    # o*-h*, oh-ho, ohs-ho (k1 kcal/mol/A^2, r0 A)
ANGLE_HOH = (45.7696, 109.47)                   # h*-o*-h* (k2 kcal/mol/rad^2, theta0 deg)
ANGLE_MOH = (30.0, 109.47)                      # Metal-oh-ho and Metal-ohs-ho
ANGLE_MOM = (800.0, 100.0)                      # M-O-M harmonic restraint (k kcal/mol/rad^2, theta0 deg) of the particle edges
                                                # (JAEA file, "analogous to Ho et al."); only for triples with an edge metal
WATER_Q = (Q["ow"], Q["hw"])                    # flexible SPC
