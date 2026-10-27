"""Helpers shared by the tests of the libraries with water: the most similar row (same element) of a library for every atom of a structure
(the SOAP of the whole structure, computed on the elements of the library), and the true types of a structure with water."""
import numpy as np

from soapff_cui import core

alias = lambda t: np.array([core.TYPE_ALIAS.get(str(x), str(x)) for x in t])


def nearest(at, lib):
    el = np.array(at.get_chemical_symbols())
    keep = np.nonzero(np.isin(el, lib["species"]))[0]
    Fk = core.soap_vectors(at[keep], lib["species"], **lib["soap"]).astype(np.float32)
    F = np.zeros((len(at), Fk.shape[1]), np.float32); F[keep] = Fk
    X = lib["X"].astype(np.float32)
    ty, sim = np.full(len(at), "", dtype=object), np.zeros(len(at))
    for e in sorted(set(el[keep])):
        m = np.nonzero(el == e)[0]; r = np.nonzero(lib["element"] == e)[0]
        S = F[m] @ X[r].T
        j = S.argmax(1)
        ty[m] = lib["type"][r][j]; sim[m] = S[np.arange(len(m)), j]
    return alias(ty), sim


def truth_rule(at):
    """bond-graph ClayFF types, the water by its shape and the ions by their element (core.roles)"""
    from soapff_cui._clayff.ff import ff_types
    role, _ = core.roles(at)
    s = at.copy()
    s.set_array("role", np.array(["water" if r == "water" else "ion" if r == "ion" else "bulk" for r in role], dtype="U12"))
    t = alias(ff_types(s))
    el = np.array(at.get_chemical_symbols())
    t = np.where(role == "water", np.where(el == "O", "ow", "hw"), t)
    return np.where(role == "ion", el, t)
