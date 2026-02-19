import numpy as np
from pymol import cmd

DEFAULT_D = 1.50  # Å


def _domain_dirs(domains: int) -> np.ndarray:
    if domains == 2:
        return np.array([[0, 0, 1], [0, 0, -1]], float)
    if domains == 3:
        return np.array(
            [[1, 0, 0], [-0.5, np.sqrt(3) / 2, 0], [-0.5, -np.sqrt(3) / 2, 0]], float
        )
    if domains == 4:
        V = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
        return V / np.linalg.norm(V[0])
    if domains == 5:
        return np.array(
            [
                [0, 0, 1],
                [0, 0, -1],
                [1, 0, 0],
                [-0.5, np.sqrt(3) / 2, 0],
                [-0.5, -np.sqrt(3) / 2, 0],
            ],
            float,
        )
    if domains == 6:
        return np.array(
            [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float
        )
    if domains == 7:
        eq = []
        for k in range(5):
            ang = 2 * np.pi * k / 5.0
            eq.append([np.cos(ang), np.sin(ang), 0.0])
        return np.array([[0, 0, 1], [0, 0, -1]] + eq, float)
    raise ValueError("Supported electron-domain counts are 2..7.")


def _lp_sites(domains: int, nlig: int, nlps: int) -> list[int]:
    if nlps <= 0:
        return []
    if domains == 3:
        pref = [0, 1, 2]
    elif domains == 4:
        pref = [0, 1, 2, 3]
    elif domains == 5:
        pref = [2, 3, 4, 0, 1]  # equatorial first
    elif domains == 6:
        # indices: 0 +x,1 -x,2 +y,3 -y,4 +z,5 -z
        if (nlig, nlps) == (5, 1):  # square pyramidal
            pref = [4, 5, 0, 1, 2, 3]
        elif (nlig, nlps) == (4, 2):  # square planar
            pref = [4, 5, 0, 1, 2, 3]
        elif (nlig, nlps) == (3, 3):  # T-shaped (deterministic)
            pref = [4, 5, 3, 2, 0, 1]
        elif (nlig, nlps) == (2, 4):  # linear
            pref = [0, 1, 2, 3, 4, 5]
        else:
            pref = [4, 5, 0, 1, 2, 3]
    elif domains == 7:
        pref = [0, 1, 2, 3, 4, 5, 6]  # axial first
    else:
        pref = list(range(domains))
    return pref[:nlps]


def _parse_lig_colon_string(lig: str, default_d: float):
    """
    lig format requested: "F:F:F:F:" (trailing ':' allowed)
    Also accepts "F:1.56:F:1.56:" etc.

    Each ligand token may be:
      - "F"                  -> uses default_d
      - "F@1.56"             -> explicit bond length
      - "F@1.56@2"           -> explicit bond length + 'order' label (not true bond order storage)
    We use '@' inside each token so we can still split ligands by ':' cleanly.
    """
    if lig is None:
        return []
    s = str(lig).strip().strip('"').strip("'")

    # split ligands on ':' and drop empties (so trailing ':' is fine)
    toks = [t.strip() for t in s.split(":") if t.strip()]
    ligs = []

    for t in toks:
        # per-ligand fields separated by '@'
        seg = [x.strip() for x in t.split("@")]
        elem = seg[0]
        if not elem:
            raise ValueError(f"Bad ligand token '{t}'")

        length = float(default_d)
        order = 1
        if len(seg) >= 2 and seg[1]:
            length = float(seg[1])
        if len(seg) >= 3 and seg[2]:
            order = int(seg[2])

        ligs.append((elem, float(length), int(order)))

    return ligs


def create_vsepr(*, obj, center, lig, lp=0, d=DEFAULT_D, _self=None, **_kw):
    """
    Intended PyMOL syntax (works in your parser mode IF you quote values):
      create_vsepr obj="SF4", center="S", lig="F:F:F:F:", lp=1, d=1.50

    lig uses ':' to separate ligands; trailing ':' allowed.
    To specify per-ligand distance/order, use '@' inside the ligand token:
      lig="F@1.56:F@1.56:F@1.56:F@1.56:"            (length)
      lig="O@1.23@2:O@1.23@2:"                       (length + order label)
    """
    _cmd = _self if _self is not None else cmd

    # sanitize object name (this doesn't create anything; it just makes it legal)
    obj = _cmd.get_legal_name(str(obj).strip().strip('"').strip("'"))
    center = str(center).strip().strip('"').strip("'")

    lp = int(lp)
    d = float(d)

    ligs = _parse_lig_colon_string(lig, d)
    nlig = len(ligs)
    domains = nlig + lp
    if not (2 <= domains <= 7):
        raise ValueError(f"electron domains = {domains}; supported 2..7")

    dirs = _domain_dirs(domains)
    lp_idx = set(_lp_sites(domains, nlig, lp))
    ligand_sites = [i for i in range(domains) if i not in lp_idx]
    if len(ligand_sites) != nlig:
        raise RuntimeError("Internal error: ligand site assignment mismatch")

    # recreate object
    if _cmd.count_atoms(obj) > 0:
        _cmd.delete(obj)

    # native construction at origin
    _cmd.pseudoatom(
        object=obj, selection="none", name="CEN", elem=center, pos=[0.0, 0.0, 0.0]
    )

    for k, site in enumerate(ligand_sites, start=1):
        elem, bl, order = ligs[k - 1]
        pos = (dirs[site] * bl).tolist()
        name = f"L{k}"
        _cmd.pseudoatom(
            object=obj, selection="none", name=name, elem=str(elem), pos=pos
        )
        _cmd.bond(f"({obj} and name CEN)", f"({obj} and name {name})")
        if order != 1:
            _cmd.alter(f"({obj} and name {name})", f'label="order {order}"')

    _cmd.show("sticks", obj)
    _cmd.show("spheres", f"({obj} and name CEN)")
    _cmd.orient(obj)


cmd.extend("create_vsepr", create_vsepr)
