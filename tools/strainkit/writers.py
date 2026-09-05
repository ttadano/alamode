# Copyright (c) 2023 Ryota Masuki (strainIFCcoupling,
#                    https://github.com/r-masuki/strainIFCcoupling)
# Copyright (c) 2026 Terumasa Tadano
# MIT license.  See LICENCE.txt of the ALAMODE package.
"""Writers (and token-stream readers) for the strain-related anphon input files.

All readers mimic the whitespace-token parsing of the C++ code (``operator>>``)
so that a file accepted here is parsed identically by anphon.  They are stricter
than anphon in a few places: only decimal numbers are recognised, and trailing
tokens or mixed SOEC/TOEC units are errors here where anphon only warns.

* strain_harmonic.in  -- anphon/ifc_derivative.cpp (calculate_delv2_delumn_finite_difference)
* strain_force.in     -- anphon/ifc_derivative.cpp (calculate_delv1_delumn_finite_difference);
                         forces in eV/Angstrom, one row per atom of the cell the strained
                         calculations used.  An optional "&reference_cell ... /" header names
                         that cell so that anphon can map the rows onto a user-defined &cell.
* elastic_constants.in-- anphon/elastic_tensor.cpp (read_elastic_constants); 81 + 729 values
                         in the full 3x3-index layout (i = 3*mu + nu): C in GPa (labels
                         "SOEC GPa" / "TOEC GPa", valid for any cell) or, in the legacy layout
                         without a unit token, V*C in Ry for one specific cell.
* C1_array.in         -- anphon/elastic_tensor.cpp (read_C1_array); 9 values, sigma in GPa
                         ("C1 GPa") or V*sigma in Ry (legacy, no unit token).
"""

import math
import re

import numpy as np

from .strain import mode_pair
from .units import BOHR_IN_ANGSTROM

_VALID_MODES = ("xx", "yy", "zz", "xy", "yz", "zx")


def _tokens(path):
    with open(path) as f:
        return f.read().split()


_NUMBER = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")


def _is_number(tok):
    """True if anphon reads ``tok`` as a number: a decimal literal consumed
    completely by strtod without overflow (so ``1_0`` and ``1e999`` are not)."""
    return bool(_NUMBER.fullmatch(tok)) and math.isfinite(float(tok))


def _float(tok, path):
    if not _is_number(tok):
        raise ValueError(f"{path}: {tok!r} is not a number")
    return float(tok)


def _check_finite(arr, what):
    arr = np.asarray(arr, dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{what}: non-finite value encountered")
    return arr


# ------------------------------------------------------------ strain_harmonic.in
def write_strain_harmonic_in(path, rows):
    """rows: iterable of (mode, smag, weight, filename)."""
    lines = []
    for mode, smag, weight, filename in rows:
        mode_pair(mode)
        if not np.isfinite(smag) or not np.isfinite(weight) or smag == 0.0:
            raise ValueError(f"strain_harmonic.in: invalid smag/weight for mode {mode}")
        if any(c.isspace() for c in str(filename)):
            raise ValueError(
                "strain_harmonic.in: filenames must not contain whitespace"
            )
        lines.append(
            "{0:4s} {1:25.15f} {2:25.15f} {3:25s}\n".format(
                mode, float(smag), float(weight), str(filename)
            )
        )
    with open(path, "w") as f:
        f.writelines(lines)


def read_strain_harmonic_in(path):
    tok = _tokens(path)
    if len(tok) % 4 != 0:
        raise ValueError(f"{path}: number of tokens is not a multiple of 4")
    rows = []
    for k in range(0, len(tok), 4):
        mode = tok[k]
        if mode not in _VALID_MODES:
            raise ValueError(f"{path}: invalid mode name {mode!r}")
        rows.append((mode, float(tok[k + 1]), float(tok[k + 2]), tok[k + 3]))
    return rows


# --------------------------------------------------------------- strain_force.in
class ReferenceCell:
    """The cell declared in the ``&reference_cell`` header of strain_force.in.

    ``lavec``: (3, 3) lattice vectors as rows in Angstrom; ``elements``: element
    symbol per atom; ``xf``: (natom, 3) fractional coordinates.  Any object with
    these three attributes (e.g. :class:`strainkit.fcsorder.AnphonPrimitive`) is
    accepted wherever a reference cell is expected.
    """

    def __init__(self, lavec, elements, xf):
        self.elements = [str(e) for e in elements]
        self.lavec = np.asarray(lavec, dtype=float).reshape(3, 3)
        self.xf = np.asarray(xf, dtype=float).reshape(len(self.elements), 3)

    @property
    def natom(self):
        return len(self.elements)


def _reference_cell_lines(cell, natmin):
    """Header lines.  No comments: anphon reads the header as a plain token stream."""
    elements = [str(e) for e in cell.elements]
    lavec = _check_finite(cell.lavec, "strain_force.in (reference cell lattice)")
    xf = _check_finite(cell.xf, "strain_force.in (reference cell coordinates)")
    if lavec.shape != (3, 3) or xf.shape != (len(elements), 3):
        raise ValueError(
            "strain_force.in: the reference cell needs lavec (3,3), elements (natom,) "
            "and xf (natom,3)"
        )
    if len(elements) != natmin:
        raise ValueError(
            f"strain_force.in: the reference cell has {len(elements)} atoms but the "
            f"force blocks have {natmin} rows"
        )
    for e in elements:
        if not e or e == "/" or any(c.isspace() for c in e):
            raise ValueError(f"strain_force.in: invalid element symbol {e!r}")
    lines = ["&reference_cell\n", f"  {1.0 / BOHR_IN_ANGSTROM:.16g}\n"]
    for v in lavec:
        lines.append("  {0:22.15f} {1:22.15f} {2:22.15f}\n".format(*v))
    lines.append(f"  {natmin}\n")
    for e, x in zip(elements, xf):
        lines.append("  {0:4s} {1:20.15f} {2:20.15f} {3:20.15f}\n".format(e, *x))
    lines.append("/\n")
    return lines


def _parse_reference_cell(tok, path):
    """Parse the header at tok[0]; returns (ReferenceCell, index of the first block token)."""
    nhead = 12  # sentinel, scale, 9 lattice components, natom
    if len(tok) < nhead:
        raise ValueError(f"{path}: truncated &reference_cell header")
    try:
        scale = _float(tok[1], path)
        lav = np.array([_float(t, path) for t in tok[2:11]]).reshape(3, 3)
        natom = int(tok[11])
    except ValueError as exc:
        raise ValueError(f"{path}: malformed &reference_cell header ({exc})") from None
    if natom < 1 or len(tok) < nhead + 4 * natom + 1:
        raise ValueError(
            f"{path}: the &reference_cell header declares {natom} atoms but is truncated"
        )
    elements, xf = [], []
    for i in range(natom):
        k = nhead + 4 * i
        elements.append(tok[k])
        try:
            xf.append([_float(t, path) for t in tok[k + 1 : k + 4]])
        except ValueError:
            raise ValueError(
                f"{path}: malformed atom line {i + 1} in the &reference_cell header"
            ) from None
    k = nhead + 4 * natom
    if tok[k] != "/":
        raise ValueError(f"{path}: the &reference_cell header is not terminated by '/'")
    return ReferenceCell(scale * lav * BOHR_IN_ANGSTROM, elements, xf), k + 1


def write_strain_force_in(path, blocks, reference_cell=None):
    """blocks: iterable of (mode, smag, weight, forces) with forces (natmin, 3) in eV/A.

    ``reference_cell`` (optional): the cell the force rows belong to, an object with
    ``lavec`` (rows, Angstrom), ``elements`` and ``xf`` attributes.  When given, the
    file starts with an ``&reference_cell ... /`` header, from which anphon maps the
    rows onto the atoms of a user-defined ``&cell`` (a nested supercell or sub-cell)
    instead of assuming that the rows already match its primitive cell.
    """
    lines = []
    natmin = None
    for mode, smag, weight, forces in blocks:
        mode_pair(mode)
        forces = _check_finite(forces, f"strain_force.in ({mode})")
        if forces.ndim != 2 or forces.shape[1] != 3:
            raise ValueError("strain_force.in: forces must have shape (natmin, 3)")
        if natmin is None:
            natmin = forces.shape[0]
        elif forces.shape[0] != natmin:
            raise ValueError(
                "strain_force.in: inconsistent number of atoms between blocks"
            )
        if not np.isfinite(smag) or smag == 0.0 or not np.isfinite(weight):
            raise ValueError(f"strain_force.in: invalid smag/weight for mode {mode}")
        lines.append(
            "{0:4s} {1:25.15f} {2:25.15f}\n".format(mode, float(smag), float(weight))
        )
        for fx, fy, fz in forces:
            lines.append("{0:25.15f} {1:25.15f} {2:25.15f}\n".format(fx, fy, fz))
    if reference_cell is not None:
        if natmin is None:
            raise ValueError("strain_force.in: no force blocks to write")
        lines = _reference_cell_lines(reference_cell, natmin) + lines
    with open(path, "w") as f:
        f.writelines(lines)


def read_strain_force_in(path, natmin=None):
    """Token-stream reader; returns (blocks, reference_cell).

    ``blocks`` is a list of (mode, smag, weight, forces(natmin, 3)) and
    ``reference_cell`` a :class:`ReferenceCell` when the file carries the header,
    else None.  ``natmin`` may be omitted for files with a header (it is checked
    against the header when given).
    """
    tok = _tokens(path)
    ref, k = None, 0
    if tok and tok[0].lower() == "&reference_cell":
        ref, k = _parse_reference_cell(tok, path)
        if natmin is not None and natmin != ref.natom:
            raise ValueError(
                f"{path}: the &reference_cell header declares {ref.natom} atoms, not {natmin}"
            )
        natmin = ref.natom
    if natmin is None:
        raise ValueError(f"{path}: no &reference_cell header, so natmin must be given")
    blocks = []
    per = 3 + 3 * natmin
    if (len(tok) - k) % per != 0:
        raise ValueError(
            f"{path}: token count {len(tok) - k} is not a multiple of {per} (natmin={natmin})"
        )
    while k < len(tok):
        mode = tok[k]
        if mode not in _VALID_MODES:
            raise ValueError(f"{path}: invalid mode name {mode!r}")
        smag, weight = _float(tok[k + 1], path), _float(tok[k + 2], path)
        vals = np.array([_float(t, path) for t in tok[k + 3 : k + per]]).reshape(
            natmin, 3
        )
        blocks.append((mode, smag, weight, vals))
        k += per
    return blocks, ref


# ------------------------------------------------------- elastic_constants.in
UNITS = ("Ry", "GPa")


def normalize_unit(unit):
    """'GPa' (cell-independent densities) or 'Ry' (V*C per cell, the legacy convention)."""
    for u in UNITS:
        if str(unit).lower() == u.lower():
            return u
    raise ValueError(f"unknown unit {unit!r} (expected one of {UNITS})")


def _label(name, unit):
    # 'Ry' files carry no unit token, exactly as before the unit token existed, so
    # that older anphon versions can still read them; only GPa files are tagged.
    return f"{name}\n" if unit == "Ry" else f"{name} {unit}\n"


def _read_section(tok, k, name, nvalues, path):
    """Read 'LABEL [unit] v1 ... vn' starting at tok[k]; returns (values, unit, k_next).

    As in anphon, the label text is not checked; the token after it is the unit
    when it is not a number, otherwise it is the first value (legacy file, 'Ry').
    """
    if k >= len(tok):
        raise ValueError(f"{path}: missing {name} section")
    k += 1
    unit = "Ry"
    if k < len(tok) and not _is_number(tok[k]):
        unit = normalize_unit(tok[k])
        k += 1
    if len(tok) < k + nvalues:
        raise ValueError(f"{path}: {name} section has fewer than {nvalues} values")
    return np.array([_float(t, path) for t in tok[k : k + nvalues]]), unit, k + nvalues


def write_elastic_constants_in(path, c2_99, c3_999, unit="Ry"):
    """c2: (9,9) and c3: (9,9,9) arrays (row-major i = 3*mu+nu).

    ``unit="GPa"``: C in GPa, valid for any anphon cell (labels ``SOEC GPa`` /
    ``TOEC GPa``).  ``unit="Ry"`` (default): V*C in Ry for one specific cell,
    written in the legacy layout without a unit token.
    """
    unit = normalize_unit(unit)
    c2 = _check_finite(c2_99, "elastic_constants.in (SOEC)")
    c3 = _check_finite(c3_999, "elastic_constants.in (TOEC)")
    if c2.shape != (9, 9) or c3.shape != (9, 9, 9):
        raise ValueError("elastic_constants.in: expected shapes (9,9) and (9,9,9)")
    with open(path, "w") as f:
        f.write(_label("SOEC", unit))
        for v in c2.reshape(-1):
            f.write(f"{v:.12e}\n")
        f.write(_label("TOEC", unit))
        for v in c3.reshape(-1):
            f.write(f"{v:.12e}\n")


def read_elastic_constants_in(path):
    """Returns (c2 (9,9), c3 (9,9,9), unit) with unit 'GPa' or 'Ry' (legacy, per cell).

    Stricter than anphon, which converts the two sections independently and only
    warns: different SOEC/TOEC units and trailing tokens are errors here.
    """
    tok = _tokens(path)
    c2, u2, k = _read_section(tok, 0, "SOEC", 81, path)
    c3, u3, k = _read_section(tok, k, "TOEC", 729, path)
    if k != len(tok):
        raise ValueError(f"{path}: {len(tok) - k} unexpected trailing tokens")
    if u2 != u3:
        raise ValueError(f"{path}: SOEC ({u2}) and TOEC ({u3}) use different units")
    return c2.reshape(9, 9), c3.reshape(9, 9, 9), u2


# ------------------------------------------------------------------ C1_array.in
def write_C1_array_in(path, sigma_3x3, unit="Ry"):
    """3x3 reference stress: sigma in GPa (``unit="GPa"``) or V*sigma in Ry (default,
    legacy layout without a unit token)."""
    unit = normalize_unit(unit)
    s = _check_finite(sigma_3x3, "C1_array.in")
    if s.shape != (3, 3):
        raise ValueError("C1_array.in: expected a 3x3 array")
    with open(path, "w") as f:
        f.write(_label("C1", unit))
        for v in s.reshape(-1):
            f.write(f"{v:.12e}\n")


def read_C1_array_in(path):
    """Returns (sigma (3,3), unit) with unit 'GPa' or 'Ry' (legacy, V*sigma per cell)."""
    tok = _tokens(path)
    s, unit, k = _read_section(tok, 0, "C1", 9, path)
    if k != len(tok):
        raise ValueError(f"{path}: {len(tok) - k} unexpected trailing tokens")
    return s.reshape(3, 3), unit


# ------------------------------------------------------- anphon-style checks
def asymmetry_rank2(s):
    """max |s_ij - s_ji| (Relaxation::set_elastic_constants warns above 1e-6)."""
    s = np.asarray(s, dtype=float).reshape(3, 3)
    return float(np.abs(s - s.T).max())


def asymmetry_rank4(c99):
    """max deviation from the intrinsic symmetries of C2 (9x9 layout)."""
    c = np.asarray(c99, dtype=float).reshape(3, 3, 3, 3)
    dev = 0.0
    for perm in ((1, 0, 2, 3), (0, 1, 3, 2), (2, 3, 0, 1)):
        dev = max(dev, float(np.abs(c - np.transpose(c, perm)).max()))
    return dev


def asymmetry_rank6(c999):
    c = np.asarray(c999, dtype=float).reshape((3,) * 6)
    dev = 0.0
    for perm in (
        (1, 0, 2, 3, 4, 5),
        (0, 1, 3, 2, 4, 5),
        (0, 1, 2, 3, 5, 4),
        (2, 3, 0, 1, 4, 5),
        (0, 1, 4, 5, 2, 3),
    ):
        dev = max(dev, float(np.abs(c - np.transpose(c, perm)).max()))
    return dev


def anphon_file_locations():
    return (
        "Where anphon expects these files:\n"
        "  STRAIN_IFC_DIR/  <- elastic_constants.in, strain_harmonic.in (+ the FC files it lists),\n"
        "                      strain_force.in\n"
        "  working dir      <- C1_array.in (read from the directory anphon runs in)\n"
    )
