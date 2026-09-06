import numpy as np
import pytest

from strainkit import elasticfit as ef
from strainkit import strain, writers


def test_elastic_constants_round_trip(tmp_path):
    rng = np.random.default_rng(0)
    c4 = ef.voigt_to_full2(rng.normal(size=21))
    c6 = ef.voigt_to_full3(rng.normal(size=56))
    p = tmp_path / "elastic_constants.in"
    writers.write_elastic_constants_in(p, ef.full2_to_9x9(c4), ef.full3_to_9x9x9(c6))
    tok = p.read_text().split()
    assert len(tok) == 1 + 81 + 1 + 729 and tok[0] == "SOEC" and tok[82] == "TOEC"
    a, b, unit = writers.read_elastic_constants_in(p)
    assert unit == "Ry"  # legacy layout: no unit token
    assert np.allclose(a, ef.full2_to_9x9(c4)) and np.allclose(b, ef.full3_to_9x9x9(c6))
    assert writers.asymmetry_rank4(a) == 0.0 and writers.asymmetry_rank6(b) == 0.0
    # GPa files carry a unit token after each label
    writers.write_elastic_constants_in(
        p, ef.full2_to_9x9(c4), ef.full3_to_9x9x9(c6), unit="GPa"
    )
    tok = p.read_text().split()
    assert len(tok) == 2 + 81 + 2 + 729
    assert tok[:2] == ["SOEC", "GPa"] and tok[83:85] == ["TOEC", "GPa"]
    a, b, unit = writers.read_elastic_constants_in(p)
    assert unit == "GPa" and np.allclose(a, ef.full2_to_9x9(c4))
    with pytest.raises(ValueError, match="unit"):
        writers.write_elastic_constants_in(
            p, ef.full2_to_9x9(c4), ef.full3_to_9x9x9(c6), unit="MPa"
        )
    p.write_text(p.read_text().replace("SOEC GPa", "SOEC", 1))
    with pytest.raises(ValueError, match="different units"):
        writers.read_elastic_constants_in(p)
    # the unit probe is as strict as anphon's strtod: "1_0" and "1e999" are not numbers
    for bad in ("1_0", "1e999"):
        p.write_text(p.read_text().replace("SOEC", f"SOEC {bad}", 1))
        with pytest.raises(ValueError, match="unit"):
            writers.read_elastic_constants_in(p)
        p.write_text(p.read_text().replace(f"SOEC {bad}", "SOEC", 1))
    with pytest.raises(ValueError):
        writers.write_elastic_constants_in(
            p, np.full((9, 9), np.nan), np.zeros((9, 9, 9))
        )
    with pytest.raises(ValueError):
        writers.write_elastic_constants_in(p, np.zeros((6, 6)), np.zeros((9, 9, 9)))


def test_c1_round_trip(tmp_path):
    s = strain.sym_from_voigt([1, 2, 3, 4, 5, 6])
    p = tmp_path / "C1_array.in"
    writers.write_C1_array_in(p, s)
    assert len(p.read_text().split()) == 10
    r, unit = writers.read_C1_array_in(p)
    assert unit == "Ry"
    assert (
        np.allclose(r, s)
        and r[1, 2] == r[2, 1] == 4
        and writers.asymmetry_rank2(r) == 0
    )
    writers.write_C1_array_in(p, s, unit="GPa")
    tok = p.read_text().split()
    assert len(tok) == 11 and tok[:2] == ["C1", "GPa"]
    r, unit = writers.read_C1_array_in(p)
    assert unit == "GPa" and np.allclose(r, s)


def test_strain_harmonic_in(tmp_path):
    rows = [
        ("xx", 0.005, 0.5, "a.xml"),
        ("xx", -0.005, 0.5, "b.xml"),
        ("yz", 0.005, 1.0, "c.h5"),
    ]
    p = tmp_path / "strain_harmonic.in"
    writers.write_strain_harmonic_in(p, rows)
    back = writers.read_strain_harmonic_in(p)
    assert (
        [r[0] for r in back] == ["xx", "xx", "yz"]
        and back[1][1] == -0.005
        and back[2][3] == "c.h5"
    )
    with pytest.raises(ValueError):
        writers.write_strain_harmonic_in(p, [("xz", 0.005, 1.0, "x.xml")])


def test_strain_force_in(tmp_path):
    rng = np.random.default_rng(1)
    pts = strain.ifc_strain_set(strain.default_modes(), central=True)
    blocks = [(q.mode, q.smag, q.weight, rng.normal(size=(4, 3))) for q in pts]
    p = tmp_path / "strain_force.in"
    writers.write_strain_force_in(p, blocks)
    back, ref = writers.read_strain_force_in(p, 4)
    assert len(back) == 12 and ref is None
    for (m1, s1, w1, f1), (m2, s2, w2, f2) in zip(blocks, back):
        assert m1 == m2 and s1 == s2 and w1 == w2 and np.allclose(f1, f2)
    assert np.allclose(strain.check_weight_sums(pts), 1.0)
    with pytest.raises(ValueError):
        writers.read_strain_force_in(p, 3)
    with pytest.raises(ValueError, match="natmin"):
        writers.read_strain_force_in(p)  # no header: natmin is required
    with pytest.raises(ValueError):
        writers.write_strain_force_in(p, [("xx", 0.005, 1.0, np.zeros((4, 2)))])


def test_strain_force_in_reference_cell(tmp_path):
    rng = np.random.default_rng(2)
    pts = strain.ifc_strain_set(strain.default_modes(), central=True)
    blocks = [(q.mode, q.smag, q.weight, rng.normal(size=(4, 3))) for q in pts]
    cell = writers.ReferenceCell(
        np.diag([3.0, 4.0, 5.0]), ["Zn", "Zn", "O", "O"], rng.uniform(size=(4, 3))
    )
    p = tmp_path / "strain_force.in"
    writers.write_strain_force_in(p, blocks, reference_cell=cell)
    tok = p.read_text().split()
    # sentinel, scale, 9 lattice components, natom, 4 x (symbol x y z), "/"
    assert tok[0] == "&reference_cell" and tok[11] == "4" and tok[12] == "Zn"
    assert tok[28] == "/" and tok[29] == blocks[0][0]
    assert "#" not in p.read_text()  # anphon reads the header as a token stream
    back, ref = writers.read_strain_force_in(p)  # natmin taken from the header
    assert ref.natom == 4 and ref.elements == ["Zn", "Zn", "O", "O"]
    assert np.allclose(ref.lavec, cell.lavec) and np.allclose(ref.xf, cell.xf)
    assert len(back) == 12 and np.allclose(back[0][3], blocks[0][3])
    assert np.allclose(writers.read_strain_force_in(p, 4)[0][5][3], blocks[5][3])
    with pytest.raises(ValueError, match="declares 4 atoms"):
        writers.read_strain_force_in(p, 3)
    with pytest.raises(ValueError, match="atoms"):
        writers.write_strain_force_in(
            p,
            blocks,
            reference_cell=writers.ReferenceCell(np.eye(3), ["Zn"], [[0, 0, 0]]),
        )
    writers.write_strain_force_in(p, blocks, reference_cell=cell)
    p.write_text(p.read_text().replace("\n/\n", "\n", 1))  # drop the terminator
    with pytest.raises(ValueError, match="terminated"):
        writers.read_strain_force_in(p)
