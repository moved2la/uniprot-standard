"""composition.py on synthetic sequences and the synthetic mass table, offline.

What these tests check, in plain terms:
- counting: a hand-computable sequence gives the expected counts, sums, and
  molecular weight (with the synthetic masses: free = 100 + index, water = 10);
- the water identity: free_mass_sum == mw + (n - 1) * water for every sequence tried;
- both fraction vectors sum to 1;
- a non-standard letter is refused;
- segment slicing: a false-flagged N-terminal run is dropped from the master set
  and the processing delta row reports exactly what was removed;
- the wide table has one row per accession per segment set.
"""
import logging
import pytest
from conftest import FIXTURES

from pipeline import common
from pipeline import composition as comp

AA = common.AMINO_ACIDS
MASSES = comp.load_masses(FIXTURES / "masses_synthetic.ini")
IDX = {a: i for i, a in enumerate(AA)}


def free(a):  # synthetic free mass
    return 100 + IDX[a]


def test_compose_hand_computable():
    seq = "GGGAAAK"
    c = comp.compose(seq, MASSES)
    assert c["n_residues"] == 7
    assert c["count"]["G"] == 3 and c["count"]["A"] == 3 and c["count"]["K"] == 1 and c["count"]["W"] == 0
    expected_free = 3 * free("G") + 3 * free("A") + free("K")
    assert c["free_mass_sum"] == pytest.approx(expected_free)
    assert c["residue_mass_sum"] == pytest.approx(expected_free - 7 * 10)
    assert c["mw"] == pytest.approx(expected_free - 7 * 10 + 10)
    assert c["free_frac"]["G"] == pytest.approx(3 * free("G") / expected_free)


@pytest.mark.parametrize("seq", ["A", "GGGAAAK", AA, AA * 3 + "WWW", "M" + AA[::-1]])
def test_water_identity_and_fraction_sums(seq):
    c = comp.compose(seq, MASSES)
    assert comp.water_identity_error(c, MASSES["water"]) < 1e-9
    assert sum(c["free_frac"].values()) == pytest.approx(1.0)
    assert sum(c["residue_frac"].values()) == pytest.approx(1.0)


def test_non_standard_letter_refused():
    with pytest.raises(ValueError):
        comp.compose("GGXAA", MASSES)
    with pytest.raises(ValueError):
        comp.compose("GGUAA", MASSES)


def _segments(acc, removed_end, length):
    seg = common.new_ini()
    seg[f"{acc}.n_terminal_removed"] = {"start": "1", "end": str(removed_end), "in_master_molecule": "false"}
    seg[f"{acc}.master"] = {"start": str(removed_end + 1), "end": str(length), "in_master_molecule": "true"}
    return seg


def test_segment_slicing():
    seq = "MRRGGGAAAK"           # 10 residues; drop the first 3
    seg = _segments("X1", 3, len(seq))
    master = comp.slice_sequence(seq, comp.segment_ranges(seg, "X1", master_only=True))
    full = comp.slice_sequence(seq, comp.segment_ranges(seg, "X1", master_only=False))
    assert master == "GGGAAAK" and full == seq
    cm, cf = comp.compose(master, MASSES), comp.compose(full, MASSES)
    assert cf["n_residues"] - cm["n_residues"] == 3
    assert cf["count"]["R"] - cm["count"]["R"] == 2 and cf["count"]["M"] - cm["count"]["M"] == 1


def test_build_rows_and_processing_delta():
    accessions = common.new_ini()
    accessions["X1"] = {"gene": "FIXA", "tier": "1", "flag_open": "false"}
    accessions["X2"] = {"gene": "FIXB", "tier": "1;2", "flag_open": "false"}
    sequences = common.new_ini()
    sequences["X1"] = {"sequence": "MRRGGGAAAK"}
    sequences["X2"] = {"sequence": "GGGAAAK"}
    seg = _segments("X1", 3, 10)
    seg["X2.master"] = {"start": "1", "end": "7", "in_master_molecule": "true"}
    log = logging.getLogger("test")
    all_rows, delta_rows, bound_rows, summary = comp.build(accessions, sequences, seg, MASSES, log)

    assert len(all_rows) == 4 and [r[3] for r in all_rows] == ["master", "metabolic", "master", "metabolic"]
    assert summary["rows_all"] == "4" and summary["entries_master_equals_full"] == "1"
    x1 = {r[0]: r for r in delta_rows}["X1"]
    n_full, n_master, n_removed = x1[3], x1[4], x1[5]
    assert (n_full, n_master, n_removed) == (10, 7, 3)
    mw_full, mw_master = float(x1[6]), float(x1[7])
    removed = float(x1[8])
    assert removed == pytest.approx((mw_full - mw_master) / mw_full)
    # X2 is unprocessed: zero delta everywhere
    x2 = {r[0]: r for r in delta_rows}["X2"]
    assert float(x2[8]) == 0.0 and float(x2[9]) == 0.0
    assert len(bound_rows) == 20
    assert float(summary["water_identity_max_abs_error_g_per_mol"]) < 1e-9


def test_open_flag_is_skipped_not_computed():
    accessions = common.new_ini()
    accessions["X1"] = {"gene": "FIXA", "tier": "1", "flag_open": "true"}
    sequences = common.new_ini()
    sequences["X1"] = {"sequence": "GGGAAAK"}
    seg = common.new_ini()
    all_rows, delta_rows, _, summary = comp.build(accessions, sequences, seg, MASSES, logging.getLogger("t"))
    assert all_rows == [] and delta_rows == [] and summary["accessions_skipped_open_flag"] == "X1"
