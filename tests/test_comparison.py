"""The comparison stage, on a synthetic paper and a synthetic standard.

The rules under test are X1-X5: one shared denominator, the set decided by the measured side's
stated method, no verdict, the transcription checked against the paper's printed sums, and a
reported zero carried as a zero.
"""
import pytest

from pipeline import common, comparison as cmp


class _Log:
    def info(self, *a):
        pass

    def warning(self, *a):
        pass

    def error(self, *a):
        pass


SYMBOLS = {"A": ("Alanine", "Ala"), "C": ("Cysteine", "Cys"), "D": ("Aspartic acid", "Asp"),
           "E": ("Glutamic acid", "Glu"), "F": ("Phenylalanine", "Phe"), "G": ("Glycine", "Gly"),
           "H": ("Histidine", "His"), "I": ("Isoleucine", "Ile"), "K": ("Lysine", "Lys"),
           "L": ("Leucine", "Leu"), "M": ("Methionine", "Met"), "N": ("Asparagine", "Asn"),
           "P": ("Proline", "Pro"), "Q": ("Glutamine", "Gln"), "R": ("Arginine", "Arg"),
           "S": ("Serine", "Ser"), "T": ("Threonine", "Thr"), "V": ("Valine", "Val"),
           "W": ("Tryptophan", "Trp"), "Y": ("Tyrosine", "Tyr")}

# a synthetic "paper": eight essential rows summing to 20.0, eight non-essential to 20.0
# Deliberately not a plausible amino acid profile: if these values ever turn up in a real config
# they must be recognisable at a glance as test data. Essential rows sum to 20.0, the rest to 20.0.
MEASURED = {"Threonine": 1.0, "Methionine": 1.0, "Phenylalanine": 1.0, "Histidine": 1.0,
            "Lysine": 4.0, "Valine": 4.0, "Isoleucine": 4.0, "Leucine": 4.0,
            "Serine": 5.0, "Glycine": 5.0, "Glutamic acid": 5.0, "Alanine": 5.0}


def _symbols(path):
    cp = common.new_ini()
    for letter, (name, three) in SYMBOLS.items():
        cp.add_section(letter)
        cp[letter]["trivial_name"] = name
        cp[letter]["three_letter"] = three
    common.write_ini(cp, path, ["synthetic"])
    return path


def _standard(path, col, values):
    lines = ["# synthetic", "\t".join(["amino_acid", "three_letter", "name", col])]
    for a in common.AMINO_ACIDS:
        lines.append("\t".join([a, SYMBOLS[a][1], SYMBOLS[a][0], f"{values[a]:.4f}"]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _config(path, measured=MEASURED, eaa_sum=20.0, neaa_sum=20.0, not_measured=True):
    lines = ["[meta]", "source_id = gorissen_2018", "decisions = D84", "",
             "[sample]", "tissue = synthetic", "n = 1", "",
             "[method]", "hydrolysis = synthetic", "",
             "[protein_content_of_the_muscle_sample]", "value = 84", "", "[values.human_muscle]"]
    lines += [f"{k} = {v}" for k, v in measured.items()]
    lines += ["", "[values.human_muscle.location]",
              "essential_rows = Threonine, Methionine, Phenylalanine, Histidine, Lysine, Valine, Isoleucine, Leucine",
              f"sums_as_printed = EAA {eaa_sum:.1f}, non-EAA {neaa_sum:.1f}", "",
              "[hydrolysis_conversions]", "Glutamic acid = Glutamic acid + Glutamine", ""]
    if not_measured:
        lines += ["[not_measured]", "Tryptophan = destroyed", "Aspartic acid = not reported",
                  "Asparagine = converts to aspartic acid",
                  "Proline = printed as zero; not a measurement",
                  "Cysteine = printed as zero; not a measurement",
                  "Tyrosine = not reported", "Arginine = not reported", ""]
    lines += ["[caveats]", "single_factor = synthetic", "", "[mass_balance]", "compare_against = synthetic", ""]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


@pytest.fixture
def comparison_repo(tmp_path, monkeypatch):
    out = tmp_path / "outputs" / "comparison"
    out.mkdir(parents=True)
    flat = {a: 5.0 for a in common.AMINO_ACIDS}                  # every amino acid equal
    hi_h = dict(flat, H=10.0)                                     # ... except histidine, doubled
    monkeypatch.setattr(cmp, "OUT_DIR", out)
    monkeypatch.setattr(cmp, "INPUTS", {
        "comparison": _config(tmp_path / "gorissen_2018_comparison.ini"),
        "sources": _sources(tmp_path / "literature_sources.ini"),
        "symbols": _symbols(tmp_path / "amino_acid_symbols.ini"),
        "standard": _standard(tmp_path / "std.tsv", "standard", flat),
        "standard_npmp": _standard(tmp_path / "std_npmp.tsv",
                                   "standard_with_non_protein_metabolite_pools", hi_h),
        "masses": _masses(tmp_path / "amino_acid_masses.ini"),
        "pools": _pools(tmp_path / "non_protein_metabolite_pools.ini"),
    })
    monkeypatch.setattr(cmp, "MANIFEST", _manifest(tmp_path, monkeypatch))
    assert cmp.INPUTS["comparison"].is_relative_to(tmp_path), \
        "the fixture did not take effect; a test must never write to the real config"
    return tmp_path


def _masses(path):
    """Free and residue masses. Real values: the factor the mass balance computes is a property of
    chemistry, and a synthetic mass would make the test assert nothing."""
    free = {"G": 75.07, "A": 89.09, "S": 105.09, "P": 115.13, "V": 117.15, "T": 119.12,
            "C": 121.16, "L": 131.17, "I": 131.17, "N": 132.12, "D": 133.10, "Q": 146.15,
            "K": 146.19, "E": 147.13, "M": 149.21, "H": 155.15, "F": 165.19, "R": 174.20,
            "Y": 181.19, "W": 204.23}
    cp = common.new_ini()
    for a in common.AMINO_ACIDS:
        cp.add_section(a)
        cp[a]["free_mass"] = f"{free[a]}"
        cp[a]["residue_mass"] = f"{free[a] - 18.02:.2f}"
        cp[a]["three_letter"] = SYMBOLS[a][1]
        cp[a]["trivial_name"] = SYMBOLS[a][0]
    common.write_ini(cp, path, ["synthetic"])
    return path


def _pools(path):
    """The second protein content the mass balance tests, as the pools config states it."""
    path.write_text("[protein_content]\nvalue = 155\nsource_id = mingrone_2001\n\n"
                    "[water_content]\nvalue = 763\nsource_id = mingrone_2001\n", encoding="utf-8")
    return path


def _sources(path):
    path.write_text("[source.gorissen_2018]\nid_short = P10\ncitation = Synthetic et al.\ndoi = 10.0000/x\n",
                    encoding="utf-8")
    return path


def _manifest(tmp_path, monkeypatch):
    """A one-file manifest, with the repo root pointed at the temp directory for the duration.

    REPO_ROOT is patched through monkeypatch, never assigned directly: a direct assignment is not
    undone when the test ends, and every later test in the session then resolves its own paths
    against this temp directory."""
    import hashlib
    pdf = tmp_path / "paper.pdf"
    pdf.write_bytes(b"%PDF-1.4 synthetic\n")
    man = tmp_path / "manifest.ini"
    rel = pdf.relative_to(tmp_path).as_posix()
    man.write_text(f"[gorissen_2018.file.1]\nsource_id = gorissen_2018\npath = {rel}\n"
                   f"sha256 = {hashlib.sha256(pdf.read_bytes()).hexdigest()}\n", encoding="utf-8")
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    return man


def _rows(text):
    body = [ln.split("\t") for ln in text.splitlines() if not ln.startswith("#")]
    return [dict(zip(body[0], r)) for r in body[1:]]


def test_build_writes_every_table(comparison_repo):
    files = cmp.build(_Log())
    for name in ("gorissen_2018_human_muscle.tsv", "calculated_vs_gorissen_2018.tsv",
                 "comparison_summary.ini"):
        assert name in files, name


def test_both_sides_are_percentages_of_the_same_set(comparison_repo):
    """X1: each column sums to 100 over the comparison rows.

    The tolerance is the printing, not the arithmetic: the table carries four decimals, so
    sixteen rows can land up to ~1e-3 from 100 when summed as printed."""
    rows = [r for r in _rows(cmp.build(_Log())["calculated_vs_gorissen_2018.tsv"])
            if r["status"] == "compared"]
    for col in ("gorissen_percent_of_measured", "calculated_protein_only_percent",
                "calculated_with_non_protein_metabolite_pools_percent"):
        assert sum(float(r[col]) for r in rows) == pytest.approx(100.0, abs=2e-3), col


def test_converted_amino_acids_are_summed_on_both_sides(comparison_repo):
    """X2: glutamine converts to glutamic acid, so E and Q are one row."""
    rows = {r["row"]: r for r in _rows(cmp.build(_Log())["calculated_vs_gorissen_2018.tsv"])}
    assert "E+Q" in rows and "E" not in rows and "Q" not in rows
    # the synthetic standard is flat at 5.0 across twenty; twelve rows cover thirteen letters,
    # so the comparison set totals 65 and E+Q, being two of them, is 10 of it
    assert float(rows["E+Q"]["calculated_protein_only_percent"]) == pytest.approx(10 / 65 * 100, abs=1e-3)


def test_unmeasured_amino_acids_are_excluded_from_both_sides_and_named(comparison_repo):
    """X2: what the method cannot report leaves the denominator on BOTH sides."""
    rows = {r["row"]: r for r in _rows(cmp.build(_Log())["calculated_vs_gorissen_2018.tsv"])}
    for letter in ("D", "N", "W", "P", "C"):
        assert rows[letter]["status"] == "not measured"
        assert rows[letter]["caveat"]
        assert rows[letter]["gorissen_percent_of_measured"] == ""


def test_an_undeclared_zero_stops_the_stage(comparison_repo):
    """X5: a zero is never compared as a measurement unless it is declared.

    This is the rule that was wrong in the first cut. A zero left in the measured values sits in
    one side's denominator and not the other's, which biases every other row."""
    bad = dict(MEASURED, Alanine=0.0, Glycine=10.0)          # keep the non-essential sum at 20.0
    _config(comparison_repo / "gorissen_2018_comparison.ini", measured=bad)
    with pytest.raises(SystemExit) as e:
        cmp.build(_Log())
    assert "X5" in str(e.value) and "Alanine" in str(e.value)


def test_a_zero_that_is_declared_as_a_measurement_is_compared(comparison_repo):
    """X5's other half: a zero the source supports is carried, with no ratio."""
    bad = dict(MEASURED, Alanine=0.0, Glycine=10.0)
    path = _config(comparison_repo / "gorissen_2018_comparison.ini", measured=bad)
    path.write_text(path.read_text(encoding="utf-8")
                    + "\n[zero_is_a_measurement]\nAlanine = the synthetic source states this is a real zero\n",
                    encoding="utf-8")
    rows = {r["row"]: r for r in _rows(cmp.build(_Log())["calculated_vs_gorissen_2018.tsv"])}
    assert rows["A"]["status"] == "compared"
    assert float(rows["A"]["gorissen_percent_of_measured"]) == 0.0
    assert rows["A"]["ratio_calculated_over_gorissen_protein_only"] == ""


def test_the_mass_balance_is_computed_for_every_protein_content_offered(comparison_repo):
    """X6: expected yield, what is measurable, and what the source reports — no verdict."""
    rows = _rows(cmp.build(_Log())["gorissen_2018_mass_balance.tsv"])
    assert rows, "no mass balance rows"
    for r in rows:
        assert float(r["free_per_protein_factor"]) > 1.0, "hydrolysis adds mass; the factor must exceed 1"
        assert float(r["expected_free_amino_acids_g"]) > float(r["protein_content_percent"])
        assert r["where_the_protein_content_comes_from"]


def test_no_verdict_is_written(comparison_repo):
    """X3: difference and ratio only — no error, no fit, no 'correct'."""
    files = cmp.build(_Log())
    text = files["calculated_vs_gorissen_2018.tsv"]
    header = [ln for ln in text.splitlines() if ln.startswith("#")]
    cols = [ln for ln in text.splitlines() if not ln.startswith("#")][0].lower()
    for word in ("error", "accuracy", "correct", "true_value", "bias"):
        assert word not in cols, f"{word} in the column names is a verdict (X3)"
    assert any("neither side is treated as the truth" in ln.lower() for ln in header)


def test_a_transcription_that_misses_the_printed_sum_stops_the_stage(comparison_repo):
    """X4: the paper's own totals are the check on the transcription."""
    bad = dict(MEASURED, Lysine=9.0)                 # essential sum now 25.0, not 20.0
    _config(comparison_repo / "gorissen_2018_comparison.ini", measured=bad)
    with pytest.raises(SystemExit) as e:
        cmp.build(_Log())
    assert "X4" in str(e.value) and "25.0" in str(e.value)


def test_an_amino_acid_that_is_neither_compared_nor_excluded_stops_the_stage(comparison_repo):
    """Every one of the twenty must be accounted for, one way or the other."""
    _config(comparison_repo / "gorissen_2018_comparison.ini", not_measured=False)
    with pytest.raises(SystemExit) as e:
        cmp.build(_Log())
    assert "neither a measured row nor" in str(e.value)


def test_the_green_light_refuses_a_source_that_is_not_hashed(comparison_repo):
    """D78: no hashed PDF, no comparison."""
    cmp.MANIFEST.write_text("[other.file.1]\nsource_id = someone_else\n", encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        cmp.build(_Log())
    assert "gorissen_2018" in str(e.value)


def test_the_config_names_the_rows_the_paper_names(comparison_repo):
    """The join is by IUPAC trivial name; a row label that is not one stops the stage."""
    bad = {("Glutamate" if k == "Glutamic acid" else k): v for k, v in MEASURED.items()}
    _config(comparison_repo / "gorissen_2018_comparison.ini", measured=bad)
    with pytest.raises(SystemExit) as e:
        cmp.build(_Log())
    assert "IUPAC trivial name" in str(e.value)
    

def test_the_printed_sum_check_reads_the_papers_own_grouping_from_config(comparison_repo):
    """X4 groups the rows as [values.human_muscle.location] essential_rows says; without it the stage stops,
    because the alternative is naming amino acids in code."""
    cfg = cmp.INPUTS["comparison"]
    text = cfg.read_text(encoding="utf-8")
    cfg.write_text("\n".join(ln for ln in text.splitlines() if not ln.startswith("essential_rows")), encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        cmp.build(_Log())
    assert "essential_rows" in str(e.value)
