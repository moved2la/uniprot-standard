"""The Match Rate stage, on synthetic references and foods — and on the four foods of the author's spreadsheet.

The transcription test is the point of this file: match_rate() must reproduce the spreadsheet's
four scores (Base_Match_Rate_Formula.xlsx, 2026-09-17, sheet "Example Match Rate Calculation",
cells G156:J156) to the last digit from the same inputs, before anything else about the stage
counts (D87). The remaining tests cover rules M1-M5.

Nothing here reads a real standard or a real archive.
"""
import csv

import pytest

from pipeline import common, match


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

# --- the spreadsheet, as it is. Column A = the reference the sheet scored against (Gorissen 2018 Table 1,
#     human muscle, g per 100 g raw material); G:J = four foods, g per 100 g raw material (rows 12-19).
SHEET_REFERENCE = {"H": 2.8, "I": 3.4, "L": 6.3, "K": 6.6, "M": 1.7, "F": 3.8, "T": 2.9, "V": 4.3}
SHEET_FOODS = {
    "Rice":            {"H": 2.83, "I": 3.41, "L": 6.57, "K": 2.79, "M": 3.20, "F": 4.18, "T": 3.16, "V": 3.74},
    "Chia":            {"H": 0.531, "I": 0.801, "L": 1.371, "K": 0.97, "M": 0.588, "F": 1.016, "T": 0.709, "V": 0.95},
    "Pea":             {"H": 1.99, "I": 4.05, "L": 6.80, "K": 6.11, "M": 0.85, "F": 3.86, "T": 3.05, "V": 4.28},
    "ChickBreastWhite": {"H": 0.717, "I": 1.219, "L": 1.732, "K": 1.962, "M": 0.639, "F": 0.916, "T": 0.975, "V": 1.145},
}
# cells G156, H156, I156, J156 ("Percent Utilized") as Excel cached them
SHEET_SCORES = {"Rice": 0.4498904709748083, "Chia": 0.673823005137884,
                "Pea": 0.5130687318489835, "ChickBreastWhite": 0.8238015780989283}
SHEET_LIMITING = {"Rice": "K", "Chia": "K", "Pea": "M", "ChickBreastWhite": "F"}
EIGHT = ["T", "M", "F", "H", "K", "V", "I", "L"]     # the paper's own essential rows, in its order


def sheet_long_road(food: dict, ref: dict, scored: list) -> float:
    """Steps 3-10b of the spreadsheet, literally: scale the food to the reference's EAA total, take
    the ratio per amino acid, divide the scaled food by the smallest ratio ("lowest common
    denominator"), sum the "need to consume", and call the excess wasted. Returns 1 - wasted."""
    ref_total = sum(ref[a] for a in scored)
    food_total = sum(food[a] for a in scored)
    scaled = {a: food[a] / (food_total / ref_total) for a in scored}            # Step 3
    ratio = {a: scaled[a] / ref[a] for a in scored}                             # Step 7
    k = min(ratio.values())                                                     # Step 7, MIN
    lcd = {a: ratio[a] / k for a in scored}                                     # Step 8
    avg = sum(lcd.values()) / len(lcd)                                          # Step 8, AVERAGE
    normalised = {a: scaled[a] * avg for a in scored}                           # Step 9
    updated = {a: normalised[a] / ref[a] for a in scored}                       # Step 10
    final = {a: normalised[a] / min(updated.values()) for a in scored}          # Step 10b
    need = sum(final.values())                                                  # "Total Need to Consume"
    wasted = (need - ref_total) / need                                          # "Percent Wasted"
    return 1 - wasted                                                           # "Percent Utilized"


# --------------------------------------------------------------------------- fixtures

def _symbols(path):
    cp = common.new_ini()
    for a in common.AMINO_ACIDS:
        cp.add_section(a)
        cp[a]["trivial_name"] = SYMBOLS[a][0]
        cp[a]["three_letter"] = SYMBOLS[a][1]
    common.write_ini(cp, path, ["synthetic"])


def _fao(path):
    path.write_text("[indispensable_amino_acids]\nheadings = His, Ile, Leu, Lys, SAA, AAA, Thr, Trp, Val\n"
                    "location = synthetic\n\n[group.SAA]\nmembers = Met, Cys\n\n[group.AAA]\nmembers = Phe, Tyr\n",
                    encoding="utf-8")


def _gorissen(path, values=None, essential=EIGHT):
    values = values or {**{SYMBOLS[a][0]: SHEET_REFERENCE[a] for a in SHEET_REFERENCE},
                        "Serine": 2.3, "Glycine": 3.1, "Glutamic acid": 13.1, "Alanine": 4.1, "Tyrosine": 2.0, "Arginine": 4.4}
    body = "[values.human_muscle]\n" + "".join(f"{k} = {v}\n" for k, v in values.items())
    body += "\n[values.human_muscle.location]\nessential_rows = " + ", ".join(SYMBOLS[a][0] for a in essential) + "\n"
    path.write_text(body, encoding="utf-8")


def _standard(path, column, values):
    lines = ["# synthetic", "\t".join(["amino_acid", "three_letter", "name", column])]
    for a in common.AMINO_ACIDS:
        lines.append("\t".join([a, SYMBOLS[a][1], SYMBOLS[a][0], f"{values[a]:.4f}"]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _usda(path, foods):
    """foods: [(fdc_id, description, protein, {letter: value or None})]"""
    cols = ["fdc_id", "data_type", "release", "description", "food_category", "publication_date",
            "protein_g_per_100g", "nitrogen_g_per_100g"] + [f"{a}_g_per_100g" for a in common.AMINO_ACIDS] + \
           ["cystine_g_per_100g", "hydroxyproline_g_per_100g", "amino_acids_present", "amino_acids_absent", "min_data_points"]
    lines = ["# synthetic", "\t".join(cols)]
    for fdc, desc, protein, vals in foods:
        row = [fdc, "sr_legacy_food", "2018-04", desc, "Synthetic", "2018-04-01", f"{protein}", ""]
        row += ["" if vals.get(a) is None else f"{vals[a]}" for a in common.AMINO_ACIDS]
        row += ["", "", str(sum(v is not None for v in vals.values())), "", "3"]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _other(path, foods):
    """foods: [(label, source, protein, {letter: value})]"""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(match.OTHER_REQUIRED)
        for label, source, protein, vals in foods:
            w.writerow([label, source, "synthetic", "", protein] +
                       ["" if vals.get(a) is None else vals[a] for a in common.AMINO_ACIDS])


def _config(path, *, primary="new", references=("new", "old")):
    refs = {
        "new": "[reference.new]\nfile = outputs/standard/new.tsv\ncolumn = standard\nscored_set = nine\n"
               "label = synthetic new reference\ndecisions = test\n\n",
        "old": "[reference.old]\nconfig = config/gorissen.ini\nsection = values.human_muscle\nscored_set = eight\n"
               "label = synthetic old reference\ndecisions = test\n\n",
    }
    text = ("[meta]\nversion = test\n\n[formula]\nprimary_reference = " + primary + "\n\n"
            "[scored_set.nine]\ndecision = D86\nkind = fao_headings_reduced\nfrom_file = config/fao.ini\n"
            "group_scored_as.SAA = Met\ngroup_scored_as.AAA = Phe\n\n"
            "[scored_set.eight]\nkind = rows_named_in_file\nfrom_file = config/gorissen.ini\n"
            "from_section = values.human_muscle.location\nfrom_key = essential_rows\n\n"
            + "".join(refs[r] for r in references) +
            "[foods.usda]\nfile = outputs/usda/foods.tsv\n\n[foods.other_sources]\nfile = config/other.csv\n\n"
            "[output]\ndecimals = 4\n")
    path.write_text(text, encoding="utf-8")


NEW_STANDARD = {"A": 6.5, "C": 1.2, "D": 5.5, "E": 12.0, "F": 4.8, "G": 4.1, "H": 3.1, "I": 6.5, "K": 11.3,
                "L": 10.0, "M": 3.2, "N": 3.5, "P": 4.5, "Q": 4.5, "R": 8.2, "S": 5.7, "T": 6.1, "V": 7.1,
                "W": 1.2, "Y": 4.6}


@pytest.fixture
def match_repo(tmp_path, monkeypatch):
    """A tiny repository: two references, a USDA table of four foods, the author's other-sources file
    holding the spreadsheet's four foods. Paths are patched through monkeypatch, never assigned."""
    (tmp_path / "config").mkdir()
    (tmp_path / "outputs" / "standard").mkdir(parents=True)
    (tmp_path / "outputs" / "usda").mkdir(parents=True)
    (tmp_path / "data" / "iupac").mkdir(parents=True)
    _symbols(tmp_path / "data" / "iupac" / "symbols.ini")
    _fao(tmp_path / "config" / "fao.ini")
    _gorissen(tmp_path / "config" / "gorissen.ini")
    _standard(tmp_path / "outputs" / "standard" / "new.tsv", "standard", NEW_STANDARD)
    full = {a: 1.0 for a in common.AMINO_ACIDS}
    _usda(tmp_path / "outputs" / "usda" / "foods.tsv", [
        ("1", "Even food", 10.0, full),                                             # every ratio equal -> 100 %
        ("2", "No tryptophan", 10.0, {**full, "W": None}),                         # unscorable on the nine, scorable on the eight
        ("3", "Zero lysine", 10.0, {**full, "K": 0.0}),                            # a published zero
        ("4", "Rice", 80.0, {**{a: None for a in common.AMINO_ACIDS}, **SHEET_FOODS["Rice"]}),  # same label as an other-sources food
    ])
    _other(tmp_path / "config" / "other.csv",
           [(name, "the author's spreadsheet, 2026-09-17", 0, vals) for name, vals in SHEET_FOODS.items()])
    _config(tmp_path / "config" / "match_rate.ini")
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(common, "AMINO_ACID_SYMBOLS_INI", tmp_path / "data" / "iupac" / "symbols.ini")
    monkeypatch.setattr(match, "CONFIG_INI", tmp_path / "config" / "match_rate.ini")
    assert match.CONFIG_INI.is_relative_to(tmp_path), "the fixture did not take effect"
    return tmp_path


def _rows(text):
    body = [ln.split("\t") for ln in text.splitlines() if not ln.startswith("#")]
    return [dict(zip(body[0], r)) for r in body[1:]]


# --------------------------------------------------------------------------- the transcription (D87)

def test_the_formula_reproduces_the_spreadsheet_to_the_last_digit():
    for name, food in SHEET_FOODS.items():
        res = match.match_rate(food, SHEET_REFERENCE, EIGHT)
        assert res["score"] == pytest.approx(SHEET_SCORES[name], abs=1e-12), name
        assert res["limiting"] == [SHEET_LIMITING[name]], name


def test_the_spreadsheets_long_road_is_the_minimum_ratio():
    """Steps 8-10b (LCD, average, normalise, need to consume, wasted) reduce to 1 - min ratio (M5)."""
    for name, food in SHEET_FOODS.items():
        long = sheet_long_road(food, SHEET_REFERENCE, EIGHT)
        short = match.match_rate(food, SHEET_REFERENCE, EIGHT)["score"]
        assert long == pytest.approx(short, abs=1e-12), name
        assert long == pytest.approx(SHEET_SCORES[name], abs=1e-12), name


def test_the_stage_scores_the_spreadsheet_foods_against_the_old_reference(match_repo):
    files = match.build(_Log())
    rows = [r for r in _rows(files["match_rate_steps.tsv"]) if r["source"] == "other" and r["reference"] == "old"]
    got = {r["food_id"]: (float(r["match_rate_percent"]), r["limiting_amino_acid"]) for r in rows}
    for name in SHEET_FOODS:
        assert got[name][0] == pytest.approx(100 * SHEET_SCORES[name], abs=5e-5), name
        assert got[name][1] == SHEET_LIMITING[name], name


def test_the_steps_table_reproduces_the_spreadsheets_cells(match_repo):
    """Pea, column I of the spreadsheet: I22 = 30.99, I23 = 0.97453, I153 = 61.98, I155 = 0.48693, I156 = 0.51307."""
    files = match.build(_Log())
    pea = [r for r in _rows(files["match_rate_steps.tsv"]) if r["food_id"] == "Pea" and r["reference"] == "old"][0]
    assert float(pea["eaa_g_per_100g"]) == pytest.approx(30.99, abs=1e-4)
    assert float(pea["reference_eaa_total"]) == pytest.approx(31.8, abs=1e-4)
    assert float(pea["percent_of_reference_total"]) == pytest.approx(30.99 / 31.8, abs=1e-4)
    assert float(pea["step3_scaled_M"]) == pytest.approx(0.85 * 31.8 / 30.99, abs=1e-4)      # I34
    assert float(pea["step7_percent_of_reference_M"]) == pytest.approx(SHEET_SCORES["Pea"], abs=1e-4)
    assert float(pea["step7_min"]) == pytest.approx(SHEET_SCORES["Pea"], abs=1e-4)
    assert pea["limiting_amino_acid"] == "M"
    assert float(pea["step10b_M"]) == pytest.approx(1.7, abs=1e-4)                           # the limiting one meets the reference
    assert float(pea["step10b_total_need_to_consume"]) == pytest.approx(31.8 / SHEET_SCORES["Pea"], abs=1e-3)
    assert float(pea["percent_wasted"]) == pytest.approx(1 - SHEET_SCORES["Pea"], abs=1e-4)
    assert float(pea["percent_utilized"]) == pytest.approx(SHEET_SCORES["Pea"], abs=1e-4)
    assert float(pea["match_rate_percent"]) == pytest.approx(100 * SHEET_SCORES["Pea"], abs=5e-3)
    assert pea["W_g_per_100g"] == "" and pea["step10b_W"] == "", "the old set has no tryptophan column filled"


def test_the_per_food_table_carries_the_eaa_columns_and_not_the_dropped_ones(match_repo):
    files = match.build(_Log())
    rows = _rows(files["match_rate_per_food.tsv"])
    even = [r for r in rows if r["food_id"] == "1"][0]
    assert float(even["eaa_g_per_100g"]) == pytest.approx(9.0)                 # nine amino acids at 1.0
    assert float(even["eaa_percent_of_protein"]) == pytest.approx(90.0)        # 9 / 10 g protein
    for gone in ("reference", "reference_label", "scored_set", "scored_amino_acids", "scored_sum_g_per_100g",
                 "min_data_points", "citation"):
        assert gone not in even, gone                   # the reference is the same on every row: it is in the header
    assert "ratio_K" in even and "ratio_W" in even
    steps = [r for r in _rows(files["match_rate_steps.tsv"]) if r["food_id"] == "1" and r["reference"] == "new"][0]
    assert all(f"{a}_g_per_100g" in steps for a in "HILKMFTWV"), "the steps table carries the grams"


def test_the_per_food_table_is_one_row_per_food_on_the_primary_reference(match_repo):
    """D127: the summary table. One row per food, the primary (`new`) only — otherwise it is
    match_rate_steps.tsv without the walk. The other references live in the by-reference table."""
    files = match.build(_Log())
    rows = _rows(files["match_rate_per_food.tsv"])
    ids = [(r["source"], r["food_id"]) for r in rows]
    assert len(ids) == len(set(ids)), "a food appears more than once"
    head = [ln for ln in files["match_rate_per_food.tsv"].splitlines() if ln.startswith("#")]
    assert any("primary reference (new) only" in ln or "reference (new) only" in ln for ln in head)
    # food 2 carries no tryptophan: unscorable on the primary's nine, so it is not in the summary,
    # but it is still scored on the eight and still has its row in the by-reference table
    assert ("usda", "2") not in ids
    assert any(r["food_id"] == "2" for r in _rows(files["match_rate_by_reference.tsv"]))
    assert any(r["food_id"] == "2" and r["reference"] == "new" for r in _rows(files["foods_not_scored.tsv"]))


def test_headers_carry_provenance_and_the_reference_values_only(match_repo):
    files = match.build(_Log())
    head = [ln for ln in files["match_rate_per_food.tsv"].splitlines() if ln.startswith("#")]
    assert any("input.config" in ln and "sha256" in ln for ln in head)
    assert any(ln.startswith("# reference.old = ") and "values: T 2.9" in ln for ln in head)
    assert not any("M3:" in ln or "M4:" in ln for ln in head), "rule prose belongs in the docs, not the header"


# --------------------------------------------------------------------------- M1

def test_protein_content_and_scale_never_enter_the_score():
    food = SHEET_FOODS["Pea"]
    a = match.match_rate(food, SHEET_REFERENCE, EIGHT)["score"]
    b = match.match_rate({k: 3.7 * v for k, v in food.items()}, SHEET_REFERENCE, EIGHT)["score"]
    c = match.match_rate(food, {k: 0.01 * v for k, v in SHEET_REFERENCE.items()}, EIGHT)["score"]
    assert a == pytest.approx(b, abs=1e-12) and a == pytest.approx(c, abs=1e-12)


def test_a_food_in_the_reference_proportion_scores_100_and_every_ratio_is_one(match_repo):
    files = match.build(_Log())
    even = [r for r in _rows(files["match_rate_steps.tsv"]) if r["food_id"] == "1" and r["reference"] == "old"]
    assert even, "the even food was not scored on the old reference"
    # the old reference is not even, so 'Even food' does not score 100 there; a food equal to the reference does
    res = match.match_rate(SHEET_REFERENCE, SHEET_REFERENCE, EIGHT)
    assert res["score"] == pytest.approx(1.0, abs=1e-12)
    assert set(res["limiting"]) == set(EIGHT)
    assert all(v == pytest.approx(1.0) for v in res["ratio"].values())


# --------------------------------------------------------------------------- M2

def test_the_nine_are_the_fao_headings_with_each_group_reduced(match_repo):
    files = match.build(_Log())
    summary = files["match_summary.ini"]
    assert "amino_acids = HILKMFTWV" in summary.replace(" = ", " = ")  # order of the FAO headings
    rows = [r for r in _rows(files["match_rate_steps.tsv"]) if r["reference"] == "new"]
    assert rows and "C_g_per_100g" not in rows[0] and "Y_g_per_100g" not in rows[0], "cysteine and tyrosine are not scored on the nine (D86)"
    assert all(r["step7_percent_of_reference_M"] != "" and r["step7_percent_of_reference_F"] != "" for r in rows)


def test_the_eight_come_from_the_papers_own_essential_rows(match_repo):
    files = match.build(_Log())
    rows = [r for r in _rows(files["match_rate_steps.tsv"]) if r["reference"] == "old"]
    assert rows and all(r["step7_percent_of_reference_W"] == "" for r in rows), "the old set has no tryptophan"
    assert all(r["step7_percent_of_reference_H"] != "" for r in rows)


def test_a_group_without_a_named_member_stops_the_stage(match_repo):
    cfg = match.CONFIG_INI
    cfg.write_text(cfg.read_text(encoding="utf-8").replace("group_scored_as.AAA = Phe\n", ""), encoding="utf-8")
    with pytest.raises(SystemExit, match="D86"):
        match.build(_Log())


# --------------------------------------------------------------------------- M3

def test_a_missing_scored_amino_acid_is_listed_not_scored(match_repo):
    files = match.build(_Log())
    ns = _rows(files["foods_not_scored.tsv"])
    hit = [r for r in ns if r["food_id"] == "2" and r["reference"] == "new"]
    assert hit and "W" in hit[0]["reason"]
    scored = [r for r in _rows(files["match_rate_steps.tsv"]) if r["food_id"] == "2"]
    assert {r["reference"] for r in scored} == {"old"}, "no tryptophan is still scorable on the eight"


def test_a_published_zero_scores_zero_and_names_the_limiting_amino_acid(match_repo):
    files = match.build(_Log())
    rows = [r for r in _rows(files["match_rate_per_food.tsv"]) if r["food_id"] == "3"]
    assert rows
    for r in rows:
        assert float(r["match_rate_percent"]) == 0.0
        assert r["limiting_amino_acid"] == "K"


# --------------------------------------------------------------------------- M4

def test_no_food_is_dropped_or_preferred_across_sources(match_repo):
    files = match.build(_Log())
    rice = [r for r in _rows(files["match_rate_by_reference.tsv"]) if r["description"] == "Rice"]
    assert {r["source"] for r in rice} == {"usda", "other"}
    assert len(rice) == 2


def test_every_row_names_its_source_and_reference(match_repo):
    files = match.build(_Log())
    for name in ("match_rate_per_food.tsv", "match_rate_by_reference.tsv", "match_rate_steps.tsv"):
        for r in _rows(files[name]):
            assert r["source"] in ("usda", "other") and r["source_detail"], name
            if r["source"] == "other":
                assert r["source_detail"].startswith("the author's spreadsheet"), name
    for r in _rows(files["match_rate_steps.tsv"]):
        assert r["reference"] in ("new", "old")


def test_a_duplicate_label_in_the_other_sources_file_stops_the_stage(match_repo):
    _other(match_repo / "config" / "other.csv",
           [("Pea", "x", 0, SHEET_FOODS["Pea"]), ("Pea", "y", 0, SHEET_FOODS["Pea"])])
    with pytest.raises(SystemExit, match="twice"):
        match.build(_Log())


def test_an_other_sources_food_without_a_source_stops_the_stage(match_repo):
    _other(match_repo / "config" / "other.csv", [("Pea", "", 0, SHEET_FOODS["Pea"])])
    with pytest.raises(SystemExit, match="source is empty"):
        match.build(_Log())


def test_a_header_only_other_sources_file_is_allowed(match_repo):
    _other(match_repo / "config" / "other.csv", [])
    files = match.build(_Log())
    assert all(r["source"] == "usda" for r in _rows(files["match_rate_by_reference.tsv"]))


# --------------------------------------------------------------------------- outputs

def test_build_writes_every_table(match_repo):
    files = match.build(_Log())
    for name in ("match_rate_per_food.tsv", "match_rate_by_reference.tsv", "match_rate_steps.tsv",
                 "new/match_rate_ranked_new.tsv", "old/match_rate_ranked_old.tsv", "limiting_amino_acid_counts.tsv",
                 "foods_not_scored.tsv", "match_summary.ini"):
        assert name in files, name


def test_the_by_reference_table_puts_old_beside_new_with_the_difference(match_repo):
    """D125: the difference is the OTHER reference minus the primary, so a positive number means
    that reference scores the food closer to 100 than the primary does."""
    files = match.build(_Log())
    rows = _rows(files["match_rate_by_reference.tsv"])
    cols = list(rows[0].keys())
    assert "match_new_percent" in cols and "match_old_percent" in cols and "difference_pp_old_minus_new" in cols
    assert "difference_pp_new_minus_old" not in cols              # `new` is the primary; it is not subtracted from itself
    even = [r for r in rows if r["food_id"] == "1"][0]          # scorable on both sets
    assert float(even["difference_pp_old_minus_new"]) == pytest.approx(
        float(even["match_old_percent"]) - float(even["match_new_percent"]), abs=2e-4)
    pea = [r for r in rows if r["food_id"] == "Pea"][0]          # the sheet's foods carry no tryptophan
    assert pea["match_old_percent"] and pea["match_new_percent"] == "" and pea["difference_pp_old_minus_new"] == ""


def test_the_steps_table_puts_the_match_rate_beside_the_food_name(match_repo):
    """D126: the headline score is column 6, not the last column; the spreadsheet's walk still ends
    at row 156 (percent_utilized), and match_rate_percent is that number x 100."""
    files = match.build(_Log())
    rows = _rows(files["match_rate_steps.tsv"])
    cols = list(rows[0].keys())
    assert cols[:7] == ["source", "source_detail", "food_id", "description", "reference",
                        "match_rate_percent", "protein_g_per_100g"]
    assert cols[-1] == "percent_utilized"
    r = rows[0]
    # both columns are written at four decimals, so the x 100 comparison carries about 5e-3 of rounding
    assert float(r["match_rate_percent"]) == pytest.approx(100 * float(r["percent_utilized"]), abs=1e-2)


def test_the_ranked_table_is_best_first(match_repo):
    files = match.build(_Log())
    vals = [float(r["match_rate_percent"]) for r in _rows(files["old/match_rate_ranked_old.tsv"])]
    assert vals == sorted(vals, reverse=True)


def test_limiting_counts_add_up(match_repo):
    files = match.build(_Log())
    rows = [r for r in _rows(files["limiting_amino_acid_counts.tsv"]) if r["reference"] == "old"]
    scored = int(rows[0]["foods_scored"])
    assert sum(int(r["foods_limited"]) for r in rows) >= scored     # ties count more than once, never fewer
