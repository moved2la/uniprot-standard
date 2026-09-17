"""The USDA stage, against a synthetic FoodData Central archive built in a temp directory.

Nothing here touches a real download: the archive is written row by row so that each rule can be
checked on its own — the membership filter (U2), the trivial-name join (U3), the unit check and the
untouched values (U4), the stated-not-decided completeness (U5), and no cross-archive dedup (U6).
"""
import csv
import io
import zipfile

import pytest

from pipeline import common, usda


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

# nutrient ids for the synthetic archive: the twenty, then the extras
NUT = {letter: str(2000 + i) for i, letter in enumerate(common.AMINO_ACIDS)}
NUT_PROTEIN, NUT_NITROGEN, NUT_CYSTINE, NUT_SODIUM = "1003", "1002", "1216", "1093"


def _csv(rows, cols):
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=cols)
    w.writeheader()
    for r in rows:
        w.writerow(r)
    return buf.getvalue()


def _archive(path, data_type="foundation_food", *, membership_file=True, sodium_unit="MG",
             foods=None, nutrient_rows=None, folder="FoodData_Central_test"):
    """A minimal FoodData Central CSV archive: food, nutrient, food_nutrient, membership, category."""
    foods = foods if foods is not None else [
        {"fdc_id": "1", "data_type": data_type, "description": "Egg, whole", "food_category_id": "1",
         "publication_date": "2026-04-01"},
        {"fdc_id": "2", "data_type": data_type, "description": "Butter", "food_category_id": "1",
         "publication_date": "2026-04-01"},
        {"fdc_id": "3", "data_type": data_type, "description": "Egg, whole", "food_category_id": "1",
         "publication_date": "2019-04-01"},          # a superseded earlier version: not in the membership file
    ]
    nutrients = ([{"id": NUT[a], "name": SYMBOLS[a][0], "unit_name": "G", "nutrient_nbr": "", "rank": ""}
                  for a in common.AMINO_ACIDS]
                 + [{"id": NUT_PROTEIN, "name": "Protein", "unit_name": "G", "nutrient_nbr": "", "rank": ""},
                    {"id": NUT_NITROGEN, "name": "Nitrogen", "unit_name": "G", "nutrient_nbr": "", "rank": ""},
                    {"id": NUT_CYSTINE, "name": "Cystine", "unit_name": "G", "nutrient_nbr": "", "rank": ""},
                    {"id": NUT_SODIUM, "name": "Sodium, Na", "unit_name": sodium_unit, "nutrient_nbr": "", "rank": ""}])
    if nutrient_rows is None:
        nutrient_rows = []
        # food 1 carries every amino acid but N and Q, plus protein and sodium
        for a in common.AMINO_ACIDS:
            if a in "NQ":
                continue
            nutrient_rows.append({"id": "", "fdc_id": "1", "nutrient_id": NUT[a], "amount": "1.5",
                                  "data_points": "5", "min": "1.4", "max": "1.6", "median": "1.5"})
        nutrient_rows += [
            {"id": "", "fdc_id": "1", "nutrient_id": NUT_PROTEIN, "amount": "12.5", "data_points": "",
             "min": "", "max": "", "median": ""},
            {"id": "", "fdc_id": "1", "nutrient_id": NUT_SODIUM, "amount": "140", "data_points": "",
             "min": "", "max": "", "median": ""},
            # food 2 carries only protein: no amino acid row at all
            {"id": "", "fdc_id": "2", "nutrient_id": NUT_PROTEIN, "amount": "0.85", "data_points": "",
             "min": "", "max": "", "median": ""},
            # food 3 is superseded: its amino acids must not appear
            {"id": "", "fdc_id": "3", "nutrient_id": NUT["H"], "amount": "99.0", "data_points": "",
             "min": "", "max": "", "median": ""},
        ]
    members = [{"fdc_id": "1", "NDB_number": "", "footnote": ""},
               {"fdc_id": "2", "NDB_number": "", "footnote": ""}]
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr(f"{folder}/food.csv", _csv(foods, ["fdc_id", "data_type", "description",
                                                       "food_category_id", "publication_date"]))
        zf.writestr(f"{folder}/nutrient.csv", _csv(nutrients, ["id", "name", "unit_name", "nutrient_nbr", "rank"]))
        zf.writestr(f"{folder}/food_nutrient.csv",
                    _csv(nutrient_rows, ["id", "fdc_id", "nutrient_id", "amount", "data_points", "min", "max", "median"]))
        zf.writestr(f"{folder}/food_category.csv",
                    _csv([{"id": "1", "code": "100", "description": "Test Category"}], ["id", "code", "description"]))
        if membership_file:
            zf.writestr(f"{folder}/{data_type}.csv", _csv(members, ["fdc_id", "NDB_number", "footnote"]))
    return path


def _symbols_ini(path):
    cp = common.new_ini()
    for letter, (name, three) in SYMBOLS.items():
        cp.add_section(letter)
        cp[letter]["trivial_name"] = name
        cp[letter]["three_letter"] = three
    common.write_ini(cp, path, ["synthetic"])
    return path


def _config_ini(path, data_types=("foundation_food",)):
    lines = ["[meta]", "decisions = D83", "", "[carried_nutrients]", "protein        = Protein",
             "nitrogen       = Nitrogen", "cystine        = Cystine", "",
             "[units]", "amino_acids = G", "protein = G", "nitrogen = G", ""]
    for i, dt in enumerate(data_types, 1):
        lines += [f"[archive.{dt}]", f"data_type      = {dt}", "release        = 2026-04",
                  f"published_name = Test archive {dt}", "url            = https://example.invalid",
                  "obtained       = manual", f"precedence     = {i}", ""]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


@pytest.fixture
def synthetic(tmp_path, monkeypatch):
    archives, out = tmp_path / "data" / "usda", tmp_path / "outputs" / "usda"
    archives.mkdir(parents=True)
    out.mkdir(parents=True)
    _archive(archives / "whatever_the_download_was_called.zip")
    monkeypatch.setattr(usda, "ARCHIVE_DIR", archives)
    monkeypatch.setattr(usda, "OUT_DIR", out)
    monkeypatch.setattr(usda, "MANIFEST", archives / "manifest.ini")
    monkeypatch.setattr(usda, "INPUTS", {"usda_config": _config_ini(tmp_path / "usda_food_data.ini"),
                                         "symbols": _symbols_ini(tmp_path / "amino_acid_symbols.ini")})
    return {"archives": archives, "out": out, "root": tmp_path}


def _rows(text):
    body = [ln.split("\t") for ln in text.splitlines() if not ln.startswith("#")]
    return body[0], [dict(zip(body[0], r)) for r in body[1:]]


def test_build_writes_every_table(synthetic):
    files, manifest = usda.build(_Log())
    for name in (usda.MAIN_TABLE, "foods_without_amino_acids.tsv", "usda_nutrient_map.tsv",
                 "usda_archive_summary.ini"):
        assert name in files, name
    assert "sha256" in manifest and "foundation_food" in manifest


def test_membership_file_drops_the_superseded_version(synthetic):
    """U2: food 3 is an earlier version of the same food and is not in foundation_food.csv."""
    files, _ = usda.build(_Log())
    _, rows = _rows(files[usda.MAIN_TABLE])
    assert [r["fdc_id"] for r in rows] == ["1"]
    assert all(r["H_g_per_100g"] != "99.0" for r in rows)


def test_falls_back_to_the_data_type_column_without_a_membership_file(synthetic):
    """U2: no <data_type>.csv — every food.csv row of that data type counts, and the log says so."""
    for p in synthetic["archives"].glob("*.zip"):
        p.unlink()
    _archive(synthetic["archives"] / "no_membership.zip", membership_file=False)
    files, _ = usda.build(_Log())
    _, rows = _rows(files[usda.MAIN_TABLE])
    assert sorted(r["fdc_id"] for r in rows) == ["1", "3"]      # the superseded row is now carried


def test_the_twenty_are_joined_by_trivial_name_not_by_id(synthetic):
    """U3: nothing is mapped by hand; the map is the audit trail."""
    files, _ = usda.build(_Log())
    _, rows = _rows(files["usda_nutrient_map.tsv"])
    carried = {r["carried_as"]: r for r in rows if r["outcome"] == "carried"}
    for letter in common.AMINO_ACIDS:
        assert carried[letter]["nutrient_name"] == SYMBOLS[letter][0]
        assert carried[letter]["nutrient_id"] == NUT[letter]
    assert carried["cystine"]["nutrient_name"] == "Cystine"


def test_values_are_carried_as_published(synthetic):
    """U4: the archive's amount, per 100 g of food, unconverted."""
    files, _ = usda.build(_Log())
    _, rows = _rows(files[usda.MAIN_TABLE])
    r = rows[0]
    assert r["H_g_per_100g"] == "1.5"
    assert r["protein_g_per_100g"] == "12.5"
    assert r["min_data_points"] == "5"


def test_a_nutrient_in_the_wrong_unit_is_flagged_and_not_carried(synthetic):
    """U4: a unit the config does not expect is a flag, never a silent conversion."""
    for p in synthetic["archives"].glob("*.zip"):
        p.unlink()
    # protein declared in milligrams: the unit check must catch it
    _archive(synthetic["archives"] / "bad_unit.zip")
    with zipfile.ZipFile(synthetic["archives"] / "bad_unit.zip") as zf:
        names = zf.namelist()
        data = {n: zf.read(n) for n in names}
    data = {n: (v.replace(b'"Protein","G"', b'"Protein","MG"').replace(b"Protein,G", b"Protein,MG")
                if n.endswith("nutrient.csv") else v) for n, v in data.items()}
    with zipfile.ZipFile(synthetic["archives"] / "bad_unit.zip", "w") as zf:
        for n, v in data.items():
            zf.writestr(n, v)
    files, _ = usda.build(_Log())
    summary = files["usda_archive_summary.ini"]
    assert "flag" in summary and "MG" in summary
    _, rows = _rows(files[usda.MAIN_TABLE])
    assert rows[0]["protein_g_per_100g"] == ""


def test_absent_amino_acids_are_stated_not_decided(synthetic):
    """U5: the row says what it lacks; no food is dropped for being incomplete."""
    files, _ = usda.build(_Log())
    _, rows = _rows(files[usda.MAIN_TABLE])
    assert rows[0]["amino_acids_present"] == "18"
    assert rows[0]["amino_acids_absent"] == "NQ"


def test_a_food_with_no_amino_acid_row_is_listed_not_silently_dropped(synthetic):
    files, _ = usda.build(_Log())
    _, rows = _rows(files["foods_without_amino_acids.tsv"])
    assert [r["fdc_id"] for r in rows] == ["2"]
    assert rows[0]["description"] == "Butter"


def test_two_archives_both_land_with_their_data_type(synthetic):
    """U6: the same food in two archives is two rows; nothing is deduped here."""
    _archive(synthetic["archives"] / "legacy.zip", data_type="sr_legacy_food")
    _config_ini(usda.INPUTS["usda_config"], data_types=("foundation_food", "sr_legacy_food"))
    files, _ = usda.build(_Log())
    _, rows = _rows(files[usda.MAIN_TABLE])
    assert sorted(r["data_type"] for r in rows) == ["foundation_food", "sr_legacy_food"]


def test_an_archive_is_identified_by_its_contents_not_its_file_name(synthetic):
    """U1: the zip may be called anything."""
    files, _ = usda.build(_Log())
    _, rows = _rows(files[usda.MAIN_TABLE])
    assert rows[0]["data_type"] == "foundation_food"
    assert "whatever_the_download_was_called.zip" in files["usda_archive_summary.ini"]


def test_a_missing_archive_stops_the_stage_and_names_the_data_type(synthetic):
    _config_ini(usda.INPUTS["usda_config"], data_types=("foundation_food", "sr_legacy_food"))
    with pytest.raises(SystemExit) as e:
        usda.build(_Log())
    assert "sr_legacy_food" in str(e.value) and "[STOP]" in str(e.value)


def test_two_archives_of_the_same_data_type_stop_the_stage(synthetic):
    _archive(synthetic["archives"] / "a_second_copy.zip")
    with pytest.raises(SystemExit) as e:
        usda.build(_Log())
    assert "two archives" in str(e.value)


def test_an_unfilled_release_stops_the_stage(synthetic):
    p = usda.INPUTS["usda_config"]
    p.write_text(p.read_text(encoding="utf-8").replace("release        = 2026-04", "release        = ___"),
                 encoding="utf-8")
    with pytest.raises(SystemExit) as e:
        usda.build(_Log())
    assert "release" in str(e.value)


def test_the_manifest_keeps_the_first_seen_time_of_an_unchanged_archive(synthetic):
    _, first = usda.build(_Log())
    usda.MANIFEST.write_text(first, encoding="utf-8")
    _, second = usda.build(_Log())
    line = [ln for ln in first.splitlines() if ln.startswith("first_seen")][0]
    assert line in second
