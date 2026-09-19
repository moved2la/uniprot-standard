"""mass_fractions_pools on a synthetic category: two pools (one with per-replicate shares, one with a
log10 abundance and runs), a shared row, a summed entry, an R5-excluded entry, immunoglobulin entries
for the D113 bound, a Hortin-shaped cross-check table (xls reader monkeypatched) and a deep dataset.
Every value is invented; the arithmetic is hand-checkable."""
import logging
import math

from openpyxl import Workbook

from pipeline import common, mass_fractions_pools as mfp

AA = list(common.AMINO_ACIDS)


def read(path):
    """The stage's tables carry '# ' provenance lines above the header."""
    lines = [l for l in path.read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    hdr = lines[0].split("\t")
    return [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l]


def _write_xlsx(path, sheet, header, rows):
    wb = Workbook(); ws = wb.active; ws.title = sheet
    ws.append(header)
    for r in rows:
        ws.append(list(r))
    wb.save(path)


def _composition_row(acc, gene, tier, n, mw, frac):
    row = [acc, gene, tier, "master", n, f"{mw:.4f}", f"{mw:.4f}", f"{mw:.4f}"] + ["0"] * 20
    row += [f"{frac[a]:.6f}" for a in AA] + [f"{frac[a]:.6f}" for a in AA]
    return row


def _frac(hi):
    f = {a: (1 - 0.6) / 19 for a in AA}
    f[hi] = 0.6
    return f


def _setup(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(common, "CONFIG_DIR", tmp_path / "config")
    monkeypatch.setattr(common, "DATA_DIR", tmp_path / "data")
    monkeypatch.setattr(common, "OUTPUTS_DIR", tmp_path / "outputs")
    monkeypatch.setattr(common, "LOGS_DIR", tmp_path / "logs")
    monkeypatch.setattr(common, "LITERATURE_SOURCES_INI", tmp_path / "config" / "literature_sources.ini")
    cfg = tmp_path / "config" / "cat"; cfg.mkdir(parents=True)
    dat = tmp_path / "data" / "cat"; dat.mkdir(parents=True)
    lit = tmp_path / "data" / "literature" / "s"; lit.mkdir(parents=True)
    comp_dir = tmp_path / "outputs" / "cat" / "intermediate" / "composition"; comp_dir.mkdir(parents=True)
    ps_dir = tmp_path / "outputs" / "cat" / "intermediate" / "protein_set"; ps_dir.mkdir(parents=True)

    # -- the two primary tables: red (4 replicate shares, percent), blue (log10 abundance, 3 runs, cv)
    red_rows = [
        ["X00001;X00001-2", "GA", 40, 44, 36, 40],        # member
        ["X00002", "GB", 20, 20, 20, 20],                 # member
        ["X00002-2", "GB", 5, 5, 5, 5],                   # a second row of X00002 -> summed
        ["S00003", "GC;GD", 10, 10, 10, 10],              # shared between X00003 and X00004 (equal: no own rows)
        ["X00005", "KRT1", 15, 11, 19, 15],               # keratin, removed by rule
        ["X00006", "GF", 10, 10, 10, 10],                 # R5-excluded entry (not in accessions.ini)
    ]
    _write_xlsx(lit / "red.xlsx", "R", ["id", "gene", "d1", "d2", "d3", "d4"], red_rows)
    blue_rows = [
        ["P00001", "IGHG1", 3.0, 3.0, 3.0, 3.0, 4.0],      # constant heavy chain; v = 1000
        ["P00002", "IGHV3-23", 2.0, 2.0, 2.0, 2.0, 2.0],   # variable heavy segment; v = 100
        ["P00003", "IGKC", 2.0, 2.0, "", 2.0, 5.0],        # constant kappa; no IGKV in the pool -> loose; one run blank
        ["P00004", "ALB", math.log10(8800), 3.9, 3.95, 3.97, 1.0],
    ]
    _write_xlsx(lit / "blue.xlsx", "B", ["id", "gene", "avg", "r1", "r2", "r3", "cv"], blue_rows)
    _write_xlsx(lit / "deep.xlsx", "D", ["id", "gene", "avg"], [["P00001;P00009", "IGHG1", 3.0], ["P00004", "ALB", 4.0], ["Q00001", "NEW", 3.0]])
    (lit / "hortin.xls").write_bytes(b"stub")
    man = common.new_ini()
    for key, name in (("s.file.1", "red.xlsx"), ("s.file.2", "blue.xlsx"), ("s.file.3", "deep.xlsx"), ("s.file.4", "hortin.xls")):
        man[f"file.{key}"] = {"path": f"data/literature/s/{name}", "sha256": common.sha256_file(lit / name)}
    common.write_ini(man, tmp_path / "data" / "literature" / "manifest.ini", ["test"])
    src = common.new_ini(); src["source.s"] = {"citation": "Synthetic S. 2026."}
    common.write_ini(src, common.LITERATURE_SOURCES_INI, ["test"])

    psd = common.new_ini()
    psd["scope"] = {"category": "cat"}
    psd["pool.red"] = {"dataset_source": "s", "dataset_file": "file.1"}
    psd["pool.blue"] = {"dataset_source": "s", "dataset_file": "file.2"}
    psd["contaminant_rules"] = {"keratin_keyword_id": "KW-0416"}
    common.write_ini(psd, cfg / "protein_set_decisions.ini", ["test"])
    dec = common.new_ini()
    dec["scope"] = {"category": "cat"}
    dec["columns.s.file.1"] = {"sheet": "R", "header_row": "1", "identity_column": "A", "identity_kind": "accession_list", "gene_column": "B",
                               "share_columns": "C, D, E, F", "share_kind": "mass_share", "share_scale": "percent"}
    dec["columns.s.file.2"] = {"sheet": "B", "header_row": "1", "identity_column": "A", "identity_kind": "accession_list", "gene_column": "B",
                               "abundance_column": "C", "abundance_scale": "log10", "share_kind": "mass_share",
                               "replicate_columns": "D, E, F", "cv_column": "G"}
    dec["columns.s.file.3"] = {"sheet": "D", "header_row": "1", "identity_column": "A", "gene_column": "B", "abundance_column": "C", "abundance_scale": "log10"}
    dec["columns.s.file.4"] = {"sheet": "H", "header_row": "1", "identity_column": "C", "gene_column": "D", "value_column": "G", "name_column": "B", "value_unit": "mg/L"}
    dec["cross_checks"] = {"hortin_source": "s", "hortin_file": "file.4", "hortin_pool": "blue", "deep_source": "s", "deep_file": "file.3", "deep_pool": "blue"}
    dec["immunoglobulin_bound"] = {"chain_prefixes": "IGH, IGK, IGL", "variable_segment_letters": "V, D, J"}
    common.write_ini(dec, cfg / "mass_fraction_decisions.ini", ["test"])

    # -- the enumeration's outputs (join tables, excluded rows, pool files)
    rows_cols = ["row_index", "identity_cell", "gene_cell", "entry", "entry_token_rank", "outcome", "tokens_checked", "published_share_percent"]
    common.write_tsv(dat / "red_dataset_rows.tsv", rows_cols, [
        [2, "X00001;X00001-2", "GA", "X00001", "1", "member", "", "40"],
        [3, "X00002", "GB", "X00002", "1", "member", "", "20"],
        [4, "X00002-2", "GB", "X00002", "1", "member_duplicate_entry", "", "5"],
        [5, "S00003", "GC;GD", "X00003;X00004", "1", "member_via_secondary_accession_shared", "", "10"],
        [6, "X00005", "KRT1", "X00005", "1", "keratin", "", "15"],
        [7, "X00006", "GF", "X00006", "1", "member", "", "10"],
    ])
    common.write_tsv(dat / "red_rows_excluded_with_share.tsv", ["row_index", "identity_cell", "gene_cell", "entry", "outcome", "published_share_percent"],
                     [[6, "X00005", "KRT1", "X00005", "keratin", "15"]])
    common.write_tsv(dat / "blue_dataset_rows.tsv", rows_cols, [
        [2, "P00001", "IGHG1", "P00001", "1", "member", "", ""],
        [3, "P00002", "IGHV3-23", "P00002", "1", "member", "", ""],
        [4, "P00003", "IGKC", "P00003", "1", "member", "", ""],
        [5, "P00004", "ALB", "P00004", "1", "member", "", ""],
    ])
    common.write_tsv(dat / "blue_rows_excluded_with_share.tsv", ["row_index", "identity_cell", "gene_cell", "entry", "outcome", "published_share_percent"], [])
    pool_cols = ["accession", "gene_primary", "protein_name", "length", "reviewed", "organism_id", "dataset_rows", "identity_cell"]
    common.write_tsv(dat / "pool_red.tsv", pool_cols, [[a, "", "", "", "reviewed", "9606", "", ""] for a in ("X00001", "X00002", "X00003", "X00004", "X00006")])
    common.write_tsv(dat / "pool_blue.tsv", pool_cols, [[a, "", "", "", "reviewed", "9606", "", ""] for a in ("P00001", "P00002", "P00003", "P00004")])

    # -- the protein set (X00006 excluded under R5, so absent) and the composition
    acc = common.new_ini()
    for a, g, t in (("X00001", "GA", "red"), ("X00002", "GB", "red"), ("X00003", "GC", "red"), ("X00004", "GD", "red"),
                    ("P00001", "IGHG1", "blue"), ("P00002", "IGHV3-23", "blue"), ("P00003", "IGKC", "blue"), ("P00004", "ALB", "red;blue")):
        acc[a] = {"gene": g, "tier": t, "flag_open": "false"}
    common.write_ini(acc, cfg / "accessions.ini", ["GENERATED by build_protein_set.py"])
    # P00004 (ALB) also in red with no red row -> w_red = 0
    comp_header = ["accession", "gene", "tier", "segment_set", "n_residues", "residue_mass_sum", "mw", "free_mass_sum"] + [f"count_{a}" for a in AA] \
        + [f"residue_frac_{a}" for a in AA] + [f"free_frac_{a}" for a in AA]
    comp_rows = [_composition_row("X00001", "GA", "red", 100, 11000.0, _frac("A")), _composition_row("X00002", "GB", "red", 100, 11000.0, _frac("L")),
                 _composition_row("X00003", "GC", "red", 100, 11000.0, _frac("K")), _composition_row("X00004", "GD", "red", 100, 11000.0, _frac("K")),
                 _composition_row("P00001", "IGHG1", "blue", 330, 36000.0, _frac("S")), _composition_row("P00002", "IGHV3-23", "blue", 110, 12000.0, _frac("S")),
                 _composition_row("P00003", "IGKC", "blue", 107, 11700.0, _frac("T")), _composition_row("P00004", "ALB", "red;blue", 585, 66000.0, _frac("E"))]
    common.write_tsv(comp_dir / "amino_acid_composition_per_protein.tsv", comp_header, comp_rows)
    common.write_tsv(comp_dir / "processing_mass_deltas.tsv", ["accession"] + [f"delta_free_frac_{a}" for a in AA],
                     [[a] + ["0.01" if a2 == "A" else "0" for a2 in AA] for a in acc.sections()])
    common.write_tsv(comp_dir / "ptm_mass_deltas.tsv", ["accession", "ptm_mass_delta_fraction_of_mw", "n_glycosylation"],
                     [[a, "0.002" if a == "P00001" else "0", "1" if a == "P00001" else "0"] for a in acc.sections()])
    common.write_tsv(ps_dir / "isoform_deltas.tsv", ["accession"] + [f"delta_{a}" for a in AA], [["X00001"] + ["0.05" if a == "A" else "0" for a in AA]])

    hortin_rows = [
        {"row_index": 2, "identity": "P00004", "gene": "ALB", "value": 43500, "name": "Albumin"},
        {"row_index": 3, "identity": "", "gene": "", "value": 4300, "name": "Ig gamma1"},
        {"row_index": 4, "identity": "P00003-2", "gene": "IGKC", "value": 2700, "name": "Ig kappa"},
        {"row_index": 5, "identity": "NP 001629", "gene": "X", "value": 100, "name": "RefSeq id"},
        {"row_index": 6, "identity": "Q99999", "gene": "NOTINPOOL", "value": 50, "name": "not in pool"},
    ]
    monkeypatch.setattr(mfp, "read_xls_by_letter", lambda path, sheet, header_row, letters, log: hortin_rows)
    return tmp_path


def test_pools_weights_bounds_and_checks(tmp_path, monkeypatch):
    _setup(tmp_path, monkeypatch)
    log = logging.getLogger("t")
    assert mfp.main(["--category", "cat"]) == 0
    per = {r["accession"]: r for r in read(tmp_path / "config" / "cat" / "mass_fractions_per_entry.tsv")}
    # red: X00001 v = mean(40,44,36,40)/100 = 0.40; X00002 summed 0.20+0.05; shared 0.10 split equally; listing 1.00, denominator 0.75 (keratin 0.15 and R5 0.10 out)
    assert math.isclose(float(per["X00001"]["v_red"]), 0.40) and math.isclose(float(per["X00001"]["w_red"]), 0.40 / 0.75)
    assert per["X00002"]["match_rule_red"] == "summed:2" and math.isclose(float(per["X00002"]["v_red"]), 0.25)
    assert per["X00003"]["match_rule_red"] == "shared_split" and math.isclose(float(per["X00003"]["v_red"]), 0.05)
    assert math.isclose(float(per["X00003"]["w_red_high"]), (0.05 + 0.0 + 0.05) / 0.75)   # sd 0 across equal replicates; the whole shared row
    assert per["X00001"]["n_red"] == "4" and math.isclose(float(per["X00001"]["sd_red"]), 0.032659863, rel_tol=1e-6)
    assert [per["X00001"][f"v_red_rep{k}"] for k in (1, 2, 3, 4)] == ["0.4", "0.44", "0.36", "0.4"]
    assert per["P00004"]["match_rule_red"] == "none" and per["P00004"]["w_red"] == "0"
    reds = [float(per[a]["w_red"]) for a in ("X00001", "X00002", "X00003", "X00004", "P00004")]
    assert math.isclose(sum(reds), 1.0)
    # blue: 10^avg; ALB 8800 of 1000+100+100+8800
    assert math.isclose(float(per["P00004"]["w_blue"]), 8800 / 10000) and per["P00003"]["n_blue"] == "2"   # one blank run
    assert per["P00001"]["cv_percent_published_blue"] == "4"
    assert math.isclose(sum(float(per[a]["w_blue"]) for a in ("P00001", "P00002", "P00003", "P00004")), 1.0)
    out = tmp_path / "outputs" / "cat" / "intermediate" / "mass_fractions"
    summ = common.read_ini(out / "mass_fractions_summary.ini")
    assert summ["pool.red"]["removed_by_rule_share_of_listing"] == "0.15" and summ["pool.red"]["completeness_gap_share_of_listing"] == "0.1"
    assert summ["pool.red"]["entries_not_in_accessions_ini"] == "1"
    exc = read(out / "excluded_rows_mass_share.tsv")
    assert {(r["entry"], r["kind"]) for r in exc} == {("X00005", "removed_by_rule"), ("X00006", "completeness_gap")}
    # D113: IGHG1 bound = w x 110/330; IGKC loose (no IGKV in the pool); IGHV3-23 measured, no bound
    ig = {r["accession"]: r for r in read(out / "immunoglobulin_variable_domain_bound.tsv")}
    assert ig["P00001"]["segment_class"] == "constant" and math.isclose(float(ig["P00001"]["variable_mass_bound"]), 0.1 * 110 / 330)
    assert ig["P00002"]["segment_class"] == "variable_segment" and ig["P00002"]["variable_mass_bound"] == "0"
    assert ig["P00003"]["how"].startswith("loose") and math.isclose(float(ig["P00003"]["variable_mass_bound"]), 0.01)
    chains = {(r["pool"], r["chain"]): r for r in read(out / "immunoglobulin_variable_domain_bound_per_chain.tsv")}
    assert math.isclose(float(chains[("blue", "all")]["variable_mass_bound"]), 0.1 / 3 + 0.01)
    # bounds: X00001's isoform delta 0.05 x w
    wb = {(r["pool"], r["kind"], r["amino_acid"]): r for r in read(out / "weighted_bounds.tsv")}
    assert math.isclose(float(wb[("red", "isoform", "A")]["max_weighted_delta"]), 0.05 * 0.40 / 0.75) and wb[("red", "isoform", "A")]["where"] == "X00001"
    assert math.isclose(float(wb[("blue", "glycosylated_entries_sum_w", "all")]["max_weighted_delta"]), 0.1)
    # Hortin: ALB and IGKC join (isoform suffix stripped); the blank, RefSeq and not-in-pool rows are listed with why
    ratio = {r["accession"]: r for r in read(out / "ratio_check_hortin_2008.tsv")}
    assert set(ratio) == {"P00004", "P00003"} and math.isclose(float(ratio["P00004"]["share_hortin"]), 43500 / 46200)
    assert math.isclose(float(ratio["P00004"]["share_geyer_joined"]), 8800 / 8900)
    nj = read(out / "hortin_2008_rows_not_joined.tsv")
    assert [r["why"] for r in nj] == ["no UniProt accession in the cell", "not a UniProt accession", "Q99999 not in the blue pool"]
    prof = {r["amino_acid"]: r for r in read(out / "profile_plasma_geyer_vs_hortin.tsv")}
    assert float(prof["E"]["profile_hortin_joined_set"]) > float(prof["T"]["profile_hortin_joined_set"])
    assert summ["cross_check.hortin_2008"]["hortin_rows_joined"] == "2" and summ["cross_check.hortin_2008"]["hortin_rows_not_joined"] == "3"
    # deep dataset: Q00001 absent from the primary; 1000 of 1000+10000+1000
    deep = read(out / "completeness_check_geyer_2016_deep_dataset.tsv")
    assert [r["identity_cell"] for r in deep] == ["Q00001"] and math.isclose(float(deep[0]["share_of_deep_signal"]), 1000 / 12000)
    # overlap: ALB in both pools
    ov = read(out / "pool_overlap_mass_share.tsv")
    assert [r["accession"] for r in ov] == ["P00004"] and math.isclose(float(summ["overlap"]["sum_w_blue"]), 0.88)
    # currency
    assert mfp.main(["--category", "cat", "--check"]) == 0
    p = tmp_path / "config" / "cat" / "mass_fractions_per_entry.tsv"
    p.write_text(p.read_text().replace("summed:2", "summed:3"))
    assert mfp.main(["--category", "cat", "--check"]) == 1


def test_molar_kind_is_refused(tmp_path, monkeypatch):
    """D105: only mass_share is built; any other kind stops the stage until a source needs it."""
    _setup(tmp_path, monkeypatch)
    p = tmp_path / "config" / "cat" / "mass_fraction_decisions.ini"
    import re
    p.write_text(re.sub(r"share_kind(\s*)=(\s*)mass_share", r"share_kind\1=\2molar", p.read_text(), count=1))
    try:
        mfp.main(["--category", "cat"])
    except SystemExit as e:
        assert "molar" in str(e) or "mass_share" in str(e)
    else:
        raise AssertionError("expected a STOP")
