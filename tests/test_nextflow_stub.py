"""Nextflow -stub-run smoke tests for both entry scripts (skipped when nextflow is not on PATH).

Runs locally with the `stub` profile; NXF_HOME defaults to a temp dir and NXF_OFFLINE=true so
no network or $HOME writes are needed."""
import os
import shutil
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pytestmark = pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow not on PATH")


def _run(entry, tmp_path, *extra):
    env = dict(os.environ)
    env.setdefault("NXF_HOME", str(tmp_path / "nxf_home"))
    env["NXF_OFFLINE"] = "true"
    env.setdefault("NXF_ANSI_LOG", "false")
    out = tmp_path / "out"
    # local checkouts are launched by script path; `-main-script` only resolves inside remote repos
    cmd = ["nextflow", "run", os.path.join(REPO, entry), "-profile", "stub", "-stub-run",
           "-work-dir", str(tmp_path / "work"), "--outdir", str(out), *extra]
    res = subprocess.run(cmd, cwd=str(tmp_path), env=env, capture_output=True, text=True, timeout=600)
    return res, out


def test_snv_entry_stub_run_publishes_final_outputs(tmp_path):
    res, out = _run("snv_genotyping.nf", tmp_path)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    for name in ("STUB_snv_genotyping.csv.gz", "STUB_snv_genotyping.csv.gz.yaml", "STUB_genotyping_qc.json",
                 "union/STUB_union_snvs.tsv.gz", "union/STUB_union_snvs.vcf.gz",
                 "pipeline_info/software_versions.yml"):
        assert (out / name).exists(), name
    assert "SNV_COUNT" in res.stdout


def test_sv_entry_stub_run_publishes_final_outputs(tmp_path):
    res, out = _run("sv_genotyping.nf", tmp_path)
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    for name in ("STUB_sv_genotyping.csv.gz", "STUB_sv_genotyping.csv.gz.yaml", "STUB_genotyping_qc.json",
                 "STUB_svtyper.vcf.gz", "STUB.lib_info.json", "union/STUB_union_svs.bedpe.gz",
                 "union/STUB_union_svs.members.tsv.gz", "union/STUB_union_svs.vcf.gz"):
        assert (out / name).exists(), name
    assert "SVTYPER" in res.stdout


def test_main_nf_refuses_to_run(tmp_path):
    env = dict(os.environ, NXF_OFFLINE="true", NXF_ANSI_LOG="false")
    env.setdefault("NXF_HOME", str(tmp_path / "nxf_home"))
    res = subprocess.run(["nextflow", "run", REPO, "-profile", "stub", "-stub-run", "-work-dir",
                          str(tmp_path / "work")], cwd=str(tmp_path), env=env, capture_output=True, text=True,
                         timeout=600)
    assert res.returncode != 0
    assert "-main-script" in (res.stdout + res.stderr)


def test_snv_entry_rejects_sv_callers(tmp_path):
    ss = tmp_path / "ss.csv"
    vcf = tmp_path / "x.vcf"
    vcf.write_text("")
    ss.write_text(f"tumor_id,normal_id,caller,vcf\nT,N,brass,{vcf}\n")
    res, _ = _run("snv_genotyping.nf", tmp_path, "--samplesheet", str(ss))
    assert res.returncode != 0
    assert "brass" in (res.stdout + res.stderr)
