import shlex
import sys

import svtyper_run


def test_full_command_with_reference_keeps_clip_flags():
    cmd, warnings = svtyper_run.build_command("svtyper", "c.vcf", "t.bam", "lib.json", "o.vcf",
                                              reference="ref.fa", extra_args="--both_sides --clip_partner_match")
    assert cmd == ["svtyper", "-i", "c.vcf", "-B", "t.bam", "-l", "lib.json", "-o", "o.vcf", "-T", "ref.fa",
                   "--both_sides", "--clip_partner_match"]
    assert warnings == []


def test_without_reference_clip_matching_flags_are_dropped_with_warning():
    cmd, warnings = svtyper_run.build_command(
        "python -m svtyper.classic", "c.vcf", "t.bam", "lib.json", "o.vcf", reference=None,
        extra_args="--keep_all_ref --clip_read_support --clip_partner_match --clip_min_base_quality 20")
    assert cmd[:3] == ["python", "-m", "svtyper.classic"]
    assert "-T" not in cmd
    assert "--clip_read_support" not in cmd and "--clip_partner_match" not in cmd
    assert cmd[-3:] == ["--keep_all_ref", "--clip_min_base_quality", "20"]
    assert any("clip" in w for w in warnings)


def test_lib_only_command_uses_empty_input():
    cmd, _ = svtyper_run.build_command("svtyper", None, "t.bam", "lib.json", None, lib_only=True)
    assert cmd == ["svtyper", "-B", "t.bam", "-l", "lib.json", "-i", "/dev/null", "-o", "/dev/null"]


def test_main_runs_command_and_creates_output(tmp_path):
    fake = tmp_path / "fake_svtyper.py"
    fake.write_text("import sys\nargs=sys.argv[1:]\nopen(args[args.index('-o')+1],'w').write(' '.join(args))\n")
    out = tmp_path / "o.vcf"
    rc = svtyper_run.main(["--svtyper-cmd", f"{shlex.quote(sys.executable)} {fake}", "--vcf", "c.vcf",
                           "--bam", "t.bam", "--lib-info", "lib.json", "--out", str(out),
                           "--extra-args=--both_sides"])
    assert rc == 0
    assert out.read_text().endswith("--both_sides")
