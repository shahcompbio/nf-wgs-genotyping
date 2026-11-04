#!/usr/bin/env python3
"""Run svtyper (marcjwilliams1/svtyper fork) on one chunk VCF, or only build the library-stats JSON.

Without --reference, the clip-matching flags (--clip_read_support, --clip_partner_match) are
removed with a warning, because matching clipped reads needs -T. Pass svtyper flags as one string
with the '=' form, e.g. --extra-args='--both_sides --keep_all_ref'.
"""
import argparse
import os
import shlex
import subprocess
import sys

CLIP_NEEDS_REF = ("--clip_read_support", "--clip_partner_match")


def build_command(svtyper_cmd, vcf, bam, lib_info, out, reference=None, extra_args="", lib_only=False):
    base = shlex.split(svtyper_cmd)
    if lib_only:
        return base + ["-B", bam, "-l", lib_info, "-i", "/dev/null", "-o", "/dev/null"], []
    args = shlex.split(extra_args or "")
    warnings = []
    if not reference:
        removed = [a for a in args if a in CLIP_NEEDS_REF]
        if removed:
            args = [a for a in args if a not in CLIP_NEEDS_REF]
            warnings.append(f"no reference FASTA: dropping clip-matching flags {removed} (they need -T)")
    cmd = base + ["-i", vcf, "-B", bam, "-l", lib_info, "-o", out]
    if reference:
        cmd += ["-T", reference]
    return cmd + args, warnings


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--svtyper-cmd", default="svtyper")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--lib-info", required=True)
    ap.add_argument("--vcf")
    ap.add_argument("--out")
    ap.add_argument("--reference", default=None)
    ap.add_argument("--extra-args", default="")
    ap.add_argument("--lib-only", action="store_true")
    a = ap.parse_args(argv)
    if not a.lib_only and not (a.vcf and a.out):
        ap.error("--vcf and --out are required unless --lib-only")
    cmd, warnings = build_command(a.svtyper_cmd, a.vcf, a.bam, a.lib_info, a.out, a.reference or None,
                                  a.extra_args, a.lib_only)
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    print("+ " + " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    if a.lib_only and not os.path.exists(a.lib_info):
        raise SystemExit(f"svtyper did not write {a.lib_info}")
    if not a.lib_only and not os.path.exists(a.out):
        open(a.out, "w").close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
