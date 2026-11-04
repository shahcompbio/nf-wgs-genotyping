"""Locate the svtyper fork for integration tests.

Set SVTYPER_SRC to a checkout of https://github.com/marcjwilliams1/svtyper (commit 27a0427), or
install it; otherwise the integration tests are skipped. The default svtyper flags tested here must
match params.svtyper_args in nextflow.config (checked by tests/test_config.py)."""
import os
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DEFAULT_SVTYPER_ARGS = ("--genotype_on_fragments --keep_all_ref --both_sides --clip_read_support "
                        "--clip_partner_match --clip_min_base_quality 20")


def svtyper_src():
    cand = os.environ.get("SVTYPER_SRC") or os.path.join(REPO, "..", "software", "svtyper")
    return os.path.abspath(cand) if os.path.isdir(os.path.join(cand, "svtyper")) else None


def svtyper_env():
    env = dict(os.environ)
    src = svtyper_src()
    if src:
        env["PYTHONPATH"] = src + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    return env


def svtyper_cmd():
    return f"{shlex.quote(sys.executable)} -m svtyper.classic"


def available():
    try:
        subprocess.run([sys.executable, "-c", "import svtyper.classic"], env=svtyper_env(), check=True,
                       capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False


def import_parsers():
    src = svtyper_src()
    if src and src not in sys.path:
        sys.path.insert(0, src)
    from svtyper import parsers  # noqa: E402

    return parsers
