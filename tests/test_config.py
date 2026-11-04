import os
import re

import svtyper_helpers as sh

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _config_text():
    with open(os.path.join(REPO, "nextflow.config")) as fh:
        return fh.read()


def test_default_svtyper_args_match_tested_flags():
    m = re.search(r"svtyper_args\s*=\s*'([^']*)'", _config_text())
    assert m, "params.svtyper_args not found in nextflow.config"
    assert m.group(1) == sh.DEFAULT_SVTYPER_ARGS


def test_both_entry_scripts_exist_and_main_is_a_usage_stub():
    for entry in ("snv_genotyping.nf", "sv_genotyping.nf"):
        assert os.path.isfile(os.path.join(REPO, entry)), entry
    main = open(os.path.join(REPO, "main.nf")).read()
    assert "-main-script snv_genotyping.nf" in main and "-main-script sv_genotyping.nf" in main


def test_manifest_version_matches_changelog():
    version = re.search(r"version\s*=\s*'([^']+)'", _config_text()).group(1)
    assert f"## v{version}" in open(os.path.join(REPO, "CHANGELOG.md")).read()
