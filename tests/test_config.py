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


def test_container_is_a_reference_singularity_and_docker_can_pull():
    # Nextflow prefixes `docker://` to any name without a `scheme://`, so 'docker:<ns>/<img>:<tag>'
    # becomes the unpullable 'docker://docker:<ns>/<img>:<tag>' and a bare '<img>:<tag>' resolves
    # to docker.io/library/<img>. Require [docker://][registry[:port]/]<ns>/<img>:<tag>.
    m = re.search(r"^\s*container\s*=\s*'([^']*)'", _config_text(), re.M)
    assert m, "params.container not found in nextflow.config"
    ref = m.group(1)
    pattern = r"(docker://)?([a-z0-9.-]+(:[0-9]+)?/)?[a-z0-9._-]+/[a-z0-9._/-]+:[A-Za-z0-9._-]+"
    assert re.fullmatch(pattern, ref), f"params.container {ref!r} is not a pullable image reference"
