# svtyper test fixtures (third-party, MIT)

Copied verbatim from the svtyper fork used by this pipeline
(https://github.com/marcjwilliams1/svtyper, commit 27a042710ddc6c2efa224b6b6cc50b8b45231c0c,
`tests/data/`). License: MIT, see `LICENSE.svtyper`. Authors: Colby Chiang, Marc J. Williams.

- `NA12878.target_loci.sorted.bam(.bai)`, `NA12878.bam.json`, `example.vcf`, `example.gt.vcf`:
  public NA12878 reads at LUMPY-called loci and svtyper's expected output.
- `clipmatch/`: de-identified clipped-read fixture (see `clipmatch/README.md`).

These VCFs follow LUMPY/svtyper conventions (a reverse-orientation BND's POS is the base
before the junction); caller VCFs handled by `sv_normalize.py` follow the VCF 4.2 spec.
