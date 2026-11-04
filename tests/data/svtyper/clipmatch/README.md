# clipmatch fixture

Three breakend pairs (`sv1`-`sv3`) for the clipped-read matching options
(`--clip_partner_match`, `--clip_min_base_quality`), used by `tests/test_clip_partner.py`.
Reads are from targeted, UMI-collapsed / uncollapsed sequencing, de-identified:
each breakend region is its own contig (`ctg1`-`ctg6`) with coordinates starting at 1,
reads are renamed, all tags except `SA` and one read group are dropped, and mates
outside the regions are marked unmapped. Only the reads each case needs are kept.

| SV | what it is | expected |
|---|---|---|
| sv1 | a real junction: clipped, split and spanning reads (consensus reads) | clips kept by both options |
| sv2 | reads clipped for other reasons; the clip is the read's own reference continuation, on the wrong side for the breakend (consensus reads) | clips counted by the default matcher, rejected by `--clip_partner_match` |
| sv3 | clips with low base quality in the clipped bases (uncollapsed reads) | removed by `--clip_min_base_quality 30` |

With `--fragment_counts`, sv1's 25 alternate reads are 19 fragments (FS 8, FP 7, FC 4,
none with reference evidence), while every sv2 / sv3 "clip" fragment also carries reference
evidence (FX 9), so they count as reference (AOF 0).

Run with `-T ref.fa -l lib.json --clip_read_support --both_sides --keep_all_ref -m 25`.
