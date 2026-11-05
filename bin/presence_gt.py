"""Allele-presence genotypes from fragment counts (shared by snv_finalize and indel_finalize).

A tumor's alt fraction depends on purity, copy number and clonality, so the diploid dosage model
(alt fraction 0, 1/2, 1) mislabels most somatic variants. GT here says which alleles are
*present* in the sample, not how many copies there are. With k alt of n = ref + alt fragments:

  0/0  only ref is present:  k   ~ BetaBinomial(n, mean eps, overdispersion rho)
  1/1  only alt is present:  n-k ~ BetaBinomial(n, mean eps, overdispersion rho)
  0/1  both are present:     k ~ Binomial(n, f), f ~ Uniform[min_vaf, 1 - min_vaf]
       L = [I_{1-min_vaf}(k+1, n-k+1) - I_{min_vaf}(k+1, n-k+1)] / ((n + 1)(1 - 2 min_vaf))

eps is the rate at which a fragment shows this specific wrong allele (sequencing/alignment error);
rho (0 = binomial) lets that rate vary between sites. PL = -10 log10(L_g / max L) (VCF order
0/0, 0/1, 1/1); with a flat prior, GQ = -10 log10(1 - P(best genotype)), capped at 99. GT is the
best genotype when GQ >= min_gq, else './.' (no-call): low depth, zero depth (PL 0,0,0) and
ambiguous counts are not forced into a genotype.
"""
import numpy as np
from scipy import stats
from scipy.special import logsumexp

GENOTYPES = ("0/0", "0/1", "1/1")
GQ_MAX = 99
PL_MAX = 10 ** 6          # cap for genotypes whose likelihood underflows
_PHRED = -10.0 / np.log(10.0)


def default_error_rate(min_baseq):
    """Worst-case per-base error allowed by the base-quality filter (10^(-BQ/10)), split evenly
    over the three wrong bases: the rate of one specific wrong base. Conservative."""
    return 10.0 ** (-float(min_baseq) / 10.0) / 3.0


def check_params(error_rate, overdispersion, min_vaf, min_gq):
    if not 0 <= min_vaf < 0.5:
        raise ValueError(f"min_vaf must be in [0, 0.5): {min_vaf}")
    if not 0 < error_rate < 1:
        raise ValueError(f"error_rate must be in (0, 1): {error_rate}")
    if min_vaf > 0 and error_rate >= min_vaf:
        raise ValueError(f"error_rate ({error_rate}) must be below min_vaf ({min_vaf})")
    if not 0 <= overdispersion < 1:
        raise ValueError(f"overdispersion must be in [0, 1): {overdispersion}")
    if min_gq < 0:
        raise ValueError(f"min_gq must be >= 0: {min_gq}")


def _log_error_only(k, n, eps, rho):
    """log P(k wrong-allele fragments of n) when the allele is absent."""
    if rho == 0:
        return stats.binom.logpmf(k, n, eps)
    a, b = eps * (1 - rho) / rho, (1 - eps) * (1 - rho) / rho
    return stats.betabinom.logpmf(k, n, a, b)


def _log_present(k, n, fmin):
    """log of the 0/1 marginal likelihood (binomial averaged over f ~ U[fmin, 1 - fmin]).
    The Beta(k+1, n-k+1) mass on [fmin, 1 - fmin] is a difference of two tails; the tails used
    are the ones that are small for these counts, so it does not cancel to 0."""
    a, b = k + 1.0, n - k + 1.0
    log_sf_lo, log_sf_hi = stats.beta.logsf(fmin, a, b), stats.beta.logsf(1 - fmin, a, b)
    log_cdf_hi, log_cdf_lo = stats.beta.logcdf(1 - fmin, a, b), stats.beta.logcdf(fmin, a, b)
    with np.errstate(divide="ignore", invalid="ignore"):
        via_sf = log_sf_lo + np.log1p(-np.exp(np.minimum(log_sf_hi - log_sf_lo, 0.0)))
        via_cdf = log_cdf_hi + np.log1p(-np.exp(np.minimum(log_cdf_lo - log_cdf_hi, 0.0)))
    log_mass = np.where(log_sf_lo <= log_cdf_hi, via_sf, via_cdf)
    return log_mass - np.log(n + 1.0) - np.log(1.0 - 2.0 * fmin)


def log_likelihoods(ref, alt, error_rate, overdispersion, min_vaf):
    """Natural-log likelihoods, shape (sites, 3), in VCF genotype order 0/0, 0/1, 1/1."""
    k = np.asarray(alt, dtype=float)
    n = k + np.asarray(ref, dtype=float)
    return np.column_stack([_log_error_only(k, n, error_rate, overdispersion),
                            _log_present(k, n, min_vaf),
                            _log_error_only(n - k, n, error_rate, overdispersion)])


def call(ref, alt, error_rate, overdispersion=0.01, min_vaf=0.05, min_gq=20):
    """Vectorized calls. Returns (gt: list of '0/0'|'0/1'|'1/1'|'./.', gq: int array,
    pl: int array of shape (sites, 3))."""
    check_params(error_rate, overdispersion, min_vaf, min_gq)
    ll = log_likelihoods(ref, alt, error_rate, overdispersion, min_vaf)
    if len(ll) == 0:
        return [], np.zeros(0, dtype=int), np.zeros((0, 3), dtype=int)
    best = np.argmax(ll, axis=1)
    rows = np.arange(len(ll))
    with np.errstate(divide="ignore", invalid="ignore"):
        pl = np.rint(_PHRED * (ll - ll[rows, best][:, None]))
    pl = np.nan_to_num(pl, nan=PL_MAX, posinf=PL_MAX).clip(0, PL_MAX).astype(int)
    others = ll.copy()
    others[rows, best] = -np.inf
    with np.errstate(divide="ignore"):
        gq = _PHRED * (logsumexp(others, axis=1) - logsumexp(ll, axis=1))
    gq = np.rint(np.nan_to_num(gq, posinf=GQ_MAX)).clip(0, GQ_MAX).astype(int)
    gt = [GENOTYPES[b] if q >= min_gq else "./." for b, q in zip(best, gq)]
    return gt, gq, pl


def format_pl(pl):
    return ",".join(str(int(x)) for x in pl)


def gt_tuple(gt):
    """'0/1' -> (0, 1); './.' -> (None, None) (pysam sample['GT'] values)."""
    return tuple(None if a == "." else int(a) for a in gt.split("/"))


def vcf_header_lines(error_rate, overdispersion, min_vaf, min_gq):
    return (
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Allele presence, not dosage: 0/0 = only ref '
        'detected (alt not above error), 0/1 = both alleles present at any fraction (heterozygous, '
        'subclonal or impure), 1/1 = only alt detected; ./. = GQ below the threshold">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred probability that GT is wrong (flat '
        'prior over 0/0, 0/1, 1/1; capped at 99)">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods of 0/0, 0/1, 1/1 '
        'under the allele-presence model (normalized to 0)">',
        f'##nf-wgs-genotyping_gt_model=<model=allele_presence,error_rate={error_rate:.6g},'
        f'overdispersion={overdispersion:g},min_vaf={min_vaf:g},min_gq={min_gq:g}>',
    )
