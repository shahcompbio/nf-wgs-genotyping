"""Allele-presence genotype model (bin/presence_gt.py): likelihoods checked against brute-force
numerical integration, and the calls/no-calls the docs promise."""
import math

import numpy as np
import pytest

import presence_gt as pg

EPS = pg.default_error_rate(20)          # 10^-2 / 3
RHO = 0.01
FMIN = 0.05


def _lbeta(a, b):
    return math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)


def _log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _brute_betabinom(k, n, mean, rho):
    """Beta-binomial pmf from its closed form with math.lgamma (independent of scipy)."""
    a, b = mean * (1 - rho) / rho, (1 - mean) * (1 - rho) / rho
    return math.exp(_log_choose(n, k) + _lbeta(k + a, n - k + b) - _lbeta(a, b))


def _brute_present(k, n, fmin, steps=200001):
    f = np.linspace(fmin, 1 - fmin, steps)
    dens = np.exp(_log_choose(n, k) + k * np.log(f) + (n - k) * np.log1p(-f))
    return np.trapz(dens, f) / (1 - 2 * fmin)


def _call(ref, alt, **kw):
    kw = {"error_rate": EPS, "overdispersion": RHO, "min_vaf": FMIN, "min_gq": 20, **kw}
    gt, gq, pl = pg.call(np.array([ref]), np.array([alt]), **kw)
    return gt[0], int(gq[0]), tuple(int(x) for x in pl[0])


def test_default_error_rate_is_the_bq_bound_split_over_three_bases():
    assert pg.default_error_rate(20) == pytest.approx(0.01 / 3)
    assert pg.default_error_rate(30) == pytest.approx(0.001 / 3)


@pytest.mark.parametrize("ref,alt", [(40, 0), (36, 4), (20, 20), (3, 0), (0, 40), (7, 2)])
def test_likelihoods_match_numerical_integration(ref, alt):
    n = ref + alt
    ll = pg.log_likelihoods(np.array([ref]), np.array([alt]), EPS, RHO, FMIN)[0]
    assert math.exp(ll[0]) == pytest.approx(_brute_betabinom(alt, n, EPS, RHO), rel=1e-6)
    assert math.exp(ll[2]) == pytest.approx(_brute_betabinom(ref, n, EPS, RHO), rel=1e-6)
    assert math.exp(ll[1]) == pytest.approx(_brute_present(alt, n, FMIN), rel=1e-4)


def test_zero_overdispersion_is_binomial():
    ll = pg.log_likelihoods(np.array([30]), np.array([2]), 0.01, 0.0, FMIN)[0]
    assert math.exp(ll[0]) == pytest.approx(math.comb(32, 2) * 0.01 ** 2 * 0.99 ** 30, rel=1e-9)


def test_min_vaf_zero_gives_the_uniform_marginal():
    ll = pg.log_likelihoods(np.array([17]), np.array([5]), EPS, RHO, 0.0)[0]
    assert math.exp(ll[1]) == pytest.approx(1 / 23, rel=1e-9)


def test_clear_reference_alternate_and_both_alleles():
    gt, gq, pl = _call(40, 0)
    assert gt == "0/0" and gq >= 20 and pl[0] == 0 and pl[1] == gq
    gt, gq, pl = _call(0, 40)
    assert gt == "1/1" and gq >= 20 and pl[2] == 0
    gt, gq, pl = _call(20, 20)
    assert gt == "0/1" and gq == 99 and pl[1] == 0


def test_subclonal_alt_is_present_not_absent():
    # 8 of 40 fragments (20%): a diploid dosage model would struggle; presence is clear
    gt, gq, _ = _call(32, 8)
    assert gt == "0/1" and gq >= 20


def test_low_depth_is_a_no_call_with_likelihoods_kept():
    gt, gq, pl = _call(3, 0)
    assert gt == "./." and gq < 20 and pl[0] == 0 and pl[1] > 0


def test_zero_depth_is_a_no_call_with_flat_likelihoods():
    # e.g. a site whose only fragments were conflicts (not in AD)
    gt, gq, pl = _call(0, 0)
    assert gt == "./." and pl == (0, 0, 0) and gq == 2      # -10 log10(2/3) = 1.76


def test_ambiguous_low_alt_count_is_a_no_call():
    gt, gq, _ = _call(38, 2)
    assert gt == "./." and gq < 20


def test_large_depth_has_no_underflow():
    for ref, alt in ((5000, 0), (0, 5000), (4000, 1000), (9990, 10)):
        gt, gq, pl = _call(ref, alt)
        assert all(0 <= x < 10 ** 7 for x in pl) and 0 <= gq <= 99
    assert _call(5000, 0)[0] == "0/0" and _call(0, 5000)[0] == "1/1" and _call(4000, 1000)[0] == "0/1"


def test_gq_is_capped_and_pl_is_normalized():
    gt, gq, pl = _call(200, 200)
    assert gq == 99 and min(pl) == 0


def test_vectorized_calls_match_single_calls():
    ref = np.array([40, 0, 20, 3, 0, 36])
    alt = np.array([0, 40, 20, 0, 0, 4])
    gts, gqs, pls = pg.call(ref, alt, EPS, RHO, FMIN, 20)
    for i in range(len(ref)):
        assert (gts[i], int(gqs[i]), tuple(int(x) for x in pls[i])) == _call(ref[i], alt[i])


@pytest.mark.parametrize("kw,word", [({"error_rate": 0.06}, "error_rate"), ({"min_vaf": 0.5}, "min_vaf"),
                                     ({"overdispersion": 1.0}, "overdispersion"),
                                     ({"overdispersion": -0.1}, "overdispersion"),
                                     ({"error_rate": 0.0}, "error_rate")])
def test_invalid_parameters_are_rejected(kw, word):
    with pytest.raises(ValueError, match=word):
        _call(10, 1, **kw)


def test_format_pl_and_gt_tuple():
    assert pg.format_pl((0, 24, 310)) == "0,24,310"
    assert pg.gt_tuple("0/1") == (0, 1) and pg.gt_tuple("1/1") == (1, 1) and pg.gt_tuple("./.") == (None, None)
