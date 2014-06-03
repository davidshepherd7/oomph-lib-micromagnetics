# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import argparse
import os
import os.path
import scipy as sp
import scipy.optimize

from os.path import join as pjoin

# Make sure *this* versions oomphpy is in the path (before any other
# versions in other places)
sys.path.insert(1, pjoin(os.path.dirname(__file__), "../"))
import oomphpy
import oomphpy.micromagnetics as mm

from fpdiff import fpdiff

def _default_label(data):
    """Get last two parts of the path"""
    b, c = os.path.split(data['-outdir'])
    a, b = os.path.split(b)
    return pjoin(b, c)


def check_error_norm(data, tol=1e-5, identifier=None):
    """Check that max error norm < tol."""

    if identifier is None:
        identifier = _default_label(data)

    max_err_norm = max(map(abs, data['error_norms']))
    norm_test = max_err_norm < tol

    if not norm_test:
        mm.badprint("FAILED in", identifier)
        mm.badprint("with max error norm of", max_err_norm)

    else:
        mm.okprint("error norm ok in", identifier)

    return norm_test


def check_m_length(data, tol=1e-8, identifier=None):
    """Check that abs(|m| - 1) < tol."""

    if identifier is None:
        identifier = _default_label(data)

    max_length_err = max(map(abs, data['m_length_error_means']))
    length_test = max_length_err < tol

    if not length_test:
        mm.badprint("FAILED in ", identifier,
                    "with |m| error of", max_length_err)

    else:
        mm.okprint("|m| max of", max_length_err, "ok in", identifier)

    return length_test


def check_ndt_less_than(data, max_ndt, identifier=None):

    if identifier is None:
        identifier = _default_label(data)

    ndt = len(data['dts'])

    passed = ndt < max_ndt

    if not passed:
        mm.badprint("FAILED in ",  identifier)
        mm.badprint("n steps is", ndt, "which is more than", max_ndt)

    else:
        mm.okprint("n steps ok in", identifier)

    return passed



def check_restarted_in_middle(restart_outdir, restart_point=20):

    last_no_soln = "soln%i.dat" % (restart_point - 1)
    first_soln = "soln%i.dat" % restart_point

    ok = (not os.path.isfile(pjoin(restart_outdir, last_no_soln))) \
      and (os.path.isfile(pjoin(restart_outdir, first_soln)))

    if ok:
        mm.okprint("Restarted at correct point in", restart_outdir)
    else:
        mm.badprint("Failed or restarted at wrong point in", restart_outdir)

    return ok


def check_solns_match(main_dir="Validation", compare_dir="validata", **kwargs):

    matches = []
    for fname in os.listdir(main_dir):
        if ("soln" in fname) and (".dat" in fname):
            match, _, _ = fpdiff(pjoin(main_dir, fname),
                                 pjoin(compare_dir, fname), **kwargs)
            matches.append(match)

    if len(matches) == 0:
        mm.badprint("No files in", main_dir)
        return False

    ok = all(matches)
    if ok:
        mm.okprint("Files match in", main_dir, compare_dir)
    else:
        mm.badprint("Files don't match in", main_dir, compare_dir)

    return ok


def check_mean_m_matches(data1, data2, tol=1e-4):


    id1 = _default_label(data1)
    id2 = _default_label(data2)

    for a, b in zip(data1['mean_mxs'], data2['mean_mxs']):
        if abs(a - b) > tol:
            mm.badprint("Failed mx comparison in", id1, id2, "error of", abs(a - b))
            return False

    for a, b in zip(data1['mean_mys'], data2['mean_mys']):
        if abs(a - b) > tol:
            mm.badprint("Failed my comparison in", id1, id2, "error of", abs(a - b))
            return False

    for a, b in zip(data1['mean_mzs'], data2['mean_mzs']):
        if abs(a - b) > tol:
            mm.badprint("Failed mz comparison in", id1, id2, "error of", abs(a - b))
            return False

    mm.okprint("m comparison ok in", id1, id2)
    return True


def check_solns_match_key(datasets, key_to_compare_over, **kwargs):
    """Given a dict of solutions split into pairs which only differ by
    key_to_compare_over then compare the pairs.

    kwargs go into underlying checking function
    """

    pairs = mm.split_to_comparable_groups(datasets, key_to_compare_over)

    ok = []
    for p in pairs:

        # Check we have pairs
        assert len(p) == 2

        # Compare
        ok.append(check_mean_m_matches(p[0], p[1], **kwargs))

    return all(ok)


def check_dicts_match(data, exact, keys, **kwargs):
    """Given two dicts and a list of keys check that those values match.

    ??ds actually temp hack: use the second value of the array given by
    data['k']...
    """

    def fpcheck(v1, v2, rtol=1e-8, zero=1e-14):
        """Check two floats are the same: if both small then true, otherwise check
        relative difference.
        """
        if v1 < zero and v2 < zero:
            return True
        else:
            return (abs((v1 - v2)/v1) < rtol)


    dataid = _default_label(data)

    ok = []
    for k in keys:
        try:
            a = fpcheck(exact[k], data[k][1], **kwargs)
            ok.append(a)

            # Print messages
            if a:
                mm.okprint(k, "comparison ok in", dataid)
            else:
                mm.badprint("Failed", k, "comparison with values of",
                            data[k][1], exact[k], "in", dataid)

        except KeyError:
            mm.badprint("Failed: missing key", k, "in", dataid)
            ok.append(False)

    return all(ok)


def check_convergence(datasets, expected_rate, tol=0.2):
    """For a set of convergence test check that the rate of convergence as dt
    and h -> 0 is as expected. dt and h are linked so it is enough to
    test only for dt.
    """

    # Check that we have a set of convergence tests
    assert all([d['-convergence-test'] == "1" for d in datasets])

    # Extract data to be fitted
    dts, errs = mm.unzip([(sp.mean(d['dts']), max(d['error_norms'])) for d in datasets])

    # fit the data to a 2nd order polynomial
    def f(x, a, b, c):
        return a*(x**b) + c
    parameters, covar = scipy.optimize.curve_fit(f, dts, errs)
    rate = parameters[1]     # rate of convergence

    # Check if near to expected value
    ok = abs(rate - expected_rate) < tol

    # Print messages, identify by first dataset for simplicity.
    dataid = _default_label(datasets[0])
    if ok:
        mm.okprint("convergence rate", rate, "ok in", dataid)
    else:
        mm.badprint("convergence rate", rate,
                    "but expected", expected_rate, "+/-", tol,
                    "in", dataid)

    # ??ds check goodness of fit?

    return ok
