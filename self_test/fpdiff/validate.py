#!/usr/bin/env python

# Python 2/3 compatability
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import sys
import os
import fpdiff
import itertools as it
import subprocess as subp

from fpdiff import fpdiff


def pytest(args, expected_result):
    """Run fpdiff as a python function and compare the result with what we
    expected."""

    print("Running as python function with args", *args)
    print("Expecting", expected_result)

    # Run within python
    actual_result = fpdiff(*args)

    # Return true if ok, false otherwise
    return actual_result == expected_result


def shelltest(args, expected_result):
    """Run fpdiff as a shell script and compare the result with what we
    expected."""

    print("Running as shell script with args", *args)
    print("Expecting", expected_result)

    # fpdiff should be accessible from $PATH so we can just call it directly
    fpdiff_path = "fpdiff.py"

    # Run as script
    with open(os.devnull, 'w') as fstdout:
        script_code = subp.call([fpdiff_path] + map(str, args))
    actual_result = (script_code == 0)

    # Return true if ok, false otherwise
    return actual_result == expected_result


def main():

    print("Hello! Watch out in this file: a report of [FAILED] doesn't actually")
    print("mean that the test failed. Since we are testing fpdiff we want to")
    print("check that it fails at the right times and so we *should* see some")
    print("[FAILED] in the output.")

    tests = [
        # Same data
        (["validata/zero", "validata/zero_2", 1, 1e-10], True),

        # Different data
        (["validata/zero", "validata/one", 1, 1e-10], False),

        # Different data, just outside numerical zero
        (["validata/zero", "validata/small", 1, 0.9e-8], False),

        # Different data, just inside numerical zero
        (["validata/zero", "validata/small", 1, 1.1e-8], True),
             ]

    # Run the tests
    passes = list(it.starmap(pytest, tests)) + list(it.starmap(shelltest, tests))

    # Print information on passes
    print("If all of the following are True then the test passed:")
    print(passes)
    print("also you can just check the exit status with 'echo $?'")

    # Check that they passed
    if all(passes):
        return 0
    else:
        return 1


if __name__ == "__main__":
    sys.exit(main())
