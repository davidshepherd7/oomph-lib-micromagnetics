#!/usr/bin/env python3

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

import io
from io import StringIO
import contextlib


def pytest(args, expected_result):
    """Run fpdiff as a python function and compare the result with what we
    expected."""

    print("Running as python function with args", *args)
    print("Expecting", expected_result)

    # Run with output redirected to a string and detailed error information
    # discarded.
    mystdout = StringIO()
    with open(os.devnull, 'w') as null:
        actual_result = fpdiff(*(args + [mystdout, null]))
    mystdout = mystdout.getvalue()

    if "[FAILED]" in mystdout:
        stdout_result = False
    elif "[OK]" in mystdout:
        stdout_result = True
    else:
        # Neither message in output: failed test
        print("Couldn't find [OK] or [FAILED] in output")
        return False

    # Return
    test_pass = (actual_result == expected_result) \
      and (stdout_result == expected_result)

    return test_pass


def shelltest(args, expected_result):
    """Run fpdiff as a shell script and compare the result with what we
    expected."""

    print("Running as shell script with args", *args)
    print("Expecting", expected_result)

    # fpdiff should be accessible from $PATH so we can just call it directly
    fpdiff_path = "fpdiff.py"

    # Run as script and get stdout
    process = subp.Popen([fpdiff_path] + list(map(str, args)),
                            stdout=subp.PIPE, stderr=subp.STDOUT)
    mystdout, _ = process.communicate()
    script_code = process.wait()

    mystdout = mystdout.decode()

    # check
    actual_result = (script_code == 0)

    if "[FAILED]" in mystdout:
        stdout_result = False
    elif "[OK]" in mystdout:
        stdout_result = True
    else:
        # Neither message in output: failed test
        print("Couldn't find [OK] or [FAILED] in output")
        return False

    # Return
    test_pass = (actual_result == expected_result) \
      and (stdout_result == expected_result)

    return test_pass

def main():

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
