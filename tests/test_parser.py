#!/usr/bin/env python
##############################################################################
#
# (c) 2025 Simon Billinge.
# All rights reserved.
#
# File coded by: Caden Myers, Simon Billinge, and members of the Billinge
#                group.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.cmipdf/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Tests for pdf package."""


import numpy as np
import pytest

from diffpy.srfit.fitbase import ProfileParser

# ----------------------------------------------------------------------------


# The tests in this file are for ProfileParser which belongs to diffpy.srfit,
# but since it is used here, we will test it here on test data with modern
# diffpy format.
def testParser1(datafile):
    filename = datafile("ni-q27r100-neutron.gr")
    parser = ProfileParser()
    parser.parse_file(filename)

    meta = parser._meta

    assert str(filename) == meta["filename"]
    assert 1 == meta["nbanks"]
    assert "N" == meta["stype"]
    assert 27 == meta["qmax"]
    assert 300 == meta.get("temperature")
    assert meta.get("qdamp") is None
    assert meta.get("qbroad") is None
    assert meta.get("spdiameter") is None
    assert meta.get("scale") is None
    assert meta.get("doping") is None

    x, y, dx, dy = parser.get_data()
    assert dx.tolist() == len(x) * [0]
    assert dy.tolist() == len(x) * [0]

    testx = np.linspace(0.01, 100, 10000)
    diff = testx - x
    res = np.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testy = np.array(
        [
            1.144,
            2.258,
            3.312,
            4.279,
            5.135,
            5.862,
            6.445,
            6.875,
            7.150,
            7.272,
        ]
    )
    diff = testy - y[:10]
    res = np.dot(diff, diff)
    assert 0 == pytest.approx(res)

    return


def testParser2(datafile):
    data = datafile("si-q27r60-xray.gr")
    parser = ProfileParser()
    parser.parse_file(data)

    meta = parser._meta

    assert str(data) == meta["filename"]
    assert 1 == meta["nbanks"]
    assert "X" == meta["stype"]
    assert 27 == meta["qmax"]
    assert 300 == meta.get("temperature")
    assert meta.get("qdamp") is None
    assert meta.get("qbroad") is None
    assert meta.get("spdiameter") is None
    assert meta.get("scale") is None
    assert meta.get("doping") is None

    x, y, dx, dy = parser.get_data()
    testx = np.linspace(0.01, 60, 5999, endpoint=False)
    diff = testx - x
    res = np.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testy = np.array(
        [
            0.1105784,
            0.2199684,
            0.3270088,
            0.4305913,
            0.5296853,
            0.6233606,
            0.7108060,
            0.7913456,
            0.8644501,
            0.9297440,
        ]
    )
    diff = testy - y[:10]
    res = np.dot(diff, diff)
    assert 0 == pytest.approx(res)

    testdy = np.array(
        [
            0.001802192,
            0.003521449,
            0.005079115,
            0.006404892,
            0.007440527,
            0.008142955,
            0.008486813,
            0.008466340,
            0.008096858,
            0.007416456,
        ]
    )
    diff = testdy - dy[:10]
    res = np.dot(diff, diff)
    assert 0 == pytest.approx(res)

    assert dx.tolist() == [0] * len(dx)
    return
