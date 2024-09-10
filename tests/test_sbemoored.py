#!/usr/bin/env python

"""Tests for `sbemoored` package."""

import pytest

import sbemoored as sbe


@pytest.fixture
def testfile_sbe56(rootdir):
    # note: rootdir has been defined as fixture in conftest.py
    file = rootdir.joinpath("data/SBE05600422_2020-10-08.csv")
    return file


@pytest.fixture
def testfile_sbe37(rootdir):
    # note: rootdir has been defined as fixture in conftest.py
    file = rootdir.joinpath("data/SBE37_12712_NISKINE.cnv")
    return file


@pytest.fixture
def testfile_sbe37_no_time(rootdir):
    # note: rootdir has been defined as fixture in conftest.py
    file = rootdir.joinpath("data/SBE37_2864_NISKINE.cnv")
    return file


@pytest.fixture
def header(testfile_sbe37):
    return sbe.sbe37.read_header(testfile_sbe37)


def test_locate_example_data(testfile_sbe37, testfile_sbe37_no_time):
    assert testfile_sbe37.exists()
    assert testfile_sbe37_no_time.exists()


def test_read_header(testfile_sbe37):
    header = sbe.sbe37.read_header(testfile_sbe37)
    assert len(header) == 291


def test_parse_header(header):
    ph = sbe.sbe37.parse_header(header)
    assert ph["sn"] == 12712


def test_read_testfile_sbe56(testfile_sbe56):
    ds = sbe.sbe56.read_csv(testfile_sbe56)
    assert ds.SN == 5600422
    assert ds[-1] == 20.1797


def test_read_testfile_sbe37(testfile_sbe37):
    ds = sbe.sbe37.read_sbe_cnv(testfile_sbe37)
    assert ds.SN == 12712
    assert ds.c[-1] == 32.54289
    assert ds.p[-1] == 2891.96


def test_read_testfile_sbe37_no_time(testfile_sbe37_no_time):
    ds = sbe.sbe37.read_sbe_cnv(testfile_sbe37_no_time)
    assert ds.SN == 2864
    assert ds.c[-1] == 35.9654
    assert ds.p[-1] == 208.475


def test_read_xmlcon(testfile_sbe37):
    ds = sbe.sbe37.proc(testfile_sbe37)
    assert ds.p.CalDate == "03-Jan-18"
    assert ds.p.SN == 2154084
