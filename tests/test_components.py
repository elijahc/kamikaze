#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for components module"""

import pytest

from kamikaze.components import CassetteFactory

@pytest.fixture
def fake_params():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    test_params = dict(
            filename='/Users/elijahc/dev/kamikaze/kamikaze/data/kras_mrna_va.fa',
            region=slice(192,762),
            lib_name='KRAS_mutagenesis',
            )
    return test_params

@pytest.fixture
def dummy_cf(fake_params):
    return CassetteFactory(**fake_params)

@pytest.fixture
def dummy_target():
    return slice(192,195)

def test_cassette_factory(fake_params):
    import Bio
    cf = CassetteFactory(**fake_params)

    # Test correct setting of attributes
    for k in fake_params.keys():
        assert k in vars(cf)
    assert fake_params['filename'] is cf.filename
    assert cf.region is fake_params['region']
    assert cf.lib_name is fake_params['lib_name']
    assert isinstance(cf.reference, Bio.SeqRecord.SeqRecord)

def test_find_pam_sites(dummy_cf):
    cf = dummy_cf
    assert 'find_pam_sites' in dir(cf)
    pam_sites = cf.find_pam_sites()
    assert isinstance(pam_sites,list)

def test_nearest_pam_site(dummy_cf):
    cf = dummy_cf
    assert 'nearest_pam_site' in dir(cf)
    with pytest.raises(AssertionError):
        cf.nearest_pam_site(120)

    k_sites = cf.nearest_pam_site(slice(192,195),n=3)
    assert isinstance(k_sites,list)
    assert len(k_sites) is 3

def test_build_slug(dummy_cf,dummy_target):
    from kamikaze.parts import Slug
    cf = dummy_cf
    tgt = dummy_target
    assert 'build_slug' in dir(cf)
    single_slug = cf.build_slug(tgt)
    assert isinstance(single_slug,list)
    assert len(single_slug) is 1

    slugs = cf.build_slug(tgt,n=3)
    assert len(slugs)<=3


def test_build_payload():
    pass

def test_build_cassette():
    pass