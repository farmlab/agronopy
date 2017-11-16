#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for `agronopy` package."""

import unittest

from agronopy import crop
from agronopy.crop import maize


class TestAgronopy(unittest.TestCase):
    """Tests for `agronopy` package."""

    def setUp(self):
        """Set up test fixtures, if any."""
        # case where ggd > base_min and < base_max
        ini = {
            'tmin': 0,
            'tmax': 20,
            'base_min': 5,
            'base_max': 30,
            'result': 5,
            'result_maize': 4
        }
        # case where gdd < base_min
        bmin = ini
        bmin['tmax'] = 8
        bmin['result'] = 0
        bmin['result_maize'] = 0
        # case where gdd > base_max
        # case where gdd < base_max for the maize
        bmax = ini
        bmax['tmin'] = 30
        bmax['tmax'] = 40
        bmax['result'] = 0
        bmin['result_maize'] = 29
        self.gdd_data = [ini, bmin, bmax]

    def tearDown(self):
        """Tear down test fixtures, if any."""

    def test_phenology_growing_degree_day(self):
        """Test something."""
        for data in self.gdd_data:
            crop_res = crop.phenology.growing_degree_day(
                data['tmin'], data['tmax'], data['base_min'], data['base_max'])
            self.assertEqual(crop_res, data['result'])

            crop_res = maize.phenology.growing_degree_day(
                data['tmin'],
                data['tmax'], )
            self.assertEqual(crop_res, data['result_maize'])
