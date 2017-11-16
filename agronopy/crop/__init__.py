# -*- coding: utf-8 -*-
"""Crop package for agronopy."""

__author__ = """Jerome Dury"""
__email__ = 'biodiversite@agrobioperigord.fr'
__version__ = '0.1.0'

from agronopy.crop import phenology
from agronopy.crop import grain

__all__ = ['phenology', 'grain']

# def normalize_weight(weight, moisture, moisture_norm):
# """Convert weigth at a given moisture. Instead of weigh, it could be a yield.

# :param weigth: grain weight (kg)
# :param moisture: grain moisture (%)
# :param moisture_norm: standart grain moisture (%, default: 0.15)

# :return weight in the same unit as the input weight
# """
# return (weight * (100 - moisture)) / (100 - moisture_norm)

# def crop_yield(weight, area, moisture, moisture_norm=15):
# """Crop yield (t/ha).

# :param area: given area (m2)
# :param w: grain weight (kg)
# :param moisture: grain moisture (%)
# :param moisture_norm: standart grain moisture (%, default: 0.15)

# :return Normalized crop yield from grain weight per area in t/ha
# """
# return converter.kg2t(normalize_weight(weight, moisture, moisture_norm)) /
# converter.square_meter2ha(area)
