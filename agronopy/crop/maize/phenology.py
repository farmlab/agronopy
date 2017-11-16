# -*- coding: utf-8 -*-
from agronopy.crop import phenology

# Constant usually used for the maize
BASE_MIN = 6
BASE_MAX = 30


def growing_degree_day(tmin, tmax, base_min=BASE_MIN, base_max=BASE_MAX):
    """Growing degree days (GDD).

        GDD, also called growing degree units (GDUs), are a heuristic tool in
        phenology and are a measure of heat accumulation used by agricultural
        scientist and farmers to predict plant growth stage.

        :param tmin: minimum temperature of the day (°C)
        :param tmax: maximun temperature of the day (°C)
        :param base_min: base temperature (°C)
        :param base_max: maximum base temperature (°C)

        :return: Growing degree days (°C)
    """
    return phenology.growing_degree_day(tmin, tmax, BASE_MIN, BASE_MAX)
