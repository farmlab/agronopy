# -*- coding: utf-8 -*-


def growing_degree_day(tmin, tmax, base_min, base_max):
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
    gdd = (tmin + tmax) / 2 - base_min
    if gdd > 0 and gdd < base_max:
        return gdd
    return 0
