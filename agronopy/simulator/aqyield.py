"""File
File: aqyield.py
Author: Jerome Dury
Email: jerome.dury@flyingheep.fr
Github: https://github.com/farmlab
Description:
http://maelia-platform.inra.fr/modeles/processus-agricoles/dynamique-sol-culture-2/dynamique-sol-culture/

Reference: Constantin, J. , Willaume, M., Murgue, C., Lacroix, B., Therond, O.
(2015). The soil-crop models STICS and AqYield predict yield and soil water
content for irrigated crops equally well with limited data. Agricultural and
Forest Meteorology (206), 55-68.
"""
import datetime
import pandas as pd
import math

from bokeh.plotting import figure
from bokeh.models.formatters import DatetimeTickFormatter
from bokeh.models import Range1d, LinearAxis

from .simulator import SimulatorMixin
from agronopy.utils.plot import COLOR_CHART


class Simulator(SimulatorMixin):
    # default_variables = {
    # ITK
    # "sowing_day" : "1996-04-15",  #15/4
    # "delta_sowing": 0,
    # "ratio_h2o": 0,
    # "ratio_n": 1,
    # }

    # default_params =  {
    # Irrigation
    # "irr_max": 45,
    # "sirr1": 0.8,
    # "sirr2": 0.9,
    # "sirr3": 0.7,

    # Crop MAIZE
    # "gdd_f": 2000,
    # "c_fm": 1.55,

    # "c_a": 3,
    # "c_root": 10,
    # "c_vig": 0.5,
    # "c_ratio": 0.5,
    # "kc_max": 1.2,

    # Soil
    # "clay": 28,
    # "ru_ini" : 50,
    # "ru_max" : 180,
    # "ru_secu" : 0.8,
    # "c_eiv": 0.9,
    # "c_eva1": 5,
    # "c_eva2": 5,
    # "layer1_depth": 5,
    # "res_perm": 30
    # }

    # def __init__(self, inputs, variables={}, params={}):
    # SimulatorMixin.__init__(self, inputs, variables, params)
    # self.sowing_day = pd.to_datetime(self.variables["sowing_day"])
    # self.harvesting_day = None

    def run(self):
        ts = self.inputs["eday"]
        tmean = self.inputs["tmean"]
        rain = self.inputs["rain"]
        pet = self.inputs["pet"]

        # Static variable
        # Sowing day
        sowing_day = pd.to_datetime(
            self.variables["sowing_day"]) - pd.datetime(1970, 1, 1)
        # growing degree day at harvesting: trigger harvesting
        gdd_h = self.params["gdd_f"] * self.params["c_fm"]
        # irrigation threshold
        s1 = self.S1(self.params["sirr1"], self.variables["ratio_h2o"])
        s2 = self.S2(self.params["sirr2"], self.variables["ratio_h2o"])
        s3 = self.S3(self.params["sirr3"], self.variables["ratio_h2o"])
        ps2 = 0.4
        ps3 = 0.8
        ps4 = 1.1

        lds = self.LimDrySoil(self.params["layer1_depth"], self.params["clay"])

        for i, d in enumerate(ts):
            # echT
            if i == 0:
                sown = [False]
                echT = [0]
                echV = [0]
                Kc0 = [0]
                Kc = [0]
                Eva = [0]
                RUr = [self.params["ru_ini"]]
                RUt = [self.params["ru_max"] / 2]
                Hs = [0]
                Hw = [self.params["ru_ini"] / 2]
                Hr = [self.params["ru_ini"] / 2]
                sirr = [0]
                doseIrr = [0]
                irr = [0]
                pir2 = [0]
                pir1 = [0]

                TRm = [0]
                TRr = [0]
                TRw = [0]
                Drainage = [0]
                TM = [1]
            else:
                echT.append(echT[i - 1] + tmean[i])
                sown.append(
                    self.is_sown(ts[i], sowing_day, echV[i - 1], gdd_h))
                # Cumulative growing degree day after sowing
                # Cumulative growing degree day after sowing
                echV.append(
                    self.GddV(ts[i], tmean[i], echV[i - 1], gdd_h, sown[i - 1],
                              self.params["c_a"]))
                # RUr: size of the soil layer accessible by the root
                RUr.append(
                    self.RUr(echV[i], self.params["gdd_f"], self.params[
                        "ru_ini"], self.params["ru_max"], self.params[
                            "c_root"], sown[i]))
                # Evaporation
                Eva.append(
                    self.Evaporation(pet[i], Kc0[i - 1], Hw[i - 1], Hs[
                        i - 1], self.params["ru_ini"], self.params[
                            "layer1_depth"], self.params["c_eiv"], self.params[
                                "c_eva1"], self.params["c_eva2"]))
                # Transpiration maxi
                TRm.append(self.Tm(pet[i], Eva[i], Kc[i - 1]))
                # Actual transpiration
                TRr.append(
                    self.Ta(TRm[i], Hr[i - 1], Hs[i - 1], RUr[i - 1],
                            self.params["clay"]))
                # Water satisfaction rate
                TRw.append(self.Tw(TRr[i], Hw[i - 1], Hr[i - 1]))

                TM.append(self.SR_h20m(TRr[i], TRm[i]))

                # Irrigation
                # Trigger irrigation
                sirr.append(
                    self.Sirr(sown[i], echV[i], self.params["gdd_f"], s1, s2,
                              s3, ps2, ps3, ps4, self.params["c_fm"]))
                # Required Dose
                doseIrr.append(
                    self.doseIrr(echV[i], echV[i - 1], RUr[i - 1], Hr[
                        i - 1], self.params['ru_secu'], self.params["gdd_f"],
                                 self.params['c_fm'], sirr[i], TM[i]))

                # Given Dose
                irr.append(self.Irr(self.params["irr_max"], doseIrr[i]))

                # Rain + irrigation
                pir2.append(
                    self.RainIrr2(pir2[i - 1], irr[i], rain[i], self.params[
                        "res_perm"]))
                pir1.append(
                    self.RainIrr1(pir2[i - 1], irr[i], rain[i], pir2[i]))

                # KC
                Kc0.append(
                    self.Kc(echV[i], self.params["gdd_f"], self.params[
                        "c_fm"], self.params["c_vig"], tmean[i], self.params[
                            "kc_max"], Kc0[i - 1], self.variables["ratio_n"]))
                Kc.append(
                    self.Kca(
                        sum(TRm), sum(TRr), Kc0[i], self.params["c_ratio"]))

                # Soil Humidity
                RUt.append(
                    self.RUt(self.params["ru_max"], RUt[i - 1], pir1[i], Eva[
                        i], TRr[i]))
                Hs.append(
                    self.Hs(self.params["layer1_depth"], lds, Hs[i - 1], pir1[
                        i], Eva[i]))
                Hw.append(
                    self.Hw(Hw[i - 1], pir1[i], Eva[i], TRw[i], self.params[
                        "ru_ini"]))
                Hr.append(
                    self.Hr(RUr[i], RUr[i - 1], self.params["ru_max"], pir1[i],
                            RUt[i - 1], Hr[i - 1], Eva[i], TRr[i]))

                Drainage.append(
                    self.Drainage(RUt[i - 1], pir1[i], Eva[i], TRr[i],
                                  self.params["ru_max"]))

        print(sum(doseIrr))

        self.res = {
            "date": ts,
            "echT": echT,
            "echV": echV,
            "tmean": tmean,
            "etp": pet,
            "rain": rain,
            "irr": irr,
            "sirr": sirr,
            "pir2": pir2,
            "pir1": pir1,
            "RUr": RUr,
            "Kc0": Kc0,
            "Kc": Kc,
            "TRm": TRm,
            "TRr": TRr,
            "Hs": Hs,
            "Hw": Hw,
            "Hr": Hr,
            "RUt": RUt,
            "Eva": Eva,
            "TM": TM,
            "Drainage": Drainage
        }

    def GddV(self, day, t_moy, gddv, gdd_h, is_sown, c_a):
        """ The function return the cumulative °D from sowing including
        chilling requirement coefficient.

            day: number of day in the year
            sowing: boolean indicatigin sowing
            t_moy: daily mean temperature
            is_sown: boolean indicating sowing
            gdd_h: Growing degree day at harversting
            c_a:  coeficient that represent the chilling seed requirements
            to trigger their subsequent emergence and growth
            (coéficient d'alternativité)
        """
        if is_sown is True:
            # before harvesting
            if gddv < gdd_h:
                doy = (datetime.datetime(1970, 1, 1, 0, 0) +
                       day).timetuple().tm_yday
                if doy > 310 or doy <= 45:
                    gddv += t_moy * c_a / 10
                else:
                    gddv += t_moy
        else:
            gddv = 0
        return gddv

    def GddT(self, sowing, echt_t1, t_moy):
        """ Cumulative °D from sowing

            sowing: boolean indicatigin sowing
            echt_t1: cumulative °D from sowing at t-1
            t_moy: daily mean temperature
        """
        if sowing is True:
            gddt = echt_t1 + t_moy
        else:
            gddt = 0
        return gddt

    # col L
    def RUr(self, gddv, gdd_f, ru_rini, ru_max, c_root, sowing):
        """ The function return the size of soil layer accessible by the root

            gddv: °D from sowing including chilling requirement coeffiecient
            gdd_f: growing degree day at flowering
            ru_rini: Initial available water content accessible by the root
            (mm)
            c_root: root growth coeficient
        """
        if sowing is False:
            ru_r = ru_rini
        else:
            if gddv < gdd_f:
                ru_r = max(ru_rini, min(ru_max, gddv / c_root))
            else:
                ru_r = max(ru_rini,
                           min(ru_max,
                               gdd_f / c_root + (gddv - gdd_f) / c_root / 2))
        return ru_r

    # col M
    def Kc(self, gddv, gdd_f, c_fm, c_vig, t_moy, kc_max, kc_0, ratio_n):
        """ kc, crop coeficient

            gddv: °D from sowing including chilling requirement coeffiecient
            gdd_f: growing degree day at flowering
            c_fm: coeficient flowering / harvesting
            c_vig: vigour coeficient
            kc_0: crop coeficient initial
            kc_max: crop coeficient maximum
            t_moy: daily mean temperature
            sowing: boolean indicatigin sowing

        """
        kc_max_flo = kc_max * (1 - (1.2 - ratio_n)**2.0)
        if gddv < gdd_f * c_fm:
            if gddv < gdd_f * c_vig:
                kc = kc_0 + 2 * kc_max_flo * (t_moy / gdd_f) * gddv / (
                    gdd_f * c_vig)
                return kc
            elif gddv < gdd_f:
                return kc_0 + 2 * kc_max_flo * (t_moy / gdd_f) * (
                    1 - gddv / gdd_f) / (1 - c_vig)
            else:
                return kc_0 + 2 * kc_max_flo * (t_moy / gdd_f) * (
                    1 - gddv / gdd_f)
        else:
            return 0

    # col N
    def Kca(self, tm_cumul, ta_cumul, kc0, c_ratio):
        """ Calculate the actual crop coeficient including water stress

            tm_cumul: cumulative maximal transpiration
            ta_cumul: cumulative actual transpiration
            kc: crop coeficient
            c_ratio: sensibility to water stress coeficient
        """
        if tm_cumul > 0:
            kc_actual = kc0 * (ta_cumul / tm_cumul)**c_ratio
        else:
            kc_actual = kc0
        return kc_actual

    # col 0
    def Tm(self, etp, eva, kc_t1):
        """ Maximal transpiration

            etp:  potential evapotranspiration (mm)
            kc_t1: crop coeficient at day-1
        """
        return (etp - eva) * kc_t1

    # col V
    def Tw(self, ta, hw, hr):
        """ Transpiration of the 30 cm horizon depanding on the hydric state
        and the root profile

            ta: actual transpiration
            hw: water content on the 30 cm layer horizon (mm)
            hr: water content of the root layer (mm)
        """
        return ta * hw / hr

    # col W
    def Ta(self, tm, hr, hs, rur_t1, clay):
        """ actual transpiration

            tm: Maximal transpiration
            hr: water content of the root layer (mm)
            hs: Shallow water reserve
            ru_rt1: size of soil layer accessible by the root at t-1
            clay: soil clay content (%)
        """
        return tm * (1 - abs(1 - (hr - hs) / (rur_t1 - 5))**(120 /
                                                             (clay + 25)))

    # COL X
    def SR_h20m(self, ta, tm):
        """ water satisfaction ratio

            tm: Maximal evapotranspiration
            ta: actual transpiration
        """
        return (ta + 0.1) / (tm + 0.1)

    def Yield(self, a, b, c, ta_cumul, tm_cumul, ratio_n, precocity, c_prec,
              y_pot):
        """ Yield function
            The yield is estimated based on a potential yield reduced by water
            nd nitrogen stress.

            a, b, c: parameter of the water function
            ta_cumul: cumulative actual transpiration
            tm_cumul: cumulative maximal evapotranspiration
            ratio_n: Nitrogen restriction coeficient (between 0 and 1)
            c_prec: coef used to represent different crop variety (precocity)
            y_pot: potential yield for mid variety
        """
        # H20
        sr_h2o_tot = ta_cumul / (tm_cumul + 0.05)
        reduc_h2o = a * sr_h2o_tot**2 + b * sr_h2o_tot + c

        # N
        if ratio_n > 1.125:
            reduc_n = 0
        else:
            reduc_n = 1.02 - (1 - 0.88 * ratio_n)**1.8

        return reduc_h2o * c_prec * reduc_n * y_pot

    # ---------------------------------------------------------
    # FARMER
    # ---------------------------------------------------------
    # def SowingDay(sowing_day, delta_sowing):
    # """ Date of sowing

    # sowing_day: date of sowing (day)
    # delta_sowing: number of day to shift sowing day (day)
    # """
    # return sowing_day + delta_sowing

    def GDDh(self, gddv, echt, gdd_f, gdd_h, sowing, c_fh):
        """ Cumulative degree at harvesting

            gddv: °D from sowing including chilling requirement coeffiecient
            echt: Cumulative °D from sowing
            gdd_f: growing degree day at flowering
            gdd_h: growing degree day at harvesting
            sowing: boolean indicatigin sowing
            c_fh: coeficient flowering / harvesting
        """
        if sowing is True:
            if gddv >= gdd_f and gdd_h is None:
                gdd_h = echt + gdd_f * (c_fh - 1)
            else:
                gdd_h = gdd_h
        else:
            gdd_h = None
        return gdd_h

    def is_sown(self, day, sowing_day, echv, gdd_h):
        """ Return true between sowing and harvesting

            day: day of the year
            sowing_day: date of sowing (day)
            echt: Cumulative °D from sowing
            gdd_h: growing degree day at harvesting
        """
        if day >= sowing_day:
            return True
        return False

    def is_growing(self, day, sowing_day, echv, gdd_h):
        """ Return true between sowing and harvesting

            day: day of the year
            sowing_day: date of sowing (day)
            echt: Cumulative °D from sowing
            gdd_h: growing degree day at harvesting
        """
        if day >= sowing_day:
            if echv > gdd_h:
                return False
            return True
        # if sowing_day > 305: # winter crops
        # if day == sowing_day:
        # sowing = True
        # elif echt > gdd_h:
        # sowing = False
        # else: # summer crops
        # if day >= sowing_day and day < 305 and echt <= gdd_h:
        # sowing = True
        # else:
        # sowing = False
        return False

    # Irrigation
    def S1(self, sirr_1, ratio_h2o):
        return sirr_1 * ratio_h2o

    def S2(self, sirr_2, ratio_h2o):
        return min(0.95, sirr_2 * ratio_h2o)

    def S3(self, sirr_3, ratio_h2o):
        return min(0.95, sirr_3 * ratio_h2o)

    def Sirr(self, sowing, gddv, gdd_f, s1, s2, s3, ps2, ps3, ps4, c_fm):
        """ Trigger irrigation

            sowing_day: date of sowing (day)
            gddv: °D from sowing including chilling requirement coeffiecient
            gdd_f: growing degree day at flowering
            s1, s2, s3, s4: Irrigation threshold
            p1, p2, p3, p4: Irrigation threshold changing phase
            c_fm: coeficient flowering / maturity
        """
        if sowing is False:
            return 0.0
        else:
            if gddv < gdd_f * 0.4:  # before 0.4 flowering
                return s1
            elif gddv < gdd_f * 0.8:  # before 0.8 flowering
                return s1 + (s2 - s1) * (gddv - gdd_f * ps2) / (gdd_f *
                                                                (ps3 - ps2))
            elif gddv < gdd_f * 1.1:  # after flowering
                return s2
            elif gddv < gdd_f * 1.5:  # after flowering
                return s2 + (s3 - s2) * (gddv - gdd_f * 1.1) / (gdd_f *
                                                                (c_fm - ps4))
        return 0.0

    def doseIrr(self, gdd_v, gddv_t1, ru_r, hr, c_securu, gdd_f, c_fm, sirr,
                tr_m):
        """
            Irrigation water dose

            ru_r: size of soil layer accessible by the root
            hr: water content of the root layer (mm)
            c_securu: security coeficient pour ru_r
            gdd_f: growing degree day at flowering
            c_fm: coeficient flowering / maturity
            gddv: °D from sowing including chilling requirement coeffiecient
            gddv_t1: °D from sowing including chilling requirement coeffiecient
            at t-1
            sirr: trigger irrigation
            tr_m: water satisfaction ratio

        """
        if gddv_t1 > 0 and gddv_t1 < gdd_f * 1.5 and tr_m < sirr:
            dose_irr = min((ru_r - hr) * c_securu,
                           max(0, (gdd_f * c_fm - gdd_v) * 3 / 25))
        else:
            dose_irr = 0
        return dose_irr

    # ---------------------------------------------------------
    # ITK
    # --------------------------------------------------------

    def Irr(self, dose_water, dose_irr):
        """ water dose

              dose_water: maximum irrigation water supply
              dose_irr: calculated irrigation water supply
          """
        return min(dose_water, dose_irr)

    # ---------------------------------------------------------
    # SOIL
    # --------------------------------------------------------
    def Evaporation(self, etp, kc0_t1, hw_t1, hs_t1, ru_rini, l_s, c_eiv,
                    c_eva1, c_eva2):
        """ Soil evaporation

            etp:  potential evapotranspiration (mm)
            c_eiv: Coeficient that reduces the evaporation
            kc: crop coefficient
            hw: water content on the 30 cm layer horizon (mm)
            hs: Shallow water reserve
            ru_rini: Initial available water content accessible by the root
            (mm)
            l_s: Shallow layer depth (mm)
            c_eva1, c_eva2: Evaporation coeficient

        """
        return etp * max(0.0, c_eiv - kc0_t1) * (
            0.3 * (hw_t1 / ru_rini)**c_eva2 + 0.7 * (hs_t1 / l_s)**c_eva1)

    # Water distribution between RainIrr1 (directly available) and
    # RainIrr2 (waiting water excess)
    def RainIrr2(self, rain_irr_2_t1, irr, rain, resperm):
        """ Water (rain + irr) not directly available

            rain_irr_2_t1: water not directly available at t-1
            irr: daily irrigation (mm)
            rain: daily rain (mm)
            resperm: reservoir permeability (mm)
        """
        return max(rain_irr_2_t1 + irr + rain - resperm, 0.0)

    def RainIrr1(self, rain_irr_2_t1, irr, rain, rain_irr_2):
        """ Water (rain + irr) directly available

            rain_irr_1_t1: water directly available at t-1
            rain_irr_2_t1: water not directly available at t-1
            irr: daily irrigation (mm)
            rain: daily rain (mm)
        """
        return rain + irr + rain_irr_2_t1 - rain_irr_2

    def LimDrySoil(self, l_s, clay):
        """
        Surface condition

            l_s: surface layer depth (mm)
            clay: soil clay content (%)
        """
        return l_s - l_s / (1 / (1 + 0.02 * clay))

    def Hs(self, l_s, lim_dry_soil, hs, rain_irr_1, eva):
        """  Shallow water reserve

            l_s: surface layer depth (mm)
            lim_dry_soil: Surface condition
            rain_irr_1: water directly available
            eva: soil evaporation

        """
        return min(l_s, max(lim_dry_soil, hs + rain_irr_1 - eva))

    def Hw(self, hw, rain_irr_1, eva, ta_w, ru_rini):
        """
            Water content on the 30 cm layer horizon (mm)

            hw: water content on the 30 cm layer horizon (mm)
            rain_irr_1: water directly available
            eva: soil evaporation
            ta_w: Transpiration of the 30 cm horizon
            ru_rini: Initial soil layer size accessible by the root (mm)
        """
        return min(hw + rain_irr_1 - eva - ta_w, ru_rini)

    def Hr(self, ru_r, ru_rt1, ru_max, rain_irr_1, rut_t1, hr_t1, eva, ta):
        """
            Water content of the root layer (mm)

            ru_r: size of soil layer accessible by the root
            ru_r_t1: size of soil layer accessible by the root at t-1
            ru_max: Maximum available water content (mm)
            rain_irr_1: water directly available
            ru_t: total water content (mm)
            hr_t1:  Water content of the root layer (mm) at t-1
            eva: soil evaporation
            ta: actual transpiration
        """
        a = hr_t1 + rain_irr_1 - eva - ta
        if (ru_rt1 < ru_max):
            b = (ru_r - ru_rt1) * (rut_t1 - hr_t1) / (ru_max - ru_rt1)
            hr = min(ru_r, a + b)
        else:
            hr = min(ru_r, a)
        return hr

    def RUt(self, ru_max, ru_t_t1, rain_irr_1, eva, ta):
        """
            Total water content (mm)

            ru_max: Maximum available water content (mm)
            ru_t_t1: total water content (mm)
            rain_irr_1: water directly available
            eva: soil evaporation
            ta: actual transpiration
        """
        return min(ru_max, ru_t_t1 + rain_irr_1 - eva - ta)

    def Drainage(self, ru_t_t1, rain_irr_1, eva, ta, ru_max):
        """
            Water loss by drainage

            ru_t_t1: total water content (mm) at t-1
            rain_irr_1: water directly available
            eva: soil evaporation
            ta: actual transpiration
            ru_max: Maximum available water content (mm)

        """
        return max(0, ru_t_t1 + rain_irr_1 - eva - ta - ru_max)

    def plot(self):

        x = [
            pd.to_datetime(d + pd.datetime(1970, 1, 1))
            for d in self.res["date"]
        ]
        y_rut = self.res["RUt"]
        y_kc = self.res["Kc"]
        y_tm = self.res["TM"]

        colors = COLOR_CHART["agp"]

        p = figure(toolbar_location="above", plot_width=900, plot_height=500)
        p.xaxis.formatter = DatetimeTickFormatter(months=['%m-%Y'])

        y_max = math.ceil(self.params["ru_max"] / 50) * 50
        p.y_range = Range1d(start=0, end=y_max)

        y_max_extra = max(max(y_kc), 1)
        p.extra_y_ranges = {"y2": Range1d(start=0, end=y_max_extra)}
        p.xaxis.axis_label = "Date"

        # croissance du maïs
        x_env = x + x[::-1]
        kc_env = y_kc + [0] * len(y_kc)
        p.patch(
            x_env,
            kc_env,
            y_range_name="y2",
            color=colors["green"][1],
            fill_alpha=0.5)

        #  taux de satisfaction en eau  du maïs
        p.add_layout(
            LinearAxis(y_range_name="y2", axis_label="ratio"), 'right')

        tm_env = y_tm + [1] * len(y_tm)

        p.patch(
            x_env,
            tm_env,
            y_range_name="y2",
            color=colors["orange"][1],
            fill_alpha=0.5)
        p.line(x, y_tm, y_range_name="y2", color=colors["purple"][0])

        # Eau dans le sol
        p.line(x, y_rut, color=colors["brown"][0])

        # irrigation
        y_irr = self.res["irr"]
        p.circle(x, y_irr, size=5, color="navy", alpha=0.5)

        # annoation: Sowing
        # sowing_arr = Span(x[25], dimension='height', line_color='red',
        # line_dash='dashed', line_width=3)
        # sowing_txt = Label(x=x[25], y=0, text='here your text')

        # p.add_layout(sowing_txt)
        # p.add_layout(sowing_arr)

        return p
