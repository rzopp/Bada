# -*- coding: utf-8 -*-
# original by (c) Daniel González Arribas 
# (c) Flightkeys / Raimund Zopp

import os
import glob
import argparse
import shutil
import math
import psycopg2
from openpyxl import load_workbook
import xml.etree.ElementTree
import numpy as np
import aero
import badadef

class EnvironmentState(object):
    def __init__(self, altitude, disa):     # altitude in m
        self.alt = altitude
        self.disa = disa
        self.delta = aero.delta(altitude)
        self.temp = aero.isaT(altitude)+disa
        self.theta = self.temp / aero.T0

    def get_sound_speed(self):
        return (aero.k * aero.R * self.temp) ** 0.5

    def P(self):
        return self.delta*aero.P0

    def __repr__(self):
        return f"EnvironmentState(altitude={self.alt}, dISA={self.disa}, Temperature={self.temp}, delta={self.delta})"

class AircraftState(object):
    def __init__(self, mach, heading, flight_path_angle, mass, mu=0):
        self.Mach = mach
        self.heading = (heading)
        self.gamma = (flight_path_angle)
        self.m = (mass)
        self.bank = (mu)

    def __repr__(self):
        return f"AircraftState(mach={self.Mach}, heading={self.heading}, flight_path_angle={self.gamma}, mass={self.m}, mu={self.bank})"

    def TAS(self, environment_state):
        return self.Mach * environment_state.get_sound_speed()

    def q(self, environment_state):
        es = environment_state
        return aero.C1*environment_state.P()*self.Mach**2

    def CAS(self, environment_state):
        return (aero.C3 * ((environment_state.delta * ((1 + aero.C5 * self.Mach ** 2) ** aero.C4 - 1) + 1) ** aero.C2 - 1)) ** 0.5

    def CAS_kt(self, environment_state):
        return self.CAS(environment_state) / aero.kt2ms

    def TAS_kt(self, environment_state):
        return self.TAS(environment_state) / aero.kt2ms

    def copy(self):
        return type(self)(self.Mach, self.heading, self.gamma, self.m, self.bank)

class VerticalControls(object):
    def __init__(self, CT: float, CL: float):
        self.CT = CT
        self.CL = CL

class AircraftPerformanceModel(object):
    def __init__(self):
        pass
    
    def v_dot(self, env: EnvironmentState, acs: AircraftState, vc: VerticalControls) -> float:
        T = self.T(env, vc.CT)
        D = self.D(env, acs, vc.CL)
        return (T-D)/acs.m - aero.g0*np.sin(acs.gamma)
    
    def m_dot(self, env: EnvironmentState, acs: AircraftState, vc: VerticalControls):
        CF = self.CF(env, acs, vc.CT)
        return -self.fc(env, CF)

class BADA4_CR(AircraftPerformanceModel):
    """
    Implementation of the BADA4 APM for jet and turboprop aircraft
    """
    def __init__(self, label):
        full_path = os.path.join(badadef.system_root, badadef.bada4_db, label, label + '.xml')
        xml_tree = xml.etree.ElementTree.parse(full_path)
        xml_root = xml_tree.getroot()
        self.ac_label = label
        self.model_label = xml_root.find('model').text
        self.ac_type = xml_root.find('type').text
        ###
        # Aerodynamic Forces and Configurations Model
        ###
        self.AFCM = AFCM = xml_root.find('AFCM')
        self.S = float(AFCM.find('S').text)
        self.configs = AFCM.findall('Configuration')
        ###
        # CR
        ###
        self.clean_cfg = config0 = [cfg for cfg in self.configs if cfg.attrib['HLid'] == '0'][0]
        self.CD_scalar = float(config0.find('LGUP').find('DPM_clean').find('scalar').text)
        self.V_FE_0 = float(config0.find('vfe').text)
        self.CD_clean = CD_clean = config0.find('LGUP').find('DPM_clean').find('CD_clean')
        self.d = [float(d.text) for d in CD_clean.findall('d')]
        LGUP = config0.find('LGUP')
        if self.ac_type == 'JET':
            BLM__ = LGUP.find('BLM_clean')
            if BLM__:
                self.M_min = float(BLM__.find('Mmin').text)
                self.M_max = float(config0.find('LGUP').find('BLM_clean').find('Mmax').text)
                self.CL_clean = CL_clean = config0.find('LGUP').find('BLM_clean').find('CL_clean')
                self.b_f = [float(bf.text) for bf in CL_clean.findall('bf')]
            else:
                self.M_min = 0.2
                self.M_max = float(config0.find('LGUP').find('DPM_clean').find('M_max').text)
                self.CL_max_clean = float(config0.find('LGUP').find('BLM').find('CL_max').text)
                self.b_f = None
        else:   # TURBOPROP
            self.M_min = 0.15
            self.M_max = float(config0.find('LGUP').find('DPM_clean').find('M_max').text)
            self.CL_max_clean = float(config0.find('LGUP').find('BLM').find('CL_max').text)
        ### takeoff CL/CD
        self.to_cfg = config1 = [cfg for cfg in self.configs if cfg.attrib['HLid'] == '1'][0]
        self.CD_to = CD_to = config1.find('LGUP').find('DPM_nonclean').find('CD_nonclean')
        self.d_to = [float(d.text) for d in CD_to.findall('d')]
        LGUP = config1.find('LGUP')
        BLM__ = LGUP.find('BLM')
        self.CL_max_to = float(BLM__.find('CL_max').text)
        ### approach CL/CD
        self.app_cfg = config3 = self.configs[-1]
        self.CD_app = CD_app = config3.find('LGDN').find('DPM_nonclean').find('CD_nonclean')
        self.d_app = [float(d.text) for d in CD_app.findall('d')]
        LGDN = config3.find('LGDN')
        BLM__ = LGDN.find('BLM')
        self.CL_max_app = float(BLM__.find('CL_max').text)
        ###
        ###
        # Propulsive Forces Model
        ###
        self.PFM = PFM = xml_root.find('PFM')
        self.m_ref = float(PFM.find('MREF').text)
        self.LHV = float(PFM.find('LHV').text)
        self.rho_fuel = float(PFM.find('rho').text)
        self.n_eng = float(PFM.find('n_eng').text)
        ##
        if self.ac_type == 'JET':
            # CT
            self.CT_ = PFM.find('TFM').find('CT')
            self.a = [float(a.text) for a in self.CT_.findall('a')]

            ##
            # CF
            CF_ = PFM.find('TFM').find('CF')
            self.f = [float(f.text) for f in CF_.findall('f')]

            ##
            # MCMB
            mclb = PFM.find('TFM').find('MCMB')
            # Flat rating
            flat_rating = mclb.find('flat_rating')
            self.b_MCMB = [float(b.text) for b in flat_rating.findall('b')]

            # Temperature rating
            temp_rating = mclb.find('temp_rating')
            self.c_MCMB = [float(c.text) for c in temp_rating.findall('c')]

            # Kink point
            self.kink_MCMB = float(mclb.find('kink').text)

            ##
            # MCRZ
            mcrz = PFM.find('TFM').find('MCRZ')
            if not mcrz:    # if no max cruise data define, use max climb
                mcrz = mclb
            # Flat rating
            flat_rating = mcrz.find('flat_rating')
            self.b_MCRZ = [float(b.text) for b in flat_rating.findall('b')]

            # Temperature rating
            temp_rating = mcrz.find('temp_rating')
            self.c_MCRZ = [float(c.text) for c in temp_rating.findall('c')]

            # Kink point
            self.kink_MCRZ = float(mcrz.find('kink').text)

            ##
            # LIDL
            CT = PFM.find('TFM').find('LIDL').find('CT')
            self.ti = [float(ti.text) for ti in CT.findall('ti')]

            CF = PFM.find('TFM').find('LIDL').find('CF')
            self.fi = [float(fi.text) for fi in CF.findall('fi')]
        else: # TURBOPROP
            # CP
            self.CP_ = PFM.find('TPM').find('CP')
            self.a = [float(a.text) for a in self.CP_.findall('a')]

            ##
            # CF
            CF_ = PFM.find('TPM').find('CF')
            self.f = [float(f.text) for f in CF_.findall('f')]

            ##
            # MCRZ
            rating = PFM.find('TPM').find('MCRZ').find('rating')
            self.p_MCRZ = [float(p.text) for p in rating.findall('p')]

            ##
            # MCMB
            rating = PFM.find('TPM').find('MCMB').find('rating')
            self.p_MCMB = [float(p.text) for p in rating.findall('p')]

            ##
            # LIDL
            CT = PFM.find('TPM').find('LIDL').find('CT')
            self.ti = [float(ti.text) for ti in CT.findall('ti')]

            CF = PFM.find('TPM').find('LIDL').find('CF')
            self.fi = [float(fi.text) for fi in CF.findall('fi')]

        ###
        # Aircraft Limitations Model
        ###
        ALM = xml_root.find('ALM')

        ##
        # GLM
        self.h_MO_ft = float(ALM.find('GLM').find('hmo').text)
        if self.ac_type == 'JET':
            self.mfa = float(ALM.find('GLM').find('mfa').text)
        else:
            self.mfa = 5000
        ##
        # KLM
        self.M_MO = float(ALM.find('KLM').find('mmo').text)
        self.V_MO_kn = float(ALM.find('KLM').find('vmo').text)
        # mlo = float(ALM.find('KLM').find('mlo').text)
        self.vloe = float(ALM.find('KLM').find('vloe').text)
        self.vlor = float(ALM.find('KLM').find('vlor').text)
        ##
        # DLM
        self.MTW = float(ALM.find('DLM').find('MTW').text)
        self.MTOW = float(ALM.find('DLM').find('MTOW').text)
        self.MLW = float(ALM.find('DLM').find('MLW').text)
        self.MZFW = float(ALM.find('DLM').find('MZFW').text)
        self.OEW = float(ALM.find('DLM').find('OEW').text)
        self.MPL = float(ALM.find('DLM').find('MPL').text)
        self.MFL = float(ALM.find('DLM').find('MFL').text)
        self.n1 = float(ALM.find('DLM').find('n1').text)
        self.n3 = float(ALM.find('DLM').find('n3').text)
        self.nf1 = float(ALM.find('DLM').find('nf1').text)
        self.nf3 = float(ALM.find('DLM').find('nf3').text)

        # speed schedule
        scheds = xml_root.find('ARPM').find('SpeedScheduleList').findall('SpeedSchedule')[0].findall('SpeedPhase')
        self.climbcas = float(scheds[0].find('CAS2').text)
        self.climbMach = float(scheds[0].find('M').text)
        self.desccas = float(scheds[2].find('CAS2').text)
        self.descMach = float(scheds[2].find('M').text)

    def get_fuel_burn_at_cruise_conditions(self,
                                           environment_state: EnvironmentState,
                                           aircraft_state: AircraftState) -> float:
        D = self.D(environment_state, aircraft_state)
        CT = self.CT(environment_state, D)
        CF = self.CF(environment_state, aircraft_state, CT)
        return self.fc(environment_state, CF)

    def get_cruise_CL(self, environment_state, aircraft_state):
        env = environment_state
        acs = aircraft_state
        q = acs.q(env)
        return acs.m * aero.g0 / (self.S * q * np.cos(acs.bank) * np.cos(acs.gamma))

    def CL(self, environment_state, aircraft_state):
        return self.get_cruise_CL(environment_state, aircraft_state)

    def CL_max(self, environment_state, aircraft_state):
        M = max(0.2, aircraft_state.Mach)
        if self.ac_type == 'JET' and self.b_f:
            CL_max =  self.b_f[0] + self.b_f[1] * M + self.b_f[2] * M ** 2 + self.b_f[3] * M ** 3 + self.b_f[4] * M ** 4
        else:
            CL_max = self.CL_max_clean
        if CL_max < 0:
            CL_max = 1
        return CL_max

    def get_bmargin(self, environment_state, aircraft_state):
        return self.CL_max(environment_state, aircraft_state) / self.CL(environment_state, aircraft_state)

    def CD_from_CL_M(self, CL, M):
        d = [self.CD_scalar] + self.d  # align 0-index with 1-index
        # C0 = d[1] + d[2]/(1-softcap(M)**2)**0.5 + d[3]/(1-M**2) + d[4]/(1-M**2)**(3/2) + d[5]/(1-M**2)**2
        C0 = d[1] + d[2] / (1 - M ** 2) ** 0.5 + d[3] / (1 - M ** 2) + d[4] / (1 - M ** 2) ** (3 / 2) + d[5] / (
                1 - M ** 2) ** 2
        C2 = d[6] + d[7] / (1 - M ** 2) ** (3 / 2) + d[8] / (1 - M ** 2) ** 3 + d[9] / (1 - M ** 2) ** (9 / 2) + d[
            10] / (1 - M ** 2) ** 6
        C6 = d[11] + d[12] / (1 - M ** 2) ** (7) + d[13] / (1 - M ** 2) ** (15 / 2) + d[14] / (1 - M ** 2) ** (8) + d[
            15] / (1 - M ** 2) ** (17 / 2)
        return d[0] * (C0 + C2 * CL ** 2 + C6 * CL ** 6)

    def CD(self, environment_state, aircraft_state, CL=None):
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        M = aircraft_state.Mach
        d = [self.CD_scalar] + self.d  # align 0-index with 1-index
        C0 = d[1] + d[2] / (1 - M ** 2) ** 0.5 + d[3] / (1 - M ** 2) + d[4] / (1 - M ** 2) ** (3 / 2) + d[5] / (
                1 - M ** 2) ** 2
        C2 = d[6] + d[7] / (1 - M ** 2) ** (3 / 2) + d[8] / (1 - M ** 2) ** 3 + d[9] / (1 - M ** 2) ** (9 / 2) + d[
            10] / (1 - M ** 2) ** 6
        C6 = d[11] + d[12] / (1 - M ** 2) ** 7 + d[13] / (1 - M ** 2) ** (15 / 2) + d[14] / (1 - M ** 2) ** 8 + d[
            15] / (1 - M ** 2) ** (17 / 2)
        return d[0] * (C0 + C2 * CL ** 2 + C6 * CL ** 6)

    def D(self, environment_state, aircraft_state, CL=None):
        # CD = self.CD(environment_state, aircraft_state, CL)
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        CD = self.CD_from_CL_M(CL, aircraft_state.Mach)
        return aircraft_state.q(environment_state) * self.S * CD

    def D_to(self, environment_state, aircraft_state, CL=None):
        # CD = self.CD(environment_state, aircraft_state, CL)
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        CD = self.d_to[0] + self.d_to[1]*CL + self.d_to[2]*CL** 2
        return aircraft_state.q(environment_state) * self.S * CD

    def D_ldg(self, environment_state, aircraft_state, CL=None):
        # CD = self.CD(environment_state, aircraft_state, CL)
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        CD = self.d_app[0] + self.d_app[1]*CL + self.d_app[2]*CL** 2
        return aircraft_state.q(environment_state) * self.S * CD

    def Vls(self, environment_state, aircraft_state):   # lowest selectable CAS in clean config
        mass = aircraft_state.m
        alt = environment_state.alt
        machls = math.sqrt(mass*aero.g0/(aero.C1*aero.P0*self.S*self.CL_max(environment_state, aircraft_state)/1.3))     # lowest selectable speed in mach
        return aero.mach2cas(machls, alt)

    def V2(self, environment_state, aircraft_state):    # takeoff speed (CAS)
        mass = aircraft_state.m
        alt = environment_state.alt
        machs = math.sqrt(mass*aero.g0/(aero.C1*aero.P0*self.S*self.CL_max_to))     # stall mach in takeoff config
        return aero.mach2cas(machs, alt)*1.2

    def Vapp(self, environment_state, aircraft_state):  # approach speed (CAS)
        mass = aircraft_state.m
        alt = environment_state.alt
        machs = math.sqrt(mass*aero.g0/(aero.C1*aero.P0*self.S*self.CL_max_app))     # stall mach in approach config
        return aero.mach2cas(machs, alt)*1.3

    def CT(self, environment_state, thrust):
        return thrust / environment_state.delta / aero.g0 / self.m_ref

    def T(self, environment_state, CT):
        return CT * environment_state.delta * aero.g0 * self.m_ref

    def CT_min(self, environment_state, aircraft_state):
        delta = environment_state.delta
        theta = environment_state.theta
        M = aircraft_state.Mach
        ti = self.ti
        if self.ac_type == 'JET':
            return (ti[0]*delta**-1 + ti[1] + ti[2] * delta + ti[3] * delta ** 2 +
                    M * (ti[4]*delta**-1 + ti[5] + ti[6] * delta + ti[7] * delta ** 2) +
                    M ** 2 * (ti[8]*delta**-1 + ti[9] + ti[10] * delta + ti[11] * delta ** 2))
        else:   # TURBOPROP
            return (ti[0]*delta**-1 + ti[1] + ti[2] * delta + ti[3] * delta ** 2 +
                    M * (ti[4]*delta**-1 + ti[5] + ti[6] * delta + ti[7] * delta ** 2) +
                    M ** 2 * (ti[8]*delta**-1 + ti[9] + ti[10] * delta + ti[11] * delta ** 2) +
                    ti[12]*theta**0.5 + ti[13]*theta + ti[14]*delta**-0.5 + ti[15]*theta**2 +
                    (ti[16]*delta**-1 + ti[17]*delta + ti[18]*delta**2 + ti[19]*M + ti[20]*M**2)*theta**-0.5 +
                    ti[21]*M**-1 + ti[22]*M**-1*delta + ti[23]*M**3 +
                    (ti[24]*M + ti[25]*M**2 + ti[26] + ti[27]*M*delta**-1)*theta**-1 +
                    ti[28]*M*delta**-1*theta**-2 + ti[29]*M**2*delta**-1*theta**-2 + ti[30]*M**2*delta**-1*theta**-0.5+ti[31]*delta*theta**-1)

    def CT_max_MCRZ(self, environment_state, aircraft_state):
        M = aircraft_state.Mach
        delta = environment_state.delta
        theta_T = environment_state.theta * (1 + (M ** 2 * aero.C5))  # total temperature ratio
        if self.ac_type == 'JET':
            delta_T_flat_MCRZ = 0
            if environment_state.disa <= self.kink_MCRZ:    # flat-rated area
                coeffs = self.b_MCRZ
                for i in range(6):
                    delta_T_flat_MCRZ += delta ** i * sum(coeffs[i * 6 + j] * M ** j for j in range(6))
            else:   # temperature-rated area
                coeffs = self.c_MCRZ
                for i in range(5):
                    delta_T_flat_MCRZ += theta_T ** i * sum(coeffs[i * 5 + j] * M ** j for j in range(5))
                for i in range(1, 5):
                    delta_T_flat_MCRZ += delta ** i * sum(coeffs[(i+4) * 5 + j] * M ** j for j in range(5))
            CT_max = 0
            for i in range(6):
                CT_max += delta_T_flat_MCRZ ** i * sum(self.a[i * 6 + j] * M ** j for j in range(6))
        else:  # TURBOPROP
            coeffs = self.p_MCRZ
            delta_T = 0
            for i in range(6):
                delta_T += theta_T ** i * sum(coeffs[i * 6 + j] * M ** j for j in range(6))
            coeffs = self.a
            CP_max = 0
            for i in range(6):
                CP_max += delta_T ** i * sum(coeffs[i * 6 + j] * M ** j for j in range(6))

            CT_max = CP_max/M

        return CT_max
    
    def CT_max_MCMB(self, environment_state, aircraft_state):
        M = aircraft_state.Mach
        delta = environment_state.delta
        theta_T = environment_state.theta * (1 + (M ** 2 * aero.C5))  # total temperature ratio
        if self.ac_type == 'JET':
            delta_T_flat_MCMB = 0
            if environment_state.disa <= self.kink_MCMB:    # flat-rated area
                coeffs = self.b_MCMB
                for i in range(6):
                    delta_T_flat_MCMB += delta ** i * sum(coeffs[i * 6 + j] * M ** j for j in range(6))
            else:   # temperature-rated area
                coeffs = self.c_MCMB
                for i in range(5):
                    delta_T_flat_MCMB += theta_T ** i * sum(coeffs[i * 5 + j] * M ** j for j in range(5))
                for i in range(1, 5):
                    delta_T_flat_MCMB += delta ** i * sum(coeffs[(i+4) * 5 + j] * M ** j for j in range(5))
            CT_max = 0
            for i in range(6):
                CT_max += delta_T_flat_MCMB ** i * sum(self.a[i * 6 + j] * M ** j for j in range(6))
        else:   #   TURBOPROP
            coeffs = self.p_MCMB
            delta_T = 0
            for i in range(6):
                delta_T += theta_T ** i * sum(coeffs[i * 6 + j] * M ** j for j in range(6))
            coeffs = self.a
            CP_max = 0
            for i in range(6):
                CP_max += delta_T ** i * sum(coeffs[i * 6 + j] * M ** j for j in range(6))
            CT_max = CP_max/M

        return CT_max

    def CF_from_CT_M(self, CT, M):
        CF = self._CF(M, CT)
        return CF

    def fc_from_thrust_P_M_T(self, thrust, P, M, T):
        env = EnvironmentState(P, T)    # ???
        CT = self.CT(env, thrust)
        CF = self._CF(M, CT)
        fc = self.fc(env, CF)
        return fc

    def _CF(self, M, CT):
        CF = 0
        if self.ac_type == 'JET':
            for i in range(5):
                CF += M ** i * sum(self.f[i * 5 + j] * CT ** j for j in range(5))
        else:   # TURBOPROP
            CP = CT*M
            for i in range(5):
                CF += M ** i * sum(self.f[i * 5 + j] * CP ** j for j in range(5))
        return CF

    def CF(self, environment_state, aircraft_state, CT):
        M = aircraft_state.Mach
        return self._CF(M, CT)

    def fc(self, environment_state, CF):
        delta = environment_state.delta
        theta = environment_state.theta
        return delta * theta ** 0.5 * self.m_ref * aero.g0 * aero.a0 * CF / self.LHV

    def CF_idle(self, environment_state, aircraft_state):
        M = aircraft_state.Mach
        delta = environment_state.delta
        theta = environment_state.theta
        CFi = 0
        for i in range(3):
            CFi += M ** i * sum(self.fi[i * 3 + j] * delta ** j for j in range(3))
        if self.ac_type != 'JET':   # TURBOPROP
            CFi += self.fi[9]*theta + self.fi[10]*theta**2 + self.fi[11]*M*theta + self.fi[12]*M*delta*theta**0.5 + self.fi[13]*M*delta*theta
        CFi *= delta**-1 * theta **-0.5
        return CFi

    def takeoff(self, disa, mass):
        alt_clean = 1500*aero.ft2m
        v2acc   = 5     # V2 acceleration to cleanup altitude in m/s
        ES = EnvironmentState(0, disa)
        AS = AircraftState(0, 0, 0, mass, 0)
        mach_v2 = aero.cas2mach(self.V2(ES, AS), ES.alt)
        AS = AircraftState(mach_v2, 0, 0, mass, 0)
        tas_v2 = AS.TAS(ES)
        D_to = self.D_to(ES, AS)
        CT = self.CT_max_MCMB(ES, AS)
        T_to = self.T(ES, CT)
        CF = self.CF(ES, AS, CT)
        fc = self.fc(ES, CF)
        # takeoff acceleration on the runway to V2
        t1 = tas_v2*mass/(T_to - D_to*0.5)
        d1 = tas_v2*t1/2
        f1 = fc*t1
        fnexc = T_to - D_to
        roc = aero.getrocacc(alt_clean/2, alt_clean, disa, mass, tas_v2, 5, fnexc)
        t2 = alt_clean/roc
        d2 = (tas_v2+v2acc/2)*t2
        f2 = fc*t2
        return alt_clean, mach_v2, tas_v2, t1+t2, d1+d2, f1+f2

    def approach(self, disa, mass):
        alt_app = 1500*aero.ft2m
        vappinc   = 5     # approach speed increment in m/s
        ES = EnvironmentState(0, disa)
        AS = AircraftState(0, 0, 0, mass, 0)
        mach_app = aero.cas2mach(self.Vapp(ES, AS), ES.alt)
        AS = AircraftState(mach_app, 0, math.radians(-3), mass, 0)
        tas_app = AS.TAS(ES)
        D_app = self.D_ldg(ES, AS)
        CT = self.CT_min(ES, AS)
        T_app = self.T(ES, CT)
        CF = self.CF(ES, AS, CT)
        fc = self.fc(ES, CF)
        # landing deceleration on the runway from Vapp
        t1 = -tas_app*mass/(T_app - D_app*0.5 - mass*aero.g0*0.3)    # deceleration at mu = 0.3
        d1 = tas_app*t1/2
        f1 = fc*t1
        roc = -tas_app*math.sin(math.radians(3))   # rate for 3° approach
        fnexc = aero.getfnexc(alt_app, disa, mass, mach_app, roc/aero.fpm2mps, 'CAS')
        CT2 = self.CT(ES, fnexc + D_app)
        CF2 = self.CF(ES, AS, CT2)
        fc2 = self.fc(ES, CF2)
        t2 = -alt_app/roc
        d2 = (tas_app+vappinc/2)*t2
        f2 = fc2*t2
        return alt_app, mach_app, tas_app, t1+t2, d1+d2, f1+f2

    def minCAS(self, mass, alt, disa):
        ES = EnvironmentState(alt, disa)
        AS = AircraftState(0, 0, 0, mass)
        return self.Vls(ES, AS)

    def clbpoint(self, mode, disa, mass, alt1, alt2, tas1, tas2, reqrate):   # altitude in m
        tasm = (tas1 + tas2) / 2
        altm = (alt1 + alt2) / 2
        dtas = tas2 - tas1
        dalt = alt2 - alt1
        machm = (aero.tas2mach(tas1, alt1, disa)+aero.tas2mach(tas2, alt2, disa))/2
        casm = (aero.tas2cas(tas1, alt1, disa) + aero.tas2cas(tas2, alt2, disa))/2
        ES = EnvironmentState(altm, disa)
        AS = AircraftState(machm, 0, 0, mass, 0)
        bmargin = self.get_bmargin(ES, AS)
        if alt2 < 10000*aero.ft2m:  # ignore buffet margin below 10000ft
            bmargin = max(2, bmargin)
        CL = self.CL(ES, AS)
        CD = self.CD(ES, AS, CL)
        CT = self.CT_max_MCMB(ES, AS)   # max climb thrust
        T = self.T(ES, CT)  # thrust
        D = self.D(ES, AS)       # drag
        reqfn = aero.getreqfn(altm, dalt, disa, mass, tasm, dtas, D, reqrate)   # thrust required to achieve requested rate
        if reqfn > T:   # can't achieve required rate
            fnexc = T - D            # excess thrust
            roc = aero.getrocacc(altm, dalt, disa, mass, tasm, dtas, fnexc)
        else:
            T = reqfn
            CT = self.CT(ES,T)
            roc = reqrate
        CF = self.CF(ES, AS, CT)
        wfe = self.fc(ES, CF)
        return machm, tasm, casm, wfe, roc, CL, CD, T, D, bmargin

    def dscpoint(self, mode, disa, mass, alt1, alt2, tas1, tas2, reqrate):   # altitude in m
        tasm = (tas1 + tas2) / 2
        altm = (alt1 + alt2) / 2
        dtas = tas2 - tas1
        dalt = alt2 - alt1
        machm = (aero.tas2mach(tas1, alt1, disa)+aero.tas2mach(tas2, alt2, disa))/2
        casm = (aero.tas2cas(tas1, alt1, disa) + aero.tas2cas(tas2, alt2, disa))/2
        ES = EnvironmentState(altm, disa)
        AS = AircraftState(machm, 0, 0, mass, 0)
        bmargin = self.get_bmargin(ES, AS)
        if alt2 < 10000*aero.ft2m:  # ignore buffet margin below 10000ft
            bmargin = max(2, bmargin)
        CL = self.CL(ES, AS)
        CD = self.CD(ES, AS, CL)
        CT = self.CT_min(ES, AS)        # idle thrust
        T = self.T(ES, CT)  # thrust
        D = self.D(ES, AS)       # drag
        reqfn = aero.getreqfn(altm, dalt, disa, mass, tasm, dtas, D, reqrate)   # thrust required to achieve requested rate
        if reqfn < T:   # can't achieve required rate with idle thrust
            fnexc = T - D            # excess thrust
            roc = aero.getrocacc(altm, dalt, disa, mass, tasm, dtas, fnexc)
            CF = self.CF_idle(ES, AS)
        else:
            T = reqfn
            CT = self.CT(ES,T)
            roc = reqrate
            CF = self.CF(ES, AS, CT)
        wfe = self.fc(ES, CF)
        return machm, tasm, casm, wfe, roc, CL, CD, T, D, bmargin

    def cruisepoint(self, mass, alt, disa, mach):
        altm = alt*aero.ft2m
        ES = EnvironmentState(altm, disa)
        AS = AircraftState(mach, 0, 0, mass, 0)
        tas = AS.TAS(ES)
        cas = AS.CAS(ES)
        wfe = 0
        bmargin = self.get_bmargin(ES, AS)
        if  cas <= self.V_MO_kn*aero.kt2ms and bmargin >= 1.0:      # the rest only makes sense if we have at least 1g margin
            CTmax = self.CT_max_MCRZ(ES, AS)
            CT = self.CT(ES, self.D(ES, AS))
            fnexc = self.T(ES, CTmax-CT)           # excess thrust
            roc = aero.getroc(alt, disa, mass, mach, fnexc, 'Mach')
            if roc >= 100*aero.fpm2mps:
                CF = self.CF(ES, AS, CT)
                wfe = self.fc(ES, CF)
        return tas, cas, wfe, bmargin

def bada4climb(sched_cas, sched_mach, reqrate):
    out = open(outfilename + '.CSV', 'w')
    out.write("BADA_CLB\n")
    out.write("Climb data created for {} from BADA4 file {}.xml\n".format(subtype, general_db))
    out.write('DISA,INITMASS,MASS,ALT,TIME,FUEL,DIST,CAS,MACH,TAS,WFE,ROC,GRAD,ALPHA,CL,CD,FN,DRAG,THR,FNEXC\n')
    out.write('DEGC,KG,KG,FT,MIN,KG,NM,KTS,,KTS,KG/H,FT/MIN,,DEG,,,DAN,DAN,,DAN\n')
    for disa in np.arange(0, 45, 5):
        for init_mass in np.arange (base_mass, apm.MTOW+mass_step, mass_step):
            mass = init_mass  # starting at takeoff mass
            out.write('{},{:.0f},{:.3f},{:.1f},{:.6f},{:.6f},{:.6f},{:.3f},{:.8f},{:.3f},{:.5f},{:.5f},{:.8f},{:.3f},{:.9f},{:.9f},{:.4f},{:.4f},{:.4f},{:.4f}\n'.format(
                    disa, init_mass, mass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            prev_alt, prev_mach, prev_tas, time, dist, fuel = apm.takeoff(disa, init_mass)
            while prev_alt < apm.h_MO_ft*aero.ft2m:
                next_alt = (int(round(prev_alt/aero.ft2m)/1000)+1)*1000*aero.ft2m
                if next_alt <= 10000*aero.ft2m:
                    req_cas = min(max(250*aero.kt2ms, apm.minCAS(mass, next_alt, disa)), sched_cas*aero.kt2ms)
                else:
                    req_cas = sched_cas*aero.kt2ms
                if prev_mach < sched_mach:    # we are still in the cas range
                    next_mach = min(sched_mach, aero.cas2mach(req_cas, next_alt))
                next_tas = aero.mach2tas(next_mach, next_alt, disa)  # tas at upper altitude
                mach, tas, cas, wfe, roc, cl, cd, fn, drag, bmargin = apm.clbpoint('CLB', disa, mass, prev_alt, next_alt, prev_tas, next_tas, reqrate*aero.fpm2mps)
                t = (next_alt - prev_alt)/roc    # time in seconds
                d = tas*t       # distance in meters
                f = wfe*t       # fuel in kg
                time = time + t
                dist = dist + d
                fuel = fuel + f
                mass = mass - f
                if t > 0 and roc >= 300*aero.fpm2mps and bmargin >= 1.3:
                    out.write(
                        '{},{:.0f},{:.3f},{:.1f},{:.6f},{:.6f},{:.6f},{:.3f},{:.8f},{:.3f},{:.5f},{:.5f},{:.8f},{:.3f},{:.9f},{:.9f},{:.4f},{:.4f},{:.4f},{:.4f}\n'.format(
                            disa, init_mass, mass, next_alt/aero.ft2m, time/60, fuel, dist/aero.nm2m, cas/aero.kt2ms, mach, tas/aero.kt2ms, wfe*3600, roc/aero.fpm2mps, 0, 0, cl, cd, fn/10, drag/10, 0, (fn-drag)/10))
                else:
                    break
                prev_alt = next_alt
                prev_tas = next_tas
                prev_mach = next_mach

    out.close()
    return

def bada4descent(sched_cas, sched_mach, reqrate):
    out = open(outfilename + '.CSV', 'w')
    out.write("BADA_DSC\n")
    out.write("Descent data created for {} from BADA4 file {}.xml\n".format(subtype, general_db))
    out.write('DISA,INITMASS,MASS,ALT,TIME,FUEL,DIST,CAS,MACH,TAS,WFE,ROC,GRAD,ALPHA,CL,CD,FN,DRAG,THR,FNEXC\n')
    out.write('DEGC,KG,KG,FT,MIN,KG,NM,KTS,,KTS,KG/H,FT/MIN,,DEG,,,DAN,DAN,,DAN\n')
    for disa in np.arange(-20, 45, 5):
        for init_mass in np.arange (base_mass, apm.MTOW+mass_step, mass_step):
            mass = init_mass  # starting at landing mass
            out.write(
                '{},{:.0f},{:.3f},{:.1f},{:.6f},{:.6f},{:.6f},{:.3f},{:.8f},{:.3f},{:.5f},{:.5f},{:.8f},{:.3f},{:.9f},{:.9f},{:.4f},{:.4f},{:.4f},{:.4f}\n'.format(
                    disa, init_mass, mass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            prev_alt, prev_mach, prev_tas, time, dist, fuel = apm.approach(disa, init_mass)
            while prev_alt < apm.h_MO_ft * aero.ft2m:
                next_alt = (int(round(prev_alt / aero.ft2m) / 1000) + 1) * 1000 * aero.ft2m
                if next_alt <= 10000 * aero.ft2m:
                    req_cas = min(max(250 * aero.kt2ms, apm.minCAS(mass, next_alt, disa)), sched_cas * aero.kt2ms)
                else:
                    req_cas = sched_cas * aero.kt2ms
                if prev_mach < sched_mach:  # we are still in the cas range
                    next_mach = min(sched_mach, aero.cas2mach(req_cas, next_alt))
                next_tas = aero.mach2tas(next_mach, next_alt, disa)  # tas at upper altitude
                mach, tas, cas, wfe, roc, cl, cd, fn, drag, bmargin = apm.dscpoint('DSC', disa, mass, prev_alt, next_alt, prev_tas, next_tas, reqrate*aero.fpm2mps)
                t = -(next_alt - prev_alt) / roc  # time in seconds
                d = tas * t  # distance in meters
                f = wfe * t  # fuel in kg
                time = time + t
                dist = dist + d
                fuel = fuel + f
                mass = mass - f
                if t > 0 and bmargin >= 1.3:
                    out.write(
                        '{},{:.0f},{:.3f},{:.1f},{:.6f},{:.6f},{:.6f},{:.3f},{:.8f},{:.3f},{:.5f},{:.5f},{:.8f},{:.3f},{:.9f},{:.9f},{:.4f},{:.4f},{:.4f},{:.4f}\n'.format(
                            disa, init_mass, mass, next_alt / aero.ft2m, time / 60, fuel, dist / aero.nm2m,
                                                   cas / aero.kt2ms, mach, tas / aero.kt2ms, wfe * 3600,
                                                   roc / aero.fpm2mps, 0, 0, cl, cd, fn / 10, drag / 10, 0,
                                                   (fn - drag) / 10))
                else:
                    break
                prev_alt = next_alt
                prev_tas = next_tas
                prev_mach = next_mach


    out.close()
    return

def bada4cruise():
    out = open(outfilename + '.CSV', 'w')
    out.write("BADA_CRZ\n")
    out.write("Cruise data created for {} from BADA4 file {}.xml\n".format(subtype, general_db))
    out.write("mass, alt, disa, mach, tas, cas, wfe, bmargin\n")
    out.write("(kg), (ft), (K), , (kt), (kt), (kg/h)\n")
    for mass in np.arange(base_mass, apm.MTOW + mass_step, mass_step):
        for alt in np.arange (0, apm.h_MO_ft+1000, 1000):
            for disa in np.arange(0, 45, 5):
                for mach in np.arange(apm.M_min, apm.M_MO+0.01, 0.01):
                    tas, cas, wfe, bmargin = apm.cruisepoint(mass, alt, disa, mach)
                    if wfe > 0:
                        out.write("{}, {}, {}, {}, {}, {}, {}, {}\n".format(mass, alt, disa, mach, tas/aero.kt2ms, cas/aero.kt2ms, wfe*3600.0, bmargin))
    out.close()
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='database connection parameters')
    parser.add_argument('subtype')
    parser.add_argument('version_nr')

    args = parser.parse_args()

    conn = psycopg2.connect(host='fkylabdb.flightkeys.com', database='fk_acperf', user='fkservice', password='ieDie4ieMi3aehaiZi5l')

    cur = conn.cursor()

    sqlcommand = '''
    select v.manufacturer, r.general_db,
    r.base_mass, r.mass_step, r.base_mach, r.mach_step, 
    r.ref_lhv/2326.009143, r.ref_cg, r.min_rate_avail_clb, r.min_rate_avail_crz 
    from acperf.aircraft_version v, acperf.aircraft_release r 
    where r.aircraft_version_uid = v.uid
    and v.subtype = '{}' and v.version_nr = {}
    and r.rel_nr = 0 
    order by 1 desc 
    limit 1
    '''.format(args.subtype, args.version_nr)

    cur.execute(sqlcommand)

    row = cur.fetchone()

    if row != None:
        manufacturer = row[0]
        general_db = row[1]
        base_mass = row[2]
        mass_step = row[3]
        base_mach = row[4]
        mach_step = row[5]
        ref_lhv = row[6]
        ref_cg = row[7]
        min_rate_avail_clb = row[8] / aero.fpm2mps
        min_rate_avail_crz = row[9] / aero.fpm2mps
    else:
        print('version not found:\n')
        quit()

    lmax = 0  # match length
    climbs = []
    descents = []
    wb = load_workbook(filename='{}perfcreation\\BADA_Params.xlsx'.format(badadef.system_root))
    sheet = wb.active
    for row in sheet.values:
        key = row[0].split('/')
        match = False

        if len(key) > 1:  # first check for exact version match
            if args.subtype == key[0] and args.version_nr == key[1]:
                match = True
                lmax = 99
            else:
                continue
        elif args.subtype.find(key[0]) == 0 and len(key[0]) > lmax:  # approximate key match
            match = True
            lmax = len(key[0])
        if match:
            bestkey = key[0]
            climbs = row[1].split(',')
            descents = row[2].split(',')

    cur.close()
    conn.close()

    target_path = badadef.target_root + '{}\\'.format(manufacturer)
    database_path = badadef.system_root + badadef.bada4_db

    apm = BADA4_CR(general_db)
    subtype = apm.model_label

    badaclimb = '{}M{}'.format(int(apm.climbcas), int(apm.climbMach*100))
    badadesc  = '{}M{}'.format(int(apm.desccas), int(apm.descMach*100))

    found = False
    for climb in climbs:
        if climb == badaclimb:
            found = True
            break
    if not found:
        climbs.append(badaclimb)

    found = False
    for desc in descents:
        if desc == badadesc:
            found = True
            break
    if not found:
        descents.append(badadesc)
    # delete all existing files for this version in work directory:
    fileList = glob.glob("{}_{:0>2}_*.*".format(args.subtype, args.version_nr))
    for filePath in fileList:
        try:
            os.remove(filePath)
        except:
            print("Error while deleting file : ", filePath)

    # create climb files:
    if badadef.do_climb:
        for climb in climbs:
            reqrate = 9E99
            if climb[0] == 'I':   # IAS-only schedule
                cas = float(climb[1:])
                mach = 0.99
            else:
                sched = climb.split('M')
                cas = float(sched[0])
                mach = float(sched[1]) / 100.0
            outfilename = "{}_{:0>2}_CLB{}_000022".format(args.subtype, args.version_nr, climb)
            print(outfilename)
            bada4climb(cas, mach, reqrate)

    if badadef.do_descent:
        for descent in descents:
            if descent[0] == 'I':   # IAS-only schedule
                cas = float(descent[1:])
                mach = 0.99
                reqrate = -1500      # ft/min
            else:
                sched = descent.split('M')
                cas = float(sched[0])
                mach = float(sched[1]) / 100.0
                reqrate = -9E99
            outfilename = "{}_{:0>2}_DSC{}_000020".format(args.subtype, args.version_nr, descent)
            print(outfilename)
            bada4descent(cas, mach, reqrate)

    if badadef.do_cruise_ae:
        outfilename = "{}_{:0>2}_CRZECON_000021".format(args.subtype, args.version_nr)
        print(outfilename)
        bada4cruise()

    # move created sourcefiles to target directory:
    fileList = glob.glob("{}_{:0>2}_*.*".format(args.subtype, args.version_nr))
    for filePath in fileList:
        dest_file = target_path + filePath
        try:
            if os.path.exists(dest_file):
                os.remove(dest_file)
            shutil.move(filePath, target_path)
        except:
            print("Error while moving file {} to {}".format(filePath, target_path))