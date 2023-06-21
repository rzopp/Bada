# -*- coding: utf-8 -*-
# (c) Flightkeys / Raimund Zopp

import os
import glob
import argparse
import shutil
import math
import psycopg2
from openpyxl import load_workbook
import numpy as np
import aero
import badadef

def getvalue(defval, val):
    if val == 0:
        return defval
    else:
        return val

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
        return f"EnvironmentState(altitude={self.altitude}, dISA={self.disa}, Temperature={self.temp}, delta={self.delta})"

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

class AircraftPerformanceModel(object):
    def __init__(self):
        pass
    
class BADA3_CR(AircraftPerformanceModel):
    """
    Implementation of the BADA3 APM for turbofan jet  and turboprop aircraft
    """
    def __init__(self, label):
        # default drag coefficients
        self.CD0 = 0.035
        self.CD2 = 0.045
        self.CD0_to = 0.035
        self.CD2_to = 0.045
        self.CD0_app = 0.039
        self.CD2_app = 0.045
        self.CD0_ldg = 0.079
        self.CD2_ldg = 0.045

        full_path = os.path.join(badadef.system_root, badadef.bada3_db, label.ljust(6,'_') + '.OPF')
        opf = open(full_path, 'r')
        for line in opf:
            if line[9:].find('Actype') == 0:
                section = 1
                i = 0
            if line[9:].find('Mass') == 0:
                section = 2
                i = 0
            if line[9:].find('Flight envelope') == 0:
                section = 3
                i = 0
            if line[9:].find('Aerodynamics') == 0:
                section = 4
                i = 0
            if line[9:].find('Engine Thrust') == 0:
                section = 5
                i = 0
            if line[9:].find('Fuel Consumption') == 0:
                section = 6
                i = 0
            if line[9:].find('Ground') == 0:
                section = 7
                i = 0
            elif line[:2] == 'CD':      # data line
                i = i + 1
                if section == 1:    # Actype
                    data = line[2:].split()
                    if i == 1:
                        self.ac_label = data[0].strip('_')
                        self.ac_type = data[3]
                    elif i == 2:
                        self.model_label = data[0]
                        self.engine_series = data[2]
                elif section == 2:  # mass
                    data = line[2:].split()
                    self.ref_mass = float(data[0])*1000
                    self.min_mass = float(data[1])*1000
                    self.max_mass = float(data[2])*1000
                elif section == 3:  # flight envelope
                    data = line[2:].split()
                    self.V_MO = float(data[0])*aero.kt2ms
                    self.M_MO = float(data[1])
                    self.h_MO_ft = float(data[2])
                elif section == 4:  # aerodynamics
                    if i == 1:
                        data = line[2:].split()
                        self.surf = float(data[1])
                        self.CL_b0 = float(data[2])
                        self.k = float(data[3])
                        self.CM16 = float(data[4])
                    elif i == 2:
                        data = line[20:].split()
                        self.Vs = float(data[0])*aero.kt2ms
                        self.CD0 = getvalue(self.CD0, float(data[1]))
                        self.CD2 = getvalue(self.CD2, float(data[2]))
                    elif i == 4:
                        data = line[20:].split()
                        self.Vs_to = float(data[0]) * aero.kt2ms
                        self.CD0_to = getvalue(self.CD0_to, float(data[1]))
                        self.CD2_to = getvalue(self.CD2_to, float(data[2]))
                    elif i == 5:
                        data = line[20:].split()
                        self.Vs_app = float(data[0]) * aero.kt2ms
                        self.CD0_app = getvalue(self.CD0_app, float(data[1]))
                        self.CD2_app = getvalue(self.CD2_app, float(data[2]))
                    elif i == 6:
                        data = line[20:].split()
                        self.Vs_ldg = float(data[0]) * aero.kt2ms
                        self.CD0_ldg = getvalue(self.CD0_ldg, float(data[1]))
                        self.CD2_ldg = getvalue(self.CD2_ldg, float(data[2]))
                elif section == 5:  # engine thrust
                    data = line[2:].split()
                    if i == 1:
                        self.Ctc1 = float(data[0])
                        self.Ctc2 = float(data[1])
                        self.Ctc3 = float(data[2])
                        self.Ctc4 = float(data[3])
                        self.Ctc5 = float(data[4])
                    if i == 2:
                        self.CTdes_lo = float(data[0])
                        self.CTdes_hi = float(data[1])
                        self.CTdes_alt = float(data[2])*aero.ft2m
                        self.CTdes_app = float(data[3])
                        self.CTdes_ld = float(data[4])
                    if i == 3:
                        self.descCas = float(data[0])
                        self.descMach = float(data[1])
                        self.climbCas = self.descCas
                        self.climbMach = self.descMach
                elif section == 6:  # fuel consumption
                    data = line[2:].split()
                    if i == 1:
                        self.Cf1 = float(data[0])
                        self.Cf2 = float(data[1])
                    if i == 2:
                        self.Cf3 = float(data[0])
                        self.Cf4 = float(data[1])
                    if i == 3:
                        self.Cfcr = float(data[1])
        opf.close()
        full_path = os.path.join(badadef.system_root, badadef.bada3_db, label.ljust(6,'_') + '.APF')
        apf = open(full_path, 'r')
        for line in apf:
            data = line[23:].split()
            if line[:2] == 'CD' and data[0] == 'AV':      # data line
                self.climbCas = float(data[2])
                self.climbMach = float(data[3])/100
                self.descCas = float(data[8])
                self.descMach = float(data[7])/100
                break
        apf.close()



        self.M_min = round(aero.cas2mach(self.Vs*1.5,0),2)
        self.Ctcr = 0.9     # assuming max cruise thrust is 90% of max climb thrust
        self.CL_max_clean = 2*self.ref_mass*aero.g0/(aero.rho0*self.Vs**2*self.surf)
        self.CL_max_to = 2.0
        self.CL_max_app = 2.5
        self.k1 = aero.kt2ms**2
        self.k2 = 1/60/1000
        self.k3 = 1/60
        self.k4 = aero.nm2m/aero.ft2m/60
        if self.ac_type == 'Jet':
            self.t0 = self.Ctc1
            self.t1 = self.Ctc1/self.Ctc2
            self.t2 = 0
            self.t3 = 0
            self.t4 = self.Ctc1*self.Ctc3
            self.f0 = 0
            self.f2 = self.k2*self.Cf1
            self.f3 = self.k2*self.Cf1/self.Cf2
            self.f4 = 0
            self.f5 = 1
        elif self.ac_type == 'Turboprop':
            self.t0 = self.Ctc3
            self.t1 = 0
            self.t2 = self.Ctc1
            self.t3 = self.Ctc1/self.Ctc2
            self.t4 = 0
            self.f0 = 0
            self.f2 = 0
            self.f3 = self.k2*self.Cf1/1000
            self.f4 = self.k2/1000*self.Cf1/self.Cf2
            self.f5 = 1
        self.t5 = self.Ctc5
        self.t6 = 1 + self.Ctc4*self.Ctc5
        self.f1 = 0

    def get_wfe(self,environment_state: EnvironmentState, aircraft_state: AircraftState, T) -> float:
        alt = environment_state.alt/aero.ft2m
        tas = aircraft_state.TAS(environment_state)/aero.kt2ms
        wfe = self.f5*(self.f0-self.f1*alt+(self.f2+self.f3*tas-self.f4*tas**2)*T)
        return wfe

    def get_cruise_wfe(self, environment_state: EnvironmentState, aircraft_state: AircraftState) -> float:
        T = D = self.D(environment_state, aircraft_state)
        alt = environment_state.alt/aero.ft2m
        tas = aircraft_state.TAS(environment_state)/aero.kt2ms
        wfe = self.Cfcr*(self.f0-self.f1*alt+(self.f2+self.f3*tas-self.f4*tas**2)*T)
        return wfe

    def get_cruise_CL(self, environment_state, aircraft_state):
        env = environment_state
        acs = aircraft_state
        q = acs.q(env)
        return acs.m * aero.g0 / (self.surf * q)

    def CL(self, environment_state, aircraft_state):
        return self.get_cruise_CL(environment_state, aircraft_state)

    def CL_max(self, environment_state, aircraft_state):
        M = aircraft_state.Mach
        return self.CL_max_clean    # self.b_f[0] + self.b_f[1] * M + self.b_f[2] * M ** 2 + self.b_f[3] * M ** 3 + self.b_f[4] * M ** 4

    def get_bmargin(self, environment_state, aircraft_state):
        return self.CL_max(environment_state, aircraft_state) / self.CL(environment_state, aircraft_state)

    def CD_from_CL_M(self, CL, M):
        return self.CD0 + self.CD2 * CL**2

    def CD(self, environment_state, aircraft_state, CL=None):
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        return self.CD0 + self.CD2 * CL**2

    def CD_app(self, environment_state, aircraft_state, CL=None):
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        return self.CD0_app + self.CD2_app * CL**2

    def CD_ldg(self, environment_state, aircraft_state, CL=None):
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        return self.CD0_ldg + self.CD2_ldg * CL**2

    def D_to(self, environment_state, aircraft_state):
        CL = self.CL(environment_state, aircraft_state)
        CD = self.CD0_to + self.CD2_to * CL**2
        return aircraft_state.q(environment_state) * self.surf * CD

    def D_app(self, environment_state, aircraft_state):
        CL = self.CL(environment_state, aircraft_state)
        CD = self.CD0_app + self.CD2_app * CL**2
        return aircraft_state.q(environment_state) * self.surf * CD

    def D_ldg(self, environment_state, aircraft_state):
        CL = self.CL(environment_state, aircraft_state)
        CD = self.CD0_ldg + self.CD2_ldg * CL**2
        return aircraft_state.q(environment_state) * self.surf * CD

    def D(self, environment_state, aircraft_state, CL=None):
        # CD = self.CD(environment_state, aircraft_state, CL)
        if CL is None:
            CL = self.CL(environment_state, aircraft_state)
        CD = self.CD_from_CL_M(CL, aircraft_state.Mach)
        return aircraft_state.q(environment_state) * self.surf * CD

    def Vls(self, environment_state, aircraft_state):   # lowest selectable CAS in clean config
        mass = aircraft_state.m
        alt = environment_state.alt
        machls = math.sqrt(mass*aero.g0/(aero.C1*aero.P0*self.surf*self.CL_max(environment_state, aircraft_state)/1.3))     # lowest selectable speed in mach
        return aero.mach2cas(machls, alt)

    def V2(self, environment_state, aircraft_state):    # takeoff speed (CAS)
        mass = aircraft_state.m
        alt = environment_state.alt
        machs = math.sqrt(mass*aero.g0/(aero.C1*aero.P0*self.surf*self.CL_max_to))     # stall mach in takeoff config
        return aero.mach2cas(machs, alt)*1.2

    def Vapp(self, environment_state, aircraft_state):  # approach speed (CAS)
        mass = aircraft_state.m
        alt = environment_state.alt
        machs = math.sqrt(mass*aero.g0/(aero.C1*aero.P0*self.surf*self.CL_max_app))     # stall mach in approach config
        return aero.mach2cas(machs, alt)*1.3

    def T_max_MCMB(self, environment_state, aircraft_state):
        tas = aircraft_state.TAS(environment_state)/aero.kt2ms
        alt = environment_state.alt/aero.ft2m
        disa = environment_state.disa

        Tisa = self.t0 - self.t1*alt + self.t2/tas - self.t3*alt/tas + self.t4*alt**2

        return Tisa*(self.t6-self.t5*disa)

    def T_max_MCRZ(self, environment_state, aircraft_state):
        return self.Ctcr*self.T_max_MCMB(environment_state, aircraft_state)

    def T_idle(self, environment_state, aircraft_state):
        if environment_state.alt >= self.CTdes_alt:
            CTidle = self.CTdes_hi
        else:
            CTidle = self.CTdes_lo
        return CTidle*self.T_max_MCMB(environment_state, aircraft_state)

    def get_wfe_idle(self, environment_state, aircraft_state):
        alt = environment_state.alt/aero.ft2m
        f0 = self.k3*self.Cf3           # local coefficient for idle calculation only
        f1 = self.k3*self.Cf3/self.Cf4  # local coefficient for idle calculation only

        F = f0-f1*alt
        return F

    def takeoff(self, disa, mass):
        alt_clean = 1500*aero.ft2m
        v2acc   = 5     # V2 acceleration to cleanup altitude in m/s
        ES = EnvironmentState(0, disa)
        AS = AircraftState(0, 0, 0, mass, 0)
        mach_v2 = aero.cas2mach(self.V2(ES, AS), ES.alt)
        AS = AircraftState(mach_v2, 0, 0, mass, 0)
        tas_v2 = AS.TAS(ES)
        D_to = self.D_to(ES, AS)
        T_to = self.T_max_MCMB(ES, AS)
        wfe = self.get_wfe(ES, AS, T_to)
        # takeoff acceleration on the runway to V2
        t1 = tas_v2*mass/(T_to - D_to*0.5)
        d1 = tas_v2*t1/2
        f1 = wfe*t1
        fnexc = T_to - D_to
        roc = aero.getrocacc(alt_clean/2, alt_clean, disa, mass, tas_v2, 5, fnexc)
        t2 = alt_clean/roc
        d2 = (tas_v2+v2acc/2)*t2
        f2 = wfe*t2
        return alt_clean, mach_v2, tas_v2, t1+t2, d1+d2, f1+f2

    def approach(self, disa, mass):
        alt_app = 1500*aero.ft2m
        vappinc   = 5     # approach speed increment in m/s
        ES = EnvironmentState(0, disa)
        AS = AircraftState(0, 0, 0, mass, 0)
        mach_app = aero.cas2mach(self.Vapp(ES, AS), ES.alt)
        AS = AircraftState(mach_app, 0, 0, mass, 0)
        tas_app = AS.TAS(ES)
        D_app = self.D_app(ES, AS)
        T_app = self.CTdes_app*self.T_max_MCMB(ES, AS)
        wfe = self.get_wfe_idle(ES, AS)
        # landing deceleration on the runway from Vapp
        t1 = -tas_app*mass/(T_app - D_app*0.5 - mass*aero.g0*0.3)    # deceleration at mu = 0.3
        d1 = tas_app*t1/2
        f1 = wfe*t1
        roc = -tas_app*math.sin(math.radians(3))   # rate for 3° approach
        fnexc = aero.getfnexc(alt_app, disa, mass, mach_app, roc, 'CAS')
        wfe2 = self.get_wfe(ES, AS, fnexc + D_app)
        t2 = -alt_app/roc
        d2 = (tas_app+vappinc/2)*t2
        f2 = wfe2*t2
        return alt_app, mach_app, tas_app, t1+t2, d1+d2, f1+f2

    def minCAS(self, mass, alt, disa):
        ES = EnvironmentState(alt, disa)
        AS = AircraftState(0, 0, 0, mass)
        return self.Vls(ES, AS)

    def clbpoint(self, mode, disa, mass, alt1, alt2, tas1, tas2, reqrate):  # altitude in m
        tasm = (tas1 + tas2) / 2
        altm = (alt1 + alt2) / 2
        dtas = tas2 - tas1
        dalt = alt2 - alt1
        machm = (aero.tas2mach(tas1, alt1, disa)+aero.tas2mach(tas2, alt2, disa))/2
        casm = (aero.tas2cas(tas1, alt1, disa) + aero.tas2cas(tas2, alt2, disa))/2
        ES = EnvironmentState(altm, disa)
        AS = AircraftState(machm, 0, 0, mass, 0)
        if alt2 <= 2000*aero.ft2m:
            bmargin = 1.5   # ignore buffet limit in takeoff climb acceleration phase
        else:
            bmargin = self.get_bmargin(ES, AS)
        CL = self.CL(ES, AS)
        CD = self.CD(ES, AS, CL)
        T = self.T_max_MCMB(ES, AS)  # max climb thrust
        D = self.D(ES, AS)  # drag
        reqfn = aero.getreqfn(altm, dalt, disa, mass, tasm, dtas, D,
                              reqrate)  # thrust required to achieve requested rate
        if reqfn > T:  # can't achieve required rate
            fnexc = T - D  # excess thrust
            roc = aero.getrocacc(altm, dalt, disa, mass, tasm, dtas, fnexc)
        else:
            T = reqfn
            roc = reqrate
        wfe = self.get_wfe(ES, AS, T)
        return machm, tasm, casm, wfe, roc, CL, CD, T, D, bmargin

    def dscpoint(self, mode, disa, mass, alt1, alt2, tas1, tas2, reqrate):  # altitude in m
        tasm = (tas1 + tas2) / 2
        altm = (alt1 + alt2) / 2
        dtas = tas2 - tas1
        dalt = alt2 - alt1
        machm = (aero.tas2mach(tas1, alt1, disa)+aero.tas2mach(tas2, alt2, disa))/2
        casm = (aero.tas2cas(tas1, alt1, disa) + aero.tas2cas(tas2, alt2, disa))/2
        ES = EnvironmentState(altm, disa)
        AS = AircraftState(machm, 0, 0, mass, 0)
        CL = self.CL(ES, AS)
        if alt2 <= 2000*aero.ft2m:
            bmargin = 1.5   # ignore buffet limit in approach deceleration phase
            CD = self.CD_ldg(ES, AS, CL)
            D = self.D_ldg(ES, AS)    # approach drag coefficient
        else:
            bmargin = self.get_bmargin(ES, AS)
            D = self.D(ES, AS)  # drag
            CD = self.CD(ES, AS, CL)
        T = self.T_idle(ES, AS)  # thrust
        reqfn = aero.getreqfn(altm, dalt, disa, mass, tasm, dtas, D,
                              reqrate)  # thrust required to achieve requested rate
        if reqfn < T:  # can't achieve required rate with idle thrust
            fnexc = T - D  # excess thrust
            roc = aero.getrocacc(altm, dalt, disa, mass, tasm, dtas, fnexc)
            wfe = self.get_wfe_idle(ES, AS)
        else:
            T = reqfn
            roc = reqrate
            wfe = self.get_wfe(ES, AS, T)
        return machm, tasm, casm, wfe, roc, CL, CD, T, D, bmargin

    def cruisepoint(self, mass, alt, disa, mach):
        altm = alt*aero.ft2m
        ES = EnvironmentState(altm, disa)
        AS = AircraftState(mach, 0, 0, mass, 0)
        tas = AS.TAS(ES)
        cas = AS.CAS(ES)
        wfe = 0
        bmargin = self.get_bmargin(ES, AS)
        if  cas <= self.V_MO and bmargin >= 1.0:      # the rest only makes sense if we have at least 1g margin
            Tmax = self.T_max_MCRZ(ES, AS)
            T = self.D(ES, AS)
            fnexc = Tmax-T          # excess thrust
            roc = aero.getroc(alt, disa, mass, mach, fnexc, 'Mach')
            if roc >= 100*aero.fpm2mps:
                wfe = self.get_wfe(ES, AS, T)
        return tas, cas, wfe, bmargin

def bada3climb(sched_cas, sched_mach, reqrate):
    out = open(outfilename + '.CSV', 'w')
    out.write("BADA_CLB\n")
    out.write("Climb data created for {} from bada3 file {}__.OPF\n".format(subtype, general_db))
    out.write('DISA,INITMASS,MASS,ALT,TIME,FUEL,DIST,CAS,MACH,TAS,WFE,ROC,GRAD,ALPHA,CL,CD,FN,DRAG,THR,FNEXC\n')
    out.write('DEGC,KG,KG,FT,MIN,KG,NM,KTS,,KTS,KG/H,FT/MIN,,DEG,,,DAN,DAN,,DAN\n')
    for disa in np.arange(-20, 45, 5):
        for init_mass in np.arange (base_mass, apm.max_mass+mass_step, mass_step):
            mass = init_mass  # starting at takeoff mass
            out.write('{},{:.0f},{:.3f},{:.1f},{:.6f},{:.6f},{:.6f},{:.3f},{:.8f},{:.3f},{:.5f},{:.5f},{:.8f},{:.3f},{:.9f},{:.9f},{:.4f},{:.4f},{:.4f},{:.4f}\n'.format(
                    disa, init_mass, mass, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            prev_alt, prev_mach, prev_tas, time, dist, fuel = apm.takeoff(disa, init_mass)
            while prev_alt < apm.h_MO_ft*aero.ft2m:
                next_alt = (int(round(prev_alt/aero.ft2m)/1000)+1)*1000*aero.ft2m
                if next_alt <= 2000 * aero.ft2m:
                    req_cas = apm.minCAS(mass, next_alt, disa)
                    min_roc = 10*aero.fpm2mps  # acceleration ROC
                elif next_alt <= 10000*aero.ft2m:
                    req_cas = min(max(250*aero.kt2ms, apm.minCAS(mass, next_alt, disa)), sched_cas*aero.kt2ms)
                    min_roc = 100*aero.fpm2mps
                else:
                    req_cas = sched_cas*aero.kt2ms
                    min_roc = 100*aero.fpm2mps
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
                if t > 0 and roc >= min_roc and bmargin >= 1.3:
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

def bada3descent(sched_cas, sched_mach, reqrate):
    out = open(outfilename + '.CSV', 'w')
    out.write("BADA_DSC\n")
    out.write("Descent data created for {} from bada3 file {}__.OPF\n".format(subtype, general_db))
    out.write('DISA,INITMASS,MASS,ALT,TIME,FUEL,DIST,CAS,MACH,TAS,WFE,ROC,GRAD,ALPHA,CL,CD,FN,DRAG,THR,FNEXC\n')
    out.write('DEGC,KG,KG,FT,MIN,KG,NM,KTS,,KTS,KG/H,FT/MIN,,DEG,,,DAN,DAN,,DAN\n')
    for disa in np.arange(-20, 45, 5):
        for init_mass in np.arange (base_mass, apm.max_mass+mass_step, mass_step):
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

def bada3cruise():
    out = open(outfilename + '.CSV', 'w')
    out.write("BADA_CRZ\n")
    out.write("Cruise data created for {} from bada3 file {}__.OPF\n".format(subtype, general_db))
    out.write("mass, alt, disa, mach, tas, cas, wfe, bmargin\n")
    out.write("(kg), (ft), (K), , (kt), (kt), (kg/h)\n")
    for mass in np.arange(base_mass, apm.max_mass + mass_step, mass_step):
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
    parser.add_argument('version_list')

    args = parser.parse_args()

    conn = psycopg2.connect(host='fkylabdb.flightkeys.com', database='fk_acperf', user='fkservice', password='ieDie4ieMi3aehaiZi5l')

    cur = conn.cursor()

    sqlcommand = '''
    select v.subtype, v.version_nr, v.manufacturer, r.general_db,
    r.base_mass, r.mass_step, r.base_mach, r.mach_step, 
    r.ref_lhv/2326.009143, r.ref_cg, r.min_rate_avail_clb, r.min_rate_avail_crz 
    from acperf.aircraft_version v, acperf.aircraft_release r 
    where r.aircraft_version_uid = v.uid
    and v.operator = 'BD3' and v.subtype like '{}%' and v.version_nr in ({})
    and r.rel_nr = 0 
    order by 1
    '''.format(args.subtype, args.version_list)

    cur.execute(sqlcommand)

    version_res = cur.fetchall()
    cur.close()
    conn.close()

    if len(version_res) == 0:
        print('no versions found:{}/{}\n'.format(args.subtype, args.version_list))
        quit()

    for version_row in version_res:
        subtype = version_row[0]
        vrnr = version_row[1]
        manufacturer = version_row[2]
        general_db = version_row[3]
        base_mass = version_row[4]
        mass_step = version_row[5]
        base_mach = version_row[6]
        mach_step = version_row[7]
        ref_lhv = version_row[8]
        ref_cg = version_row[9]
        min_rate_avail_clb = version_row[10] / aero.fpm2mps
        min_rate_avail_crz = version_row[11] / aero.fpm2mps

        lmax = 0  # match length

        target_path = badadef.target_root + '{}\\'.format(manufacturer)
        database_path = badadef.system_root + badadef.bada3_db

        if len(glob.glob(target_path)) == 0:    # manufacturer not found
            target_path = badadef.target_root
        if len(glob.glob(database_path)) == 0:    # manufacturer not found
            print("BADA3 database not found: {}".format(database_path))
            quit()

        apm = BADA3_CR(general_db)

        badaclimb = '{}M{}'.format(int(apm.climbCas), int(apm.climbMach*100))
        badadesc  = '{}M{}'.format(int(apm.descCas), int(apm.descMach*100))

        # delete all existing files for this version in work directory:
        fileList = glob.glob("{}_{:0>2}_*.*".format(subtype, vrnr))
        for filePath in fileList:
            try:
                os.remove(filePath)
            except:
                print("Error while deleting file : ", filePath)

        # create climb files:
        if badadef.do_climb:
            reqrate = 9E99      # max climb rate
            cas = apm.climbCas
            mach = apm.climbMach
            outfilename = "{}_{:0>2}_CLB{}M{}_000022".format(subtype, vrnr, int(cas), int(mach*100))
            print(outfilename)
            bada3climb(cas, mach, reqrate)

        if badadef.do_descent:
            reqrate = -9E99     # max descent rate
            cas = apm.descCas
            mach = apm.descMach
            outfilename = "{}_{:0>2}_DSC{}M{}_000020".format(subtype, vrnr, int(cas), int(mach * 100))
            print(outfilename)
            bada3descent(cas, mach, reqrate)

        if badadef.do_cruise_ae:
            outfilename = "{}_{:0>2}_CRZECON_000021".format(subtype, vrnr)
            print(outfilename)
            bada3cruise()

        # move created sourcefiles to target directory:
        fileList = glob.glob("{}_{:0>2}_*.*".format(subtype, vrnr))
        for filePath in fileList:
            dest_file = target_path + filePath
            try:
                if os.path.exists(dest_file):
                    os.remove(dest_file)
                shutil.move(filePath, target_path)
            except:
                print("Error while moving file {} to {}".format(filePath, target_path))