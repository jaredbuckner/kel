#!/bin/env python3

import math
import sys

g_kerbin = 9.81
g_mun    = 1.63
g_minmus = 0.491

class Tank:
    def __init__(self, *, empty_mass, full_mass, is_solid=False):
        self.full_mass = full_mass
        self.empty_mass = empty_mass
        self.is_solid = is_solid


tank_tiny = Tank(full_mass = 0.225, empty_mass = 0.025)
tank_small = Tank(full_mass = 0.5625, empty_mass = 0.0625)
tank_large = Tank(full_mass = 4.5, empty_mass = 0.5)

class Rocket:
    def __init__(self, *, mass,
                 isp_atm, thrust_atm,
                 isp_vac, thrust_vac,
                 tank):
        self._mass = mass
        self._isp_atm = isp_atm
        self._thrust_atm = thrust_atm
        self._isp_vac = isp_vac
        self._thrust_vac = thrust_vac
        self._default_tank = tank

    def liftable_mass_atm(self, *, twr=1.15, g_body=g_kerbin,
                          tank=None, tanks=1, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks = m

        return self._thrust_atm * m / twr / g_body - self._mass * m - tanks * tank.full_mass

    def liftable_mass_vac(self, *, twr=1.15, g_body=g_mun,
                          tank=None, tanks=1, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks=m
        
        return self._thrust_vac * m / twr / g_body - self._mass * m - tanks * tank.full_mass

    def twr_loaded_atm(self, payload, *, g_body=g_kerbin,
                       tank=None, tanks=1, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks=m

        return self._thrust_atm * m / (payload + self._mass * m + tanks * tank.full_mass) / g_body
    
    def twr_loaded_vac(self, payload, *, g_body=g_mun,
                       tank=None, tanks=1, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks=m

        return self._thrust_vac * m / (payload + self._mass * m + tanks * tank.full_mass) / g_body
    
    
    def dv_atm(self, *, tank=None, tanks=1, payload=0, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks = m
        
        return g_kerbin * self._isp_atm * math.log((self._mass * m + payload + tanks * tank.full_mass) / (self._mass * m + payload + tanks * tank.empty_mass))

    def dv_vac(self, *, tank=None, tanks=1, payload=0, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks = m
        
        return g_kerbin * self._isp_vac * math.log((self._mass * m + payload + tanks * tank.full_mass) / (self._mass * m + payload + tanks * tank.empty_mass))
            

    def vtm_atm(self, payload, *, g_body=g_kerbin, tank=None, tanks=1, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks = m

        return (self.dv_atm(tank=tank, tanks=tanks, payload=payload, m=m),
                self.twr_loaded_atm(payload=payload, g_body=g_body,
                                    tank=tank, tanks=tanks, m=m),
                self._mass * m + tanks * tank.full_mass + payload)
    
    def vtm_vac(self, payload, *, g_body=g_kerbin, tank=None, tanks=1, m=1):
        if tank is None:   tank = self._default_tank
        if tank.is_solid:  tanks = m

        return (self.dv_vac(tank=tank, tanks=tanks, payload=payload, m=m),
                self.twr_loaded_vac(payload=payload, g_body=g_body,
                                    tank=tank, tanks=tanks, m=m),
                self._mass * m + tanks * tank.full_mass + payload)
    
    
    def _print_et_flow_base(self, payload, *, g_body, max_tanks = 999, max_engines=9,
                          min_dv=0, min_twr=0, tank=None, subcall):
        print("dV mtot twr")
        if tank is None:  tank = self._default_tank

        for tanks in range(1, max_tanks + 1):
            if tank.is_solid and tanks > max_engines:
                break
            
            bits = list();
            for m in range(1,max_engines + 1):
                if not tank.is_solid or tanks == m:
                    dv, twr, mtot = subcall(payload=payload, g_body=g_body, tank=tank, tanks=tanks, m=m)
                    if dv >= min_dv and twr >= min_twr:
                        bits.append(f"{dv:6.1f} {twr:5.2f} {mtot:6.2f}")
                        continue

                bits.append(f"{'  -  ':19s}")
            print(f"{tanks:3d}:"," | ".join(bits))
    
    def print_et_flow_atm(self, *args, g_body=g_kerbin, **kwargs):        
        return self._print_et_flow_base(*args, g_body=g_body, subcall=self.vtm_atm, **kwargs)

    def print_et_flow_vac(self, *args, g_body=g_mun, **kwargs):
        return self._print_et_flow_base(*args, g_body=g_body, subcall=self.vtm_vac, **kwargs)
    

def stage_configurations_atm(rocket_tank_choices, *, payload, g_body=g_kerbin,
                             max_tanks = 999, max_engines=9, min_dv=0, min_twr=0, max_dv=None):
    for name, rocket, tank in rocket_tank_choices:
        is_solid = rocket._default_tank.is_solid if tank is None else tank.is_solid
        
        for m in range(1, max_engines+1):
            for tanks in (m,) if is_solid else range(1,max_tanks+1):
                dv, twr, tmass = rocket.vtm_atm(payload=payload, g_body=g_body,
                                                tank=tank, tanks=tanks, m=m)
                if twr < min_twr:
                    break
                if max_dv is not None and dv > max_dv:
                    break
                
                if dv >= min_dv:
                    yield((name, m, tanks), (dv, twr, tmass))

def stage_configurations_vac(rocket_tank_choices, *, payload, g_body=g_mun,
                             max_tanks = 999, max_engines=9, min_dv=0, min_twr=0, max_dv=None):
    for name, rocket, tank in rocket_tank_choices:
        is_solid = rocket._default_tank.is_solid if tank is None else tank.is_solid
        
        for m in range(1, max_engines+1):
            for tanks in (m,) if is_solid else range(1,max_tanks+1):
                dv, twr, tmass = rocket.vtm_vac(payload=payload, g_body=g_body,
                                                tank=tank, tanks=tanks, m=m)
                if twr < min_twr:
                    break
                if max_dv is not None and dv > max_dv:
                    break
                
                if dv >= min_dv:
                    yield((name, m, tanks), (dv, twr, tmass))


def multi_configurations_atm(rocket_tank_choices, *, payload, stages, min_dv=0, min_twr=0,
                             g_body=g_kerbin, max_tanks = 999, max_engines=9, max_dv=None):
    for ident_a, data_a in pareto(stage_configurations_atm(rocket_tank_choices,
                                                           payload = payload,
                                                           min_dv = min_dv if stages==1 else 0,
                                                           min_twr = min_twr,
                                                           g_body=g_body,
                                                           max_tanks=max_tanks,
                                                           max_engines=max_engines,
                                                           max_dv = max_dv)):
        if stages == 1:
            yield ((ident_a,), data_a)
            continue
        
        dv_a, twr_a, tmass_a = data_a
        stage_min_dv = 0 if dv_a > min_dv else min_dv - dv_a
        stage_max_dv = None if max_dv is None else max_dv - dv_a
        stage_payload = payload + tmass_a
        for ident_b, data_b in pareto(multi_configurations_atm(rocket_tank_choices,
                                                               payload = stage_payload,
                                                               stages = stages - 1,
                                                               min_dv = stage_min_dv,
                                                               min_twr = min_twr,
                                                               g_body=g_body,
                                                               max_tanks=max_tanks,
                                                               max_engines=max_engines,
                                                               max_dv = stage_max_dv)):
            dv_b, twr_b, tmass_b = data_b
            ident = (ident_a,) + ident_b
            dv = dv_a + dv_b
            twr = twr_a if twr_a < twr_b else twr_b
            tmass = tmass_b
            yield (ident, (dv, twr, tmass))


def pareto(items, *, consider_twr = False, debug=False):
    values = list()

    for ident, data in items:
        global_pareto = True
        for idx in range(len(values) - 1, -1, -1):
            local_pareto = False
            is_dominated = True
            cdata = values[idx][1]
            
            if data[0] > cdata[0]:
                local_pareto = True
            else:
                is_dominated = False

            if consider_twr:
                if data[1] > cdata[1]:
                    local_pareto = True
                else:
                    is_dominated = False

            if data[2] < cdata[2]:
                local_pareto = True
            else:
                is_dominated = False
                
            if not local_pareto:
                global_pareto = False
                break

            if is_dominated:
                del values[idx]
                if debug: print("-", file=sys.stderr, end="")

        if global_pareto:
            values.append((ident, data))
            if debug: print("+", file=sys.stderr, end="")

        
    yield from values

    
##--------
rocket_spark = Rocket(mass=0.13, tank=tank_tiny,
                      isp_atm=265, isp_vac=320, thrust_atm=16.56, thrust_vac=20)

##--------
rocket_terrier = Rocket(mass=0.50, tank=tank_small,
                        isp_atm= 85, isp_vac=345, thrust_atm=14.78, thrust_vac=60)
rocket_dart    = Rocket(mass=1.00, tank=tank_small,
                        isp_atm=260, isp_vac=340, thrust_atm=153.53, thrust_vac=180)
rocket_swivel  = Rocket(mass=1.50, tank=tank_small,
                        isp_atm=250, isp_vac=320,  thrust_atm=167.97, thrust_vac=215)
rocket_reliant = Rocket(mass=1.25, tank=tank_small,
                        isp_atm=265, isp_vac=310, thrust_atm=205.16, thrust_vac=240)
rocket_vector  = Rocket(mass=4.00, tank=tank_small,
                        isp_atm=295, isp_vac=315, thrust_atm=936.51, thrust_vac=1000)

##--------
rocket_poodle   = Rocket(mass=1.75, tank=tank_large,
                         isp_atm= 90, isp_vac=350, thrust_atm=64.29, thrust_vac=250)
rocket_skipper  = Rocket(mass=3.00, tank=tank_large,
                         isp_atm=280, isp_vac=320, thrust_atm=568.75, thrust_vac=650)
rocket_mainsail = Rocket(mass=6.00, tank=tank_large,
                         isp_atm=285, isp_vac=310, thrust_atm=1379.03, thrust_vac=1500)

##--------
rocket_mite    = Rocket(mass=0, tank=Tank(full_mass=0.375, empty_mass=0.075, is_solid=True),
                        isp_atm=185, isp_vac=210, thrust_atm=11.012, thrust_vac=12.5)
rocket_shrimp  = Rocket(mass=0, tank=Tank(full_mass=0.875, empty_mass=0.155, is_solid=True),
                        isp_atm=190, isp_vac=215, thrust_atm=26.512, thrust_vac=30)

##--------
rocket_flea     = Rocket(mass=0, tank=Tank(full_mass=1.5, empty_mass=0.45, is_solid=True),
                         isp_atm=140, isp_vac=165, thrust_atm=162.91, thrust_vac=192)
rocket_hammer   = Rocket(mass=0, tank=Tank(full_mass=3.56, empty_mass=0.75, is_solid=True),
                         isp_atm=170, isp_vac=195, thrust_atm=197.9, thrust_vac=227)
rocket_thumper  = Rocket(mass=0, tank=Tank(full_mass=7.65, empty_mass=1.5, is_solid=True),
                         isp_atm=175, isp_vac=210, thrust_atm=250, thrust_vac=300)
rocket_kickback = Rocket(mass=0, tank=Tank(full_mass=24, empty_mass=4.5, is_solid=True),
                         isp_atm=195, isp_vac=220, thrust_atm=593.86, thrust_vac=670)

##--------
rocket_thoroughbred = Rocket(mass=0, tank=Tank(full_mass=70, empty_mass=10, is_solid=True),
                             isp_atm=205, isp_vac=230, thrust_atm=1515.217, thrust_vac=1700)
rocket_clydesdale   = Rocket(mass=0, tank=Tank(full_mass=144, empty_mass=21, is_solid=True),
                             isp_atm=210, isp_vac=235, thrust_atm=2948.936, thrust_vac=3300)


if __name__ == '__main__':
    rockets = (('swivel', rocket_swivel, None),
               #('reliant', rocket_reliant, None),
               ('flea', rocket_flea, None),
               ('hammer', rocket_hammer, None),
               ('thumper', rocket_thumper, None),
               ('kickback', rocket_kickback, None),
    )

    for stages in range(1,4):
        print(f"==== Stages: {stages} ====");
        configurator = multi_configurations_atm(rockets, stages=stages, payload=1, min_dv=3400, min_twr=1.15, max_dv=4000);
        collector = pareto(configurator, debug=False)
        scol = sorted(collector, key=lambda d: (d[1][2], -d[1][0], -d[1][1], d[0]))
        
        for reply in scol:
            print (repr(reply));
            
        print("");
