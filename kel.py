#!/bin/env python3

import math
import sys

## Surface gravity
g_kerbin = 9.81
g_mun    = 1.63
g_minmus = 0.491

## Surface pressure
p_kerbin = 101.325
p_mun    = 0
p_minmus = 0
p_space  = 0

class Tank:
    def __init__(self, *, empty_mass, full_mass, is_solid=False):
        self.full_mass = full_mass
        self.empty_mass = empty_mass
        self.is_solid = is_solid


tank_tiny = Tank(full_mass = 0.225, empty_mass = 0.025)
tank_small = Tank(full_mass = 0.5625, empty_mass = 0.0625)
tank_large = Tank(full_mass = 4.5, empty_mass = 0.5)
tank_xlarge = Tank(full_mass = 20.25, empty_mass = 2.25)


class Rocket:
    def __init__(self, *, name, mass,
                 isp_atm, thrust_atm,
                 isp_vac, thrust_vac,
                 tank, does_gimbal=True, is_radial=False):
        self._name = name
        self._does_gimbal = does_gimbal
        self._is_radial = is_radial
        self._mass = mass
        self._isp_atm = isp_atm
        self._thrust_atm = thrust_atm
        self._isp_vac = isp_vac
        self._thrust_vac = thrust_vac
        self._tank = tank

    # Change name and tank type
    def clone(self, *, name, tank):
        return Rocket(name=name, mass=self._mass, tank=tank,
                      does_gimbal=self._does_gimbal, is_radial=self._is_radial,
                      isp_atm=self._isp_atm, thrust_atm=self._thrust_atm,
                      isp_vac=self._isp_vac, thrust_vac=self._thrust_vac)
    
    
    def name(self):
        return self._name

    def isp(self, p_body):
        return self._isp_atm if p_body else self._isp_vac

    def thrust(self, p_body):
        return self._thrust_atm if p_body else self._thrust_vac
    
    def tank(self):
        return self._tank

    def does_gimbal(self):
        return self._does_gimbal

    def is_radial(self):
        return self._is_radial
    
    def liftable_mass(self, *, twr=1.15, g_body=g_kerbin, p_body = p_kerbin,
                      tanks=1, m=1):
        if self._tank.is_solid:  tanks = m
        
        return self.thrust(p_body) * m / twr / g_body - self._mass * m - tanks * self._tank.full_mass
    
    def twr_loaded(self, payload, *, g_body=g_kerbin, p_body=p_kerbin,
                   tanks=1, m=1):
        if self._tank.is_solid:  tanks=m

        return self.thrust(p_body) * m / (payload + self._mass * m + tanks * self._tank.full_mass) / g_body
    
    def dv(self, *, tanks=1, payload=0, m=1, p_body=p_kerbin):
        if self._tank.is_solid:  tanks = m

        mass_ratio = ((self._mass * m + payload + tanks * self._tank.full_mass) /
                      (self._mass * m + payload + tanks * self._tank.empty_mass))
        
        ## g_kerbin is the unit correction, not the local g
        return g_kerbin * self.isp(p_body) * math.log(mass_ratio)            
    
    def vtm(self, payload, *, g_body=g_kerbin, p_body=p_kerbin, tanks=1, m=1):
        if self._tank.is_solid:  tanks = m
        
        return (self.dv(p_body=p_body, tanks=tanks, payload=payload, m=m),
                self.twr_loaded(payload=payload, g_body=g_body, p_body=p_body,
                                tanks=tanks, m=m),
                self._mass * m + tanks * self._tank.full_mass + payload)    
    
    def print_et_flow(self, payload, *, g_body=g_kerbin, p_body=p_kerbin,
                      max_tanks = 999, max_engines=9,
                      min_dv=0, min_twr=0):
        print("dV mtot twr")
        
        for tanks in range(1, max_tanks + 1):
            if self._tank.is_solid and tanks > max_engines:
                break
            
            bits = list();
            for m in range(1,max_engines + 1):
                if not self._tank.is_solid or tanks == m:
                    dv, twr, mtot = self.vtm(payload=payload, g_body=g_body, p_body=p_body,
                                             tanks=tanks, m=m)
                    if dv >= min_dv and twr >= min_twr:
                        bits.append(f"{dv:6.1f} {twr:5.2f} {mtot:6.2f}")
                        continue

                bits.append(f"{'  -  ':19s}")
            print(f"{tanks:3d}:"," | ".join(bits))
    

def stage_configurations(rocket_choices, *, payload, g_body=g_kerbin, p_body=p_kerbin,
                         max_tanks = 999, max_engines=9, min_dv=0, min_twr=0,
                         max_dv=None, max_mass=None):
    for rocket in rocket_choices:
        me = 8 if rocket.is_radial() and max_engines > 8 else max_engines
        for m in range(2 if rocket.is_radial() else 1, me+1):
            for tanks in (m,) if rocket.tank().is_solid else range(1,max_tanks+1):
                dv, twr, tmass = rocket.vtm(payload=payload, g_body=g_body, p_body=p_body,
                                            tanks=tanks, m=m)
                if twr < min_twr:
                    break
                if max_dv is not None and dv > max_dv:
                    break
                if max_mass is not None and tmass > max_mass:
                    break
                
                if dv >= min_dv:
                    yield((rocket, m, tanks), (dv, twr, tmass))


def multi_configurations(rocket_choices, *, payload, stages, min_dv=0, min_twr=0,
                         g_body=g_kerbin, p_body=p_kerbin,
                         max_tanks = 999, max_engines=9,
                         max_dv=None, max_mass=None):
    for ident_a, data_a in pareto(stage_configurations(rocket_choices,
                                                       payload = payload,
                                                       min_dv = min_dv if stages==1 else 0,
                                                       min_twr = min_twr,
                                                       g_body=g_body,
                                                       p_body=p_body,
                                                       max_tanks=max_tanks,
                                                       max_engines=max_engines,
                                                       max_dv = max_dv,
                                                       max_mass = max_mass)):
        if stages == 1:
            yield ((ident_a, ), data_a)
            continue
        
        dv_a, twr_a, tmass_a = data_a
        stage_min_dv = 0 if dv_a > min_dv else min_dv - dv_a
        stage_max_dv = None if max_dv is None else max_dv - dv_a
        stage_max_mass = None if max_mass is None else max_mass - tmass_a
        stage_payload = tmass_a
        for ident_b, data_b in pareto(multi_configurations(rocket_choices,
                                                           payload = stage_payload,
                                                           stages = stages - 1,
                                                           min_dv = stage_min_dv,
                                                           min_twr = min_twr,
                                                           g_body=g_body,
                                                           p_body=p_body,
                                                           max_tanks=max_tanks,
                                                           max_engines=max_engines,
                                                           max_dv = stage_max_dv,
                                                           max_mass = stage_max_mass)):
            dv_b, twr_b, tmass_b = data_b
            ident = (ident_a, ) + ident_b
            dv = dv_a + dv_b
            twr = twr_a if twr_a < twr_b else twr_b
            tmass = tmass_b
            yield (ident, (dv, twr, tmass))


def print_stages(stage_iter, *, payload, g_body=g_kerbin, p_body=p_kerbin):
    tdv = 0
    for rocket, m, tanks in stage_iter:
        dv, twr, tmass = rocket.vtm(payload, g_body=g_body, p_body=p_body, tanks=tanks, m=m)
        tdv += dv
        print(f'  {payload:6.2f}T + {m:1d}*{rocket.name():<12s} +{tanks:3d} fuel => {tmass:6.2f}T@{twr:4.1f}TWR {tdv:6.1f}dV ({dv:+6.1f})')
        payload = tmass
    

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
rocket_ant    = Rocket(mass=0.02, name='ant', tank=tank_tiny, does_gimbal=False,
                       isp_atm= 80, isp_vac=315, thrust_atm=0.51, thrust_vac=2)
rocket_spider = Rocket(mass=0.02, name='spider', tank=tank_tiny, is_radial=True,
                       isp_atm=260, isp_vac=290, thrust_atm=1.79, thrust_vac=2)
rocket_twitch = Rocket(mass=0.08, name='twitch', tank=tank_tiny, is_radial=True,
                       isp_atm=275, isp_vac=290, thrust_atm=15.17, thrust_vac=16.00)
rocket_spark  = Rocket(mass=0.13, name='spark', tank=tank_tiny,
                       isp_atm=265, isp_vac=320, thrust_atm=16.56, thrust_vac=20)

##--------
rocket_terrier = Rocket(mass=0.50, name='terrier', tank=tank_small,
                        isp_atm= 85, isp_vac=345, thrust_atm=14.78, thrust_vac=60)
rocket_thud    = Rocket(mass=0.90, name='thud', tank=tank_small, is_radial=True,
                        isp_atm=275, isp_vac=305, thrust_atm=108.197, thrust_vac=120);
rocket_dart    = Rocket(mass=1.00, name='dart', tank=tank_small, does_gimbal=False,
                        isp_atm=260, isp_vac=340, thrust_atm=153.53, thrust_vac=180)
rocket_swivel  = Rocket(mass=1.50, name='swivel', tank=tank_small,
                        isp_atm=250, isp_vac=320,  thrust_atm=167.97, thrust_vac=215)
rocket_reliant = Rocket(mass=1.25, name='reliant', tank=tank_small, does_gimbal=False,
                        isp_atm=265, isp_vac=310, thrust_atm=205.16, thrust_vac=240)
rocket_vector  = Rocket(mass=4.00, name='vector', tank=tank_small,
                        isp_atm=295, isp_vac=315, thrust_atm=936.51, thrust_vac=1000)

##--------
rocket_poodle   = Rocket(mass=1.75, name='poodle', tank=tank_large,
                         isp_atm= 90, isp_vac=350, thrust_atm=64.29, thrust_vac=250)
rocket_skipper  = Rocket(mass=3.00, name='skipper', tank=tank_large,
                         isp_atm=280, isp_vac=320, thrust_atm=568.75, thrust_vac=650)
rocket_mainsail = Rocket(mass=6.00, name='mainsail', tank=tank_large,
                         isp_atm=285, isp_vac=310, thrust_atm=1379.03, thrust_vac=1500)

##--------
rocket_rhino    = Rocket(mass=9.00, name='rhino', tank=tank_xlarge,
                         isp_atm=205, isp_vac=340, thrust_atm=1205.88, thrust_vac=2000)
rocket_mammoth  = Rocket(mass=15, name='mammoth', tank=tank_xlarge,
                         isp_atm=295, isp_vac=315, thrust_atm=3746.03, thrust_vac=4000)

##--------
rocket_sepratron = Rocket(mass=0, name='sepratron', tank=Tank(full_mass=0.07, empty_mass=0.01, is_solid=True),
                          does_gimbal=False, is_radial=True,
                          isp_atm=118, isp_vac=154, thrust_atm=13.79, thrust_vac=18)

##--------
rocket_mite    = Rocket(mass=0, name='mite', tank=Tank(full_mass=0.375, empty_mass=0.075, is_solid=True),
                        does_gimbal=False, isp_atm=185, isp_vac=210, thrust_atm=11.012, thrust_vac=12.5)
rocket_shrimp  = Rocket(mass=0, name='shrimp', tank=Tank(full_mass=0.875, empty_mass=0.155, is_solid=True),
                        does_gimbal=False, isp_atm=190, isp_vac=215, thrust_atm=26.512, thrust_vac=30)

##--------
rocket_flea     = Rocket(mass=0, name='flea', tank=Tank(full_mass=1.5, empty_mass=0.45, is_solid=True),
                         does_gimbal=False, isp_atm=140, isp_vac=165, thrust_atm=162.91, thrust_vac=192)
rocket_hammer   = Rocket(mass=0, name='hammer', tank=Tank(full_mass=3.56, empty_mass=0.75, is_solid=True),
                         does_gimbal=False, isp_atm=170, isp_vac=195, thrust_atm=197.9, thrust_vac=227)
rocket_thumper  = Rocket(mass=0, name='thumper', tank=Tank(full_mass=7.65, empty_mass=1.5, is_solid=True),
                         does_gimbal=False, isp_atm=175, isp_vac=210, thrust_atm=250, thrust_vac=300)
rocket_kickback = Rocket(mass=0, name='kickback', tank=Tank(full_mass=24, empty_mass=4.5, is_solid=True),
                         does_gimbal=False, isp_atm=195, isp_vac=220, thrust_atm=593.86, thrust_vac=670)

##--------
rocket_thoroughbred = Rocket(mass=0, name='thoroughbred', tank=Tank(full_mass=70, empty_mass=10, is_solid=True),
                             isp_atm=205, isp_vac=230, thrust_atm=1515.217, thrust_vac=1700)
rocket_clydesdale   = Rocket(mass=0, name='clydesdale', tank=Tank(full_mass=144, empty_mass=21, is_solid=True),
                             isp_atm=210, isp_vac=235, thrust_atm=2948.936, thrust_vac=3300)

##--------
rocket_nerv = Rocket(mass=3, name='nerv', tank=tank_small, does_gimbal=False,
                     isp_atm=185, isp_vac=800, thrust_atm=13.88, thrust_vac=60)

##--------
rocket_dawn = Rocket(mass=0.25, name='dawn', tank=Tank(full_mass=0.10, empty_mass=0.02), does_gimbal=False,
                     isp_atm=100, isp_vac=4200, thrust_atm=0.05, thrust_vac=2)
                     


rockets_start = (rocket_flea,)
rockets_basic = rockets_start + (rocket_swivel, rocket_hammer)
rockets_general = rockets_basic + (rocket_reliant, rocket_thumper)
rockets_advanced = rockets_general + (rocket_terrier, rocket_thud)
rockets_heavy = rockets_advanced + (rocket_poodle, rocket_skipper, rocket_kickback)
rockets_heavier = rockets_heavy + (rocket_mainsail, rocket_thoroughbred)
rockets_very_heavy = rockets_heavier + (rocket_vector, rocket_rhino, rocket_clydesdale, rocket_mammoth)

rockets_add_prop = (rocket_spark, rocket_ant, rocket_mite, rocket_shrimp)
rockets_add_prec = (rocket_sepratron, rocket_twitch, rocket_spider)
rockets_add_nerv = (rocket_nerv,)
rockets_add_dart = (rocket_dart,)
rockets_add_dawn = (rocket_dawn,)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Kerbel Engine Lift Calculator')

    parser.add_argument('--payload', type=float, metavar='tons', required=True)
    rocketsel = parser.add_mutually_exclusive_group(required=True)
    rocketsel.add_argument('--start', dest='rockets', const=rockets_start, action='store_const')
    rocketsel.add_argument('--basic', dest='rockets', const=rockets_basic, action='store_const')
    rocketsel.add_argument('--general', dest='rockets', const=rockets_general, action='store_const')
    rocketsel.add_argument('--advanced', dest='rockets', const=rockets_advanced, action='store_const')
    rocketsel.add_argument('--heavy', dest='rockets', const=rockets_heavy, action='store_const')
    rocketsel.add_argument('--heavier', dest='rockets', const=rockets_heavier, action='store_const')
    rocketsel.add_argument('--very_heavy', dest='rockets', const=rockets_very_heavy, action='store_const')
    parser.add_argument('--with_prop', action='store_true')
    parser.add_argument('--with_prec', action='store_true')
    parser.add_argument('--with_nerv', action='store_true')
    parser.add_argument('--with_dart', action='store_true')
    parser.add_argument('--with_dawn', action='store_true')
    bodysel = parser.add_mutually_exclusive_group(required=True)
    bodysel.add_argument('--kerbin', dest='body', const=g_kerbin, action='store_const')
    bodysel.add_argument('--mun', dest='body', const=g_mun, action='store_const')
    bodysel.add_argument('--minmus', dest='body', const=g_minmus, action='store_const')
    psel = parser.add_mutually_exclusive_group(required=True)
    psel.add_argument('--vac', dest='press', const=p_space, action='store_const')
    psel.add_argument('--atm', dest='press', const=p_kerbin, action='store_const')
    parser.add_argument('--twr', type=float, metavar='value', default=1.15)
    dvsel = parser.add_mutually_exclusive_group(required=True)
    dvsel.add_argument('--dvmin', type=float, metavar='m/s', default=None)
    dvsel.add_argument('--dvtgt', type=float, metavar='m/s', default=None)
    parser.add_argument('--dvmax', type=float, metavar='m/s', default=None)
    parser.add_argument('--gimbal', action='store_true')
    parser.add_argument('--show', type=int, default=None)
    parser.add_argument('--stagesmax', type=int, default=5)
    
    args = parser.parse_args()

    payload = args.payload
    g_body  = args.body
    p_body  = args.press
    min_twr = args.twr
    show    = args.show
    
    tgt_dv  = args.dvtgt if args.dvmin is None else args.dvmin * 1.1
    max_dv  = (args.dvmax if args.dvmax is not None else
               None if tgt_dv == 0 else
               tgt_dv * 2.0)

    rockets = args.rockets
    if args.with_prop:
        rockets += rockets_add_prop
    if args.with_prec:
        rockets += rockets_add_prec
    if args.with_nerv:
        rockets += rockets_add_nerv
    if args.with_dart:
        rockets += rockets_add_dart
    if args.with_dawn:
        rockets += rockets_add_dawn
    
    if args.gimbal:
        rockets = tuple(r for r in rockets if r.does_gimbal())

    best_prev = None
    for stages in range(1,args.stagesmax + 1):
        print(f"==== Stages: {stages} ====");
        configurator = multi_configurations(rockets,
                                            stages=stages,
                                            payload=payload,
                                            g_body=g_body,
                                            p_body = p_body,
                                            min_dv=tgt_dv,
                                            min_twr=min_twr,
                                            max_dv=max_dv);
        collector = pareto(configurator, debug=False)
        scol = sorted(collector, key=lambda d: (d[1][2], -d[1][0], -d[1][1], d[0]))

        best = None
        for stages, totals in scol if show is None else scol[:show]:
            dv, twr, tmass = totals
            if best is None or tmass < best:
                best = tmass
            print(f'{tmass:6.2f}T @ {dv:6.1f}dV : {twr:4.1f} TWR')
            print_stages(stages, payload=payload, g_body=g_body, p_body=p_body);

        if best is not None and (best_prev is None or best < best_prev):            
            print("");
            best_prev = best
        elif best is not None:
            break
