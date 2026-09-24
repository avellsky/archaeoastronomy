import os,sys,math,json,numpy as np
import _setup
from astrarium import *
from astrarium.observer import Site
from astrarium.riseset import rise_set
from astrarium.frames import ecliptic_to_icrs, airmass_kasten_young
from astrarium.meteors import solar_longitude_j2000
eph=Ephemeris(); site=Site(35.311770,139.489411,18,"Asia/Tokyo","")
JD0=2185579.5
def jst(h,dts=0.0): return Time.from_ut1_jd(JD0,(h-9+24)/24+dts/86400)
def hm(t,dts=0.0):
    x=((t.jd_ut1-dts/86400-JD0)*24+9)%24; return f"{int(x):02d}:{(x%1)*60:04.1f}"
VJ=(133.0,-(46+22/60+1/3600))
print("t2 deltaT",jst(2).delta_t)
# DeltaT sensitivity: shift the instant by dts in TT relative to UT (equivalent to deltaT error)
# emulate by evaluating at UT1 shifted and reporting against original clock -> shows the effect of a UT offset on positions
for dts in (-600,-300,0,300,600):
    t0=jst(0,dts); t1=jst(12,dts)
    ms=[e for e in rise_set(eph,"moon",site,t0,t1).events if e.kind=="set"][0]
    vr=[e for e in rise_set(eph,"star",site,t0,t1,star=VJ).events if e.kind=="rise"][0]
    p=observe_star(*VJ,jst(3.75,dts),site)
    print("dUT(s)",dts,"moonset",hm(ms.time,dts),"VJ rise",hm(vr.time,dts),"VJ alt 03:45 %.2f"%float(p.alt))
# note: shifting UT1 by dts with fixed sky mostly shifts clock; better: print sidereal effect
# moonset with alternative h0
for h0,lab in [(None,"astrarium default (topocentric, upper limb, refraction)"),(-0.833,"h0=-0.833 (no parallax)"),(-0.566,"centre, refraction only"),(0.0,"geometric centre 0")]:
    ev=[e for e in rise_set(eph,"moon",site,jst(0),jst(24),h0=h0).events]
    print(lab,[(e.kind,hm(e.time)) for e in ev if e.kind!="transit"])
# direction criterion: SE component of radiant unit vector
def svec(az,alt):
    a,h=math.radians(az),math.radians(alt); return np.array([math.cos(h)*math.sin(a),math.cos(h)*math.cos(a),math.sin(h)])
SE=svec(135,0)
sh={s["iau_code"]:s for s in json.load(open(_setup.shower_file()))["showers"]}
def rad(c,t):
    s=sh[c]; d=float(solar_longitude_j2000(eph,t))-s["peak_sollon"]
    return s["radiant_ra_deg"]+s["drift_ra_deg_per_sollon"]*d, s["radiant_dec_deg"]+s["drift_dec_deg_per_sollon"]*d
out={}
for h in (3+37/60,3.75,4.0,4.25,4.5):
    t=jst(h); sl=float(solar_longitude_j2000(eph,t)); row={}
    srcs={c:rad(c,t) for c in ("ORI","STA","NTA")}
    srcs["apex"]=tuple(float(x) for x in ecliptic_to_icrs(sl-90,0.0,None))
    srcs["antihelion"]=tuple(float(x) for x in ecliptic_to_icrs(sl+180-15+30,0.0,None))  # placeholder, replaced below
    srcs["antihelion"]=tuple(float(x) for x in ecliptic_to_icrs(sl+195-360,0.0,None))
    for k,(ra,de) in srcs.items():
        p=observe_star(ra,de,t,site); v=svec(float(p.az),float(p.alt))
        row[k]=(round(float(p.az),1),round(float(p.alt),1),round(float(v@SE),2))
    out[f"{h:.3f}"]=row
for k,v in out.items(): print(k,v)
# SN apparent magnitude
for M,lab in [(-19.3,"Ia"),(-18.0,"CC bright"),(-16.8,"IIP mean"),(-16.0,"CC faint")]:
    for d in (385,750):
        m=M+5*math.log10(d/10); print(lab,d,"pc m=%.1f"%m)
for h in (3.0,3+37/60,4.0,4.5,5+57/60):
    p=observe_star(*VJ,jst(h),site); X=float(airmass_kasten_young(float(p.alt)))
    print("h %.2f alt %.2f X %.1f ext(k=0.25) %.1f sun %.1f"%(h,float(p.alt),X,0.25*X,float(observe_planet(eph,"sun",jst(h),site).alt)))
