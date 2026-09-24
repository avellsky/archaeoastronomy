import os, sys, json, math, numpy as np, erfa
import _setup
from astrarium import *
from astrarium.observer import Site
from astrarium.riseset import rise_set, twilight_times
from astrarium.bodies import moon_age_days
from astrarium.phenomena import previous_moon_phase
from astrarium.meteors import solar_longitude_j2000, MeteorCalendar
from astrarium.frames import icrs_to_galactic, icrs_to_ecliptic, airmass_kasten_young
eph=Ephemeris()
site=Site(35.311770,139.489411,18,"Asia/Tokyo","龍口寺")
JD0=2185579.5  # Julian cal 1271-10-17 12:00 UT? -> JD .5 = 0h UT on 10-17
def jst(h):  # JST hour on 1271-10-18 (Julian) ; JST h -> UT h-9 on 10-17 when h<9
    return Time.from_ut1_jd(JD0, (h-9+24)/24.0) if True else None
def hm(t):
    x=(t.jd_ut1 - JD0)*24 + 9 - 24   # JST hours on Oct18
    x=x % 24 if x<0 or x>=24 else x
    h=int(x); m=round((x-h)*60)
    if m==60: h,m=h+1,0
    return f"{h:02d}:{m:02d}"
R={}
t2=jst(2.0)
R["calendar"]=dict(jd_ut1_0h_jst=float(jst(0).jd_ut1), deltaT=float(t2.delta_t),
  greg=[int(x) for x in erfa.jd2cal(jst(2).jd_ut1+9/24,0)[:3]],
  weekday=["月","火","水","木","金","土","日"][int(math.floor(jst(2).jd_ut1+9/24+0.5))%7])
# new moon before
nm=previous_moon_phase(eph,t2,0); fm=None
R["new_moon"]=dict(jd_ut1=float(nm.jd_ut1), jst_date_offset_days=float((nm.jd_ut1+9/24)-(JD0+0.5)), )
nm_jst=nm.jd_ut1+9/24+0.5
y,mo,d,f=erfa.jd2cal(nm.jd_ut1+9/24,0)
# Julian calendar date of new moon
def jul(jd):
    Z=int(jd+0.5);F=jd+0.5-Z;A=Z;B=A+1524;C=int((B-122.1)/365.25);D=int(365.25*C);E=int((B-D)/30.6001)
    day=B-D-int(30.6001*E)+F; mon=E-1 if E<14 else E-13; yr=C-4716 if mon>2 else C-4715
    return yr,mon,day
R["new_moon"]["julian_jst"]=jul(nm.jd_ut1+9/24)
R["moon_age_2h"]=float(moon_age_days(eph,t2)); R["moon_age_0h_jst"]=float(moon_age_days(eph,jst(0)))
R["moon_age_21h_jst_prev(=12h UT?)"]=float(moon_age_days(eph,Time.from_ut1_jd(JD0,0.5)))
# rise/set local day
t0=jst(0.0); t1=jst(24.0)
rs={}
for b in ["sun","moon","mercury","venus","mars","jupiter","saturn"]:
    arc=rise_set(eph,b,site,t0,t1)
    rs[b]={e.kind:hm(e.time)+(f" (az {e.az:.1f})" if e.kind!="transit" else f" (alt {e.alt:.1f})") for e in arc.events}
    bv=observe_planet(eph,b,t2,site); rs[b]["mag"]=None if b in("sun","moon") else round(float(bv.magnitude),1)
VJ=(133.0, -(46+22/60+1/3600))
arc=rise_set(eph,"star",site,t0,t1,star=VJ)
rs["VelaJr"]={e.kind:hm(e.time)+(f" (az {e.az:.1f})" if e.kind!="transit" else f" (alt {e.alt:.2f})") for e in arc.events}
R["riseset"]=rs
tw=twilight_times(eph,site,t0,t1); R["twilight"]={k:[hm(x) for x in v] for k,v in tw.items()}
# Vela coords
p=observe_star(*VJ,t2,site)
lon,lat=p.ecliptic(apparent=True); gl,gb=icrs_to_galactic(*VJ)
el_j2000=icrs_to_ecliptic(*VJ,None)
b1950=erfa.fk524(math.radians(VJ[0]),math.radians(VJ[1]),0,0,0,0)
R["vela_coords"]=dict(ra_app_h=float(p.ra_apparent)/15, dec_app=float(p.dec_apparent), ecl_app=(float(lon),float(lat)),
  ecl_j2000=(float(el_j2000[0]),float(el_j2000[1])), gal=(float(gl),float(gb)), b1950=(math.degrees(b1950[0])/15, math.degrees(b1950[1])))
# time series
sl=float(solar_longitude_j2000(eph,t2)); R["sollon_j2000"]=sl
cal=json.load(open(_setup.shower_file()))["showers"]
sh={s["iau_code"]:s for s in cal}
def radiant(code, sollon):
    s=sh[code]; d=sollon-s["peak_sollon"]
    return s["radiant_ra_deg"]+s["drift_ra_deg_per_sollon"]*d, s["radiant_dec_deg"]+s["drift_dec_deg_per_sollon"]*d
R["radiants_j2000"]={c:radiant(c,sl) for c in ("ORI","STA","NTA")}
R["showers"]={c:dict(start=sh[c]["activity_start_sollon"],end=sh[c]["activity_end_sollon"],peak=sh[c]["peak_sollon"],v=sh[c]["vinf_km_s"]) for c in ("ORI","STA","NTA")}
rows=[]
for h in np.arange(0,6.01,0.5):
    t=jst(h); r={"jst":f"{int(h):02d}:{int(round((h%1)*60)):02d}"}
    pv=observe_star(*VJ,t,site); r["vj"]=(float(pv.az),float(pv.alt))
    m=observe_planet(eph,"moon",t,site); r["moon"]=(float(m.az),float(m.alt))
    s=observe_planet(eph,"sun",t,site); r["sun_alt"]=float(s.alt)
    for c in ("ORI","STA"):
        rp=observe_star(*radiant(c,float(solar_longitude_j2000(eph,t))),t,site); r[c]=(float(rp.az),float(rp.alt))
    # local apparent solar time
    r["LAT"]=float((observe_body(eph,"sun",t,site).hour_angle_h+12)%24); r["LMT"]=(h+ (139.489411-135)/15)%24
    rows.append(r)
R["series"]=rows

# Enoshima: island centre & 江島神社辺津宮 (approx.)
from math import radians,degrees,atan2,sin,cos,acos
def bearing(lat2,lon2,lat1=35.311770,lon1=139.489411):
    dl=radians(lon2-lon1);p1,p2=radians(lat1),radians(lat2)
    return degrees(atan2(sin(dl)*cos(p2),cos(p1)*sin(p2)-sin(p1)*cos(p2)*cos(dl)))%360, 6371*acos(sin(p1)*sin(p2)+cos(p1)*cos(p2)*cos(dl))
R["enoshima"]={"island_centre":bearing(35.2992,139.4800),"hetsu_miya":bearing(35.3004,139.4806)}
# Vela transit dark? sun alt at transit
json.dump(R,open(_setup.path("results","results.json"),"w"),ensure_ascii=False,indent=1,default=str)
print(json.dumps(R,ensure_ascii=False,indent=1,default=str))
