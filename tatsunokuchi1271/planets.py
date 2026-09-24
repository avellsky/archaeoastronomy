import os,sys
import _setup
from astrarium import *
from astrarium.observer import Site
from astrarium.riseset import rise_set
eph=Ephemeris(); site=Site(35.311770,139.489411,18,"Asia/Tokyo","")
def jst(h): return Time.from_ut1_jd(2185579.5,(h-9+24)/24)   # h relative to Oct18 00 JST (negative = Oct 17)
def hm(t):
    x=(t.jd_ut1-2185579.5)*24+9-24; d="10/18" if x>=0 else "10/17"; x%=24; return f"{d} {int(x):02d}:{round((x%1)*60):02d}"
for h in (-7.0,-6.5,-6.0,-5.5,2.0,3.75,5.95):
    j=observe_planet(eph,"jupiter",jst(h),site); s=observe_planet(eph,"sun",jst(h),site)
    print(hm(jst(h)),"Jup Az %.1f El %.1f mag %.1f elong %.1f | Sun El %.1f"%(float(j.az),float(j.alt),float(j.magnitude),float(j.elongation_deg),float(s.alt)))
for b in ("jupiter","venus","mars","saturn","mercury"):
    print(b,[(e.kind,hm(e.time)) for e in rise_set(eph,b,site,jst(-24),jst(24)).events if e.kind!="transit"])
for h in (2.0,3.75):
    p=observe_planet(eph,"saturn",jst(h),site); print("Saturn",hm(jst(h)),"Az %.1f El %.1f mag %.1f"%(float(p.az),float(p.alt),float(p.magnitude)))
p=observe_planet(eph,"mercury",jst(5.0),site); print("Mercury 05:00 El %.1f"%float(p.alt))
