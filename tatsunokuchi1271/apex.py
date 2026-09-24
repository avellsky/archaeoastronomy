import os,sys,math,numpy as np
import _setup
from astrarium import *
from astrarium.observer import Site
from astrarium.frames import ecliptic_to_icrs
from astrarium.meteors import solar_longitude_j2000
eph=Ephemeris(); site=Site(35.311770,139.489411,18,"Asia/Tokyo","")
for h in [2.0,3.0,3.75,4.25]:
    t=Time.from_ut1_jd(2185579.5,(h-9+24)/24)
    sl=float(solar_longitude_j2000(eph,t))
    ra,de=ecliptic_to_icrs(sl-90,0.0,None)
    p=observe_star(float(ra),float(de),t,site)
    ra2,de2=ecliptic_to_icrs(sl-180+15,0.0,None)   # antihelion source approx (λ-λ☉≈195°)
    q=observe_star(float(ra2),float(de2),t,site)
    print(h, "apex az %.1f alt %.1f"%(float(p.az),float(p.alt)), "antihelion az %.1f alt %.1f"%(float(q.az),float(q.alt)))
