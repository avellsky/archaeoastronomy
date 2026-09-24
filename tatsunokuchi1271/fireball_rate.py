# Fireball flux from the European Fireball Network catalogue (Borovicka et al. 2022, A&A 667, A157)
# and the expected rate of Moon-bright fireballs seen from one site.
import math, numpy as np
recs=[]
import os, urllib.request, _setup
CAT=_setup.path("data","en_catalog.dat")
if not os.path.exists(CAT):
    urllib.request.urlretrieve("https://cdsarc.cds.unistra.fr/ftp/J/A+A/667/A157/catalog.dat", CAT)
for l in open(CAT):
    try: m=float(l[539:544])
    except: continue
    if m==0: continue
    recs.append((l[:15].strip(),l[17:27],float(l[60:68]),m,l[611:614].strip()))
code=np.array([r[0] for r in recs]); date=np.array([r[1] for r in recs])
L=np.array([r[2] for r in recs]); M=np.array([r[3] for r in recs]); S=np.array([r[4] for r in recs])
YEARS,AREA=2.0,7e5
out={}
print("cumulative counts / density (per km^2 per year of network operation):")
for lim in (-10,-11,-12,-13,-14):
    n=int((M<=lim).sum()); out[lim]=n; print(f"  M<={lim}: N={n:3d}  rho={n/YEARS/AREA:.2e}")
ms=np.array([-10,-11,-12,-13]); ns=np.array([(M<=m).sum() for m in ms]); b,a=np.polyfit(ms,np.log10(ns),1)
print(f"fit: log10 N(<=M) = {a:.2f} + {b:.3f} M")
def rho(Mlim,Mcut):
    if Mlim<Mcut: Mlim=Mcut   # no fireballs brighter than Mcut
    return 10**(a+b*Mlim)/YEARS/AREA - (10**(a+b*Mcut)/YEARS/AREA if Mcut>-99 else 0)
RE=6371.0; K=0.25
def airmass(el):
    return 1/(math.sin(math.radians(el))+0.50572*(el+6.07995)**-1.6364)
def site_rate(m0,Hkm=80.0,Mcut=-14.5,Rmax=1200,dmax=None):
    tot=0.0; dR=0.5
    for R in np.arange(dR/2,Rmax,dR):
        d=math.hypot(R,Hkm); el=math.degrees(math.atan2(Hkm-R**2/(2*RE),R))
        if el<=0.5: break
        if dmax and d>dmax: break
        Mlim=m0-5*math.log10(d/100.0)-K*airmass(el)
        if Mlim<Mcut: break
        tot+=2*math.pi*R*dR*(10**(a+b*Mlim)-10**(a+b*Mcut))/YEARS/AREA
    return tot
for Mcut,lab in ((-14.5,"truncated at M=-14.5 (brightest observed -14.1)"),(-17.0,"power law extrapolated to M=-17")):
    r=site_rate(-12,Mcut=Mcut); rl=site_rate(-12,Mcut=Mcut,dmax=150)
    print(f"apparent m<=-12, {lab}: {r:.2f}/yr (1 per {1/r:.1f} yr); within d<=150 km: {rl:.3f}/yr (1 per {1/rl:.0f} yr)")
sel=(L>=200)&(L<=245)
print("lambda 200-245 (12.5% of year): M<=-10",int(((M<=-10)&sel).sum()),"/",int((M<=-10).sum()),"; M<=-12",int(((M<=-12)&sel).sum()),"/",int((M<=-12).sum()))
tau=np.isin(S,["STA","NTA","TAU","LTA","FTA","TAT"])
print("Taurid-assigned: all",int(tau.sum()),"; M<=-10",int(((M<=-10)&tau).sum()),"; M<=-12",int(((M<=-12)&tau).sum()))
for i in np.where((M<=-12)&(abs(L-221.5)<3))[0]: print("near 221.5:",code[i],date[i],L[i],S[i] or "sporadic",M[i])
