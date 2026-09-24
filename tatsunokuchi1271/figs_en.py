import os, sys, json, math, numpy as np
import _setup
from astrarium import *
from astrarium.observer import Site
from astrarium.catalog import StarCatalog
from astrarium.meteors import solar_longitude_j2000
from astrarium.frames import airmass_kasten_young
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon, Circle
import matplotlib.transforms as mtrans
def draw_moon(ax, fig, t, R_in=0.14):
    """Moon disk with the illuminated fraction and bright-limb orientation
    of the instant, drawn in inches around its (Az, El) on the chart."""
    mi=moon_info(eph,t,site)
    az,el=float(mi.alt*0+mi.az),float(mi.alt)
    k=float(mi.illuminated_fraction)
    chi=math.radians(float(mi.bright_limb_pa_deg))          # from N through E
    H=float(mi.position.hour_angle_h)*math.pi/12; dec=math.radians(float(mi.position.dec_apparent))
    phi=math.radians(site.lat_deg)
    q=math.atan2(math.sin(H),math.tan(phi)*math.cos(dec)-math.sin(dec)*math.cos(H))  # parallactic angle
    th=chi-q                                  # from zenith, counter-clockwise (towards east = left)
    ux,uy=-math.sin(th),math.cos(th); vx,vy=-uy,ux
    pts=[]
    for tt in np.linspace(-90,90,60):
        a,b=R_in*math.cos(math.radians(tt)),R_in*math.sin(math.radians(tt)); pts.append((a*ux+b*vx,a*uy+b*vy))
    for tt in np.linspace(90,270,60):
        a,b=R_in*(2*k-1)*math.cos(math.radians(tt)),R_in*math.sin(math.radians(tt)); pts.append((a*ux+b*vx,a*uy+b*vy))
    tr=fig.dpi_scale_trans+mtrans.ScaledTranslation(az,el,ax.transData)
    ax.add_patch(Circle((0,0),R_in,transform=tr,fc="#5c5e78",ec="#c9b25a",lw=0.7,zorder=5))
    ax.add_patch(Polygon(pts,closed=True,transform=tr,fc="#ffd24d",ec="none",zorder=5.1))
    return mi
plt.rcParams["font.family"]=["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"]=False
eph=Ephemeris(); site=Site(35.311770,139.489411,18,"Asia/Tokyo","龍口寺")
JD0=2185579.5
def jst(h): return Time.from_ut1_jd(JD0,(h-9+24)/24.0)
VJ=(133.0,-(46+22/60+1/3600))
sh={s["iau_code"]:s for s in json.load(open(_setup.shower_file()))["showers"]}
def radiant(c,t):
    s=sh[c]; d=float(solar_longitude_j2000(eph,t))-s["peak_sollon"]
    return s["radiant_ra_deg"]+s["drift_ra_deg_per_sollon"]*d, s["radiant_dec_deg"]+s["drift_dec_deg_per_sollon"]*d
EOFF=1271.8-2000.0
# JST offset of local apparent time (from analysis): LAT = JST + 33.4 min
LATOFF=0.557
C=dict(vj="#d62728",moon="#e6a700",sun="#555555",ori="#1f77b4",sta="#2ca02c")

# ---------- Fig 1: altitude vs time ----------
hs=np.arange(-2.0,6.51,1/12)
def series(fn):
    return np.array([fn(jst(h)) for h in hs])
alt_vj=series(lambda t: float(observe_star(*VJ,t,site).alt))
alt_moon=series(lambda t: float(observe_planet(eph,"moon",t,site).alt))
alt_sun=series(lambda t: float(observe_planet(eph,"sun",t,site).alt))
alt_ori=series(lambda t: float(observe_star(*radiant("ORI",t),t,site).alt))
alt_sta=series(lambda t: float(observe_star(*radiant("STA",t),t,site).alt))
fig,ax=plt.subplots(figsize=(9,5.2),dpi=200)
# 辰刻 bands (定時法, 地方真太陽時)
kokus=[("I",21),("Ne",23),("Ushi",1),("Tora",3),("U",5)]
for i,(n,s) in enumerate(kokus):
    a=(s-LATOFF*1.0); a=a-24 if a>12 else a
    b=a+2
    ax.axvspan(a,b,color=("#f2f2f7" if i%2 else "#e6e9f2"),zorder=0)
    ax.text(max((a+b)/2,-1.6),78,f"{n}",ha="center",va="center",fontsize=10,color="#333")
for x,lab in [(4+32/60,"astron. dawn 04:32"),(5+57/60,"sunrise 05:57")]:
    ax.axvline(x,color="#888",ls=":",lw=1); ax.text(x+0.05,-26,lab,rotation=90,fontsize=7.5,color="#555",va="bottom")
ax.plot(hs,alt_moon,color=C["moon"],lw=2,label="Moon (age 12.0 d)")
ax.plot(hs,alt_vj,color=C["vj"],lw=2.4,label="Vela Jr. (RX J0852.0$-$4622)")
ax.plot(hs,alt_ori,color=C["ori"],lw=1.6,ls="--",label="Orionid radiant")
ax.plot(hs,alt_sta,color=C["sta"],lw=1.6,ls="--",label="S. Taurid radiant")
ax.plot(hs,np.where(alt_sun<-30,np.nan,alt_sun),color=C["sun"],lw=1.2,ls="-.",label="Sun")
ax.axhline(0,color="k",lw=0.8)
ax.annotate("moonset 03:37",xy=(3+37/60,0),xytext=(2.2,-18),arrowprops=dict(arrowstyle="->",color=C["moon"]),fontsize=8.5,color="#8a6500")
ax.annotate("Vela Jr. rises 02:43",xy=(2+43/60,0),xytext=(1.2,-27),arrowprops=dict(arrowstyle="->",color=C["vj"]),fontsize=8.5,color=C["vj"])
ax.annotate("transit 05:57\nalt. 10.9°",xy=(5+57/60,10.9),xytext=(4.6,22),arrowprops=dict(arrowstyle="->",color=C["vj"]),fontsize=8.5,color=C["vj"])
ax.set_xlim(-2,6.5); ax.set_ylim(-30,85)
ax.set_xticks(range(-2,7)); ax.set_xticklabels([f"{h%24:02d}:00" for h in range(-2,7)])
ax.set_xlabel("JST (UT+9 h), night of 1271 Oct 17/18 (Julian)"); ax.set_ylabel("Altitude [deg]")

ax.legend(loc="upper left",fontsize=8,framealpha=0.9,bbox_to_anchor=(0.0,0.93))
ax.grid(alpha=0.3)
fig.tight_layout(); fig.savefig(_setup.path("figures/en","fig1_altitude.png")); plt.close(fig)

# ---------- panorama ----------
cat=StarCatalog()
lines=cat.constellation_lines(EOFF)
cons=cat.constellations(EOFF)
def altaz(ra,dec,t):
    p=observe_star(np.asarray(ra),np.asarray(dec),t,site); return np.atleast_1d(p.az),np.atleast_1d(p.alt)
def gc_path(raz,ralt,paz,palt,ext=40,n=60):
    # great circle from radiant through point, extend beyond point by ext deg
    def v(az,alt):
        a,h=math.radians(az),math.radians(alt); return np.array([math.cos(h)*math.sin(a),math.cos(h)*math.cos(a),math.sin(h)])
    r=v(raz,ralt); p=v(paz,palt)
    ax_=np.cross(r,p); ax_/=np.linalg.norm(ax_)
    d0=math.degrees(math.acos(np.clip(np.dot(r,p),-1,1)))
    out=[]
    for d in np.linspace(d0-15,d0+ext,n):
        th=math.radians(d); w=r*math.cos(th)+np.cross(ax_,r)*math.sin(th)
        out.append((math.degrees(math.atan2(w[0],w[1]))%360, math.degrees(math.asin(np.clip(w[2],-1,1)))))
    return np.array(out)
DIR={0:"N",45:"NE",90:"E",135:"SE\n(Tatsumi)",180:"S",225:"SW",270:"W",315:"NW\n(Inui)",360:"N"}
def panorama(h,fname,title,meteors=False):
    t=jst(h); az0,az1=60,345
    fig,ax=plt.subplots(figsize=(10,4.6),dpi=200)
    sun_alt=float(observe_planet(eph,"sun",t,site).alt)
    ax.set_facecolor("#0b1026" if sun_alt<-18 else "#1b2550")
    stars=cat.apparent_stars(t,site,vmag_limit=4.8)
    for (abbr,r1,d1,r2,d2) in lines:
        a,b=altaz([r1,r2],[d1,d2],t)
        if b.max()<-12 or abs(a[0]-a[1])>180: continue
        ax.plot(a,b,color="#4a6a9a",lw=0.6,zorder=1)
    sa=np.array([s["az"] for s in stars]); sal=np.array([s["alt"] for s in stars]); sm=np.array([s["vmag"] for s in stars])
    ax.scatter(sa,sal,s=np.clip((5.2-sm)**2*2.2,0.3,40),c="white",lw=0,zorder=2)
    for c in cons:
        a,b=altaz(c["ra_deg"],c["dec_deg"],t)
        if -8<b[0]<72 and az0<a[0]<az1: ax.text(a[0],b[0],c["name_latin"],color="#8fb0e0",fontsize=6.5,ha="center",zorder=2)
    # Vela Jr
    p=observe_star(*VJ,t,site); ax.scatter([p.az],[p.alt],marker="*",s=160,c=C["vj"],edgecolors="w",lw=0.6,zorder=5)
    ax.text(float(p.az),float(p.alt)+3.5,f"Vela Jr.\nAz={float(p.az):.1f}°, El={float(p.alt):+.1f}°",color="#ff8080",fontsize=8,ha="center",va="bottom",zorder=6,bbox=dict(fc="#0b1026",ec="none",alpha=0.75,pad=1.5))
    sat=observe_planet(eph,"saturn",t,site)
    if float(sat.alt)>-10:
        ax.scatter([sat.az],[sat.alt],s=70,c="#f2d9a0",edgecolors="#a07830",lw=0.8,zorder=6)
        ax.text(float(sat.az),float(sat.alt)+2.2,f"Saturn\nAz={float(sat.az):.1f}°, El={float(sat.alt):+.1f}°",color="#f2d9a0",fontsize=8,ha="center",va="bottom",zorder=6,bbox=dict(fc="#0b1026",ec="none",alpha=0.75,pad=1.5))
    can=observe_star(95.98796,-52.69566,t,site,pmra_mas_yr=19.93,pmdec_mas_yr=23.24,parallax_mas=10.55,rv_km_s=20.3)
    if float(can.alt)>-10:
        ax.text(float(can.az)+2.0,float(can.alt)+0.6,f"Canopus\nAz={float(can.az):.1f}°, El={float(can.alt):+.1f}°",color="#e8e8f0",fontsize=7.5,ha="left",va="bottom",zorder=6,bbox=dict(fc="#0b1026",ec="none",alpha=0.75,pad=1.2))
    m=observe_planet(eph,"moon",t,site)
    if float(m.alt)>-10:
        draw_moon(ax,fig,t)
        ax.text(float(m.az),float(m.alt)+5.5,f"Moon (age 12.0 d)\nAz={float(m.az):.1f}°, El={float(m.alt):+.1f}°",color="#ffd24d",fontsize=8,ha="center",va="bottom",zorder=6,bbox=dict(fc="#0b1026",ec="none",alpha=0.75,pad=1.5))
    for code,name,col in [("ORI","ORI radiant",C["ori"]),("STA","STA radiant",C["sta"])]:
        ra,de=radiant(code,t); a,b=altaz(ra,de,t)
        ax.scatter(a,b,marker="+",s=180,c="#7fc4ff" if code=="ORI" else "#7fe07f",lw=2,zorder=5)
        if code=="ORI":
            ax.text(a[0],b[0]-3.0,name,color="#7fc4ff",fontsize=8,ha="center",va="top",zorder=5)
        else:
            ax.text(a[0]+2,b[0]+1.5,name,color="#7fe07f",fontsize=8,zorder=5)
        if meteors:
            for (pz,pa) in [(212,15),(250,30),(180,35),(135,25)]:
                g=gc_path(a[0],b[0],pz,pa,ext=12,n=30)[15:]
                ok=(g[:,1]>0)
                if ok.sum()<3: continue
                g=g[ok]; ax.annotate("",xy=g[-1],xytext=g[0],arrowprops=dict(arrowstyle="->",color="#7fc4ff" if code=="ORI" else "#7fe07f",lw=1.1,alpha=0.85))
    ax.axvline(212,color="#ffb000",ls="--",lw=1); ax.text(213,62,"Enoshima\n(Az=212°, 1.6 km)",color="#ffb000",fontsize=8)
    if meteors:
        ax.annotate("",xy=(315,8),xytext=(135,8),arrowprops=dict(arrowstyle="->",color="#ff9ad5",lw=1.2,ls=(0,(4,3))))
        ax.text(222,10,"記述：辰巳(南東) → 戌亥(北西)",color="#ff9ad5",fontsize=8)
    ax.axhspan(-15,0,color="#2b2418",zorder=3)
    ax.set_xlim(az0,az1); ax.set_ylim(-12,75)
    ticks=[90,135,180,225,270,315]; ax.set_xticks(ticks); ax.set_xticklabels([f"{DIR[k]}\n{k}°" for k in ticks],fontsize=8)
    ax.set_ylabel("Altitude [deg]"); ax.set_title(title,fontsize=10.5)
    fig.tight_layout(); fig.savefig(_setup.path("figures/en",fname)); plt.close(fig)
panorama(2.0,"fig2a_sky0200.png","(a) 02:00 JST (LAT 02:33)")
panorama(3+45/60,"fig2b_sky0345.png","(b) 03:45 JST (LAT 04:18), after moonset")

# ---------- Fig 4: compass ----------
fig=plt.figure(figsize=(5.2,5.2),dpi=200); ax=fig.add_subplot(projection="polar")
ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_ylim(0,1); ax.set_yticks([])
eto=["Ne","Ushi","Tora","U","Tatsu","Mi","Uma","Hitsuji","Saru","Tori","Inu","I"]
ax.set_xticks(np.radians(np.arange(0,360,30))); ax.set_xticklabels([f"{e}\n{i*30}°" for i,e in enumerate(eto)],fontsize=8)
def ray(az,col,lab,r=0.95,ls="-"):
    ax.plot([0,math.radians(az)],[0,r],color=col,lw=2,ls=ls); ax.text(math.radians(az),r*0.62,lab,color=col,fontsize=8,ha="center",
        bbox=dict(fc="white",ec="none",alpha=0.8,pad=1))
ray(212,"#d08a00","Enoshima 212°")
ray(135,"#c2185b","Tatsumi 135°",ls="--"); ray(315,"#c2185b","Inui 315°",ls="--")
ray(147.2,C["vj"],"Vela Jr. rise 147°",r=0.8)
ray(269.4,"#b08000","moonset 269°",r=0.8)

fig.tight_layout(); fig.savefig(_setup.path("figures/en","fig4_bearings.png")); plt.close(fig)


# ---------- Fig 5: all-sky with meteor directions ----------
def xy(az,alt):
    r=(90-np.asarray(alt))/90; a=np.radians(np.asarray(az))
    return -r*np.sin(a), r*np.cos(a)   # N up, E left
t=jst(3+45/60)
fig,ax=plt.subplots(figsize=(6.4,6.4),dpi=200); ax.set_aspect("equal"); ax.axis("off")
ax.add_patch(plt.Circle((0,0),1,color="#0b1026"))
for alt in (30,60): ax.add_patch(plt.Circle((0,0),(90-alt)/90,fill=False,ec="#34466e",lw=0.6,ls=":"))
for (abbr,r1,d1,r2,d2) in lines:
    a,b=altaz([r1,r2],[d1,d2],t)
    if b.min()<0: continue
    x,y=xy(a,b); ax.plot(x,y,color="#3d5c8c",lw=0.5)
stars=cat.apparent_stars(t,site,vmag_limit=4.5)
st=[(s["az"],s["alt"],s["vmag"]) for s in stars if s["alt"]>0]
x,y=xy([s[0] for s in st],[s[1] for s in st]); ax.scatter(x,y,s=np.clip((5-np.array([s[2] for s in st]))**2*1.6,0.3,30),c="w",lw=0)
eto=["Ne","Ushi","Tora","U","Tatsu","Mi","Uma","Hitsuji","Saru","Tori","Inu","I"]
for i,e in enumerate(eto):
    x,y=xy(i*30,-6); ax.text(x,y,e,ha="center",va="center",fontsize=10,color="#222")
for az,lab in [(0,"N"),(90,"E"),(180,"S"),(270,"W")]:
    x,y=xy(az,-14); ax.text(x,y,lab,ha="center",va="center",fontsize=9,color="#666")
cols={"ORI":"#6cb8ff","STA":"#6fdc6f"}
for code in ("ORI","STA"):
    ra,de=radiant(code,t); a,b=altaz(ra,de,t)
    X,Y=xy(a,b); ax.scatter(X,Y,marker="+",s=220,c=cols[code],lw=2.2,zorder=6)
    ax.text(X[0]+0.04,Y[0]+0.05,{"ORI":"ORI radiant","STA":"STA radiant"}[code],color=cols[code],fontsize=8,zorder=6)
    # apparent path length: trail of height extent DH seen at range r(El_P)
    # L_atm = DH / sin(El_R); projected angular length = L_atm * sin(D) / r
    DH,HM,RE=20.0,90.0,6371.0
    L_atm=DH/max(math.sin(math.radians(b[0])),0.05)
    for paz in (np.arange(0,360,45)+(0 if code=="ORI" else 22.5)):
        for palt in (12,30,55):
            g0=gc_path(a[0],b[0],paz,palt,ext=0,n=2)
            rv=np.array([math.cos(math.radians(b[0]))*math.sin(math.radians(a[0])),math.cos(math.radians(b[0]))*math.cos(math.radians(a[0])),math.sin(math.radians(b[0]))])
            pv=np.array([math.cos(math.radians(palt))*math.sin(math.radians(paz)),math.cos(math.radians(palt))*math.cos(math.radians(paz)),math.sin(math.radians(palt))])
            D=math.acos(np.clip(rv@pv,-1,1))
            if D<math.radians(3): continue
            el=math.radians(palt); rng=math.sqrt((RE+HM)**2-(RE*math.cos(el))**2)-RE*math.sin(el)
            ell=min(math.degrees(L_atm*math.sin(D)/rng),35.0)
            if ell<1.5: continue
            d0=math.degrees(D)
            def pt(d):
                ax_=np.cross(rv,pv); ax_/=np.linalg.norm(ax_); th=math.radians(d); w=rv*math.cos(th)+np.cross(ax_,rv)*math.sin(th)
                return math.degrees(math.atan2(w[0],w[1]))%360, math.degrees(math.asin(np.clip(w[2],-1,1)))
            (az0_,el0_),(az1_,el1_)=pt(d0-ell/2),pt(d0+ell/2)
            if min(el0_,el1_)<1: continue
            (x0_,y0_),(x1_,y1_)=xy(az0_,el0_),xy(az1_,el1_)
            ax.annotate("",xy=(x1_,y1_),xytext=(x0_,y0_),arrowprops=dict(arrowstyle="->",color=cols[code],lw=1.1,alpha=0.9))
# SE->NW path through zenith
x0,y0=xy(135,5); x1,y1=xy(315,5)
ax.annotate("",xy=(x1,y1),xytext=(x0,y0),arrowprops=dict(arrowstyle="->",color="#ff8ccf",lw=1.6,ls=(0,(5,3))))
ax.text(x0+0.02,y0-0.1,"Tatsumi→Inui\n(as recorded)",color="#ff8ccf",fontsize=8)
x,y=xy(212,0); ax.plot([0,x],[0,y],color="#ffb000",lw=0.9,ls="--"); ax.text(x*0.8+0.03,y*0.8,"Enoshima 212°",color="#ffb000",fontsize=8)
sat=observe_planet(eph,"saturn",t,site); X,Y=xy(float(sat.az),float(sat.alt)); ax.scatter([X],[Y],s=55,c="#f2d9a0",edgecolors="#a07830",lw=0.8,zorder=6); ax.text(X[0] if hasattr(X,"__len__") else X,(Y[0] if hasattr(Y,"__len__") else Y)+0.04,"Saturn",color="#f2d9a0",fontsize=8,ha="center",zorder=6)
p=observe_star(*VJ,t,site); X,Y=xy(float(p.az),float(p.alt)); ax.scatter([X],[Y],marker="*",s=150,c="#ff5050",edgecolors="w",lw=0.5,zorder=6)
ax.text(X+0.03,Y+0.03,"Vela Jr.",color="#ff8080",fontsize=8)
ax.set_xlim(-1.2,1.2); ax.set_ylim(-1.2,1.2)

fig.tight_layout(); fig.savefig(_setup.path("figures/en","fig3_allsky.png")); plt.close(fig)
