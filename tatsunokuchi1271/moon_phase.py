"""Moon on the night of 1271 October 17/18 (Julian): phase, age and the
surrounding principal phases, as seen from Ryukoji (Tatsunokuchi)."""
import _setup  # noqa: F401  (sets up the Astrarium core and ephemeris)
from astrarium import Ephemeris, Time, moon_info
from astrarium.observer import Site
from astrarium.bodies import moon_age_days
from astrarium.phenomena import find_moon_phases

eph = Ephemeris()
site = Site(35.311770, 139.489411, 18, "Asia/Tokyo", "Ryukoji")
JD0 = 2185579.5            # 1271-10-17 00:00 UT1 (Julian calendar)


def jst(h):
    """JST hour h on 1271 Oct 18 (Julian); negative = evening of Oct 17."""
    return Time.from_ut1_jd(JD0, (h - 9 + 24) / 24)


def julian_date_jst(t):
    """(month, day, 'hh:mm') of an instant in JST, Julian calendar."""
    jd = t.jd_ut1 + 9 / 24 + 0.5
    z = int(jd); f = jd - z
    b = z + 1524; c = int((b - 122.1) / 365.25); d = int(365.25 * c); e = int((b - d) / 30.6001)
    day = b - d - int(30.6001 * e); mon = e - 1 if e < 14 else e - 13
    hm = f * 24
    return mon, day, f"{int(hm):02d}:{round((hm % 1) * 60):02d}"


for h in (2.0, 3 + 37 / 60, 3.75, 12.0):
    t = jst(h); m = moon_info(eph, t, site)
    print(f"{h:5.2f} JST  age {moon_age_days(eph, t):5.2f} d  illum {m.illuminated_fraction:.3f}  "
          f"phase angle {m.phase_angle_deg:5.1f}  mag {m.magnitude:6.2f}  {m.phase_name_en}  "
          f"Az {float(m.az):6.1f}  El {float(m.alt):+5.1f}  bright-limb PA {m.bright_limb_pa_deg:5.1f}")

print("principal phases (JST, Julian calendar):")
for p in find_moon_phases(eph, Time.from_ut1_jd(2185565.5), Time.from_ut1_jd(2185600.5)):
    mon, day, hm = julian_date_jst(p.time)
    print(f"  {p.kind:14s} 1271-{mon:02d}-{day:02d} {hm}")
