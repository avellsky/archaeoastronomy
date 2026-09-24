# archaeoastronomy

Positional-astronomy reconstructions of historical celestial events.

| Directory | Study |
|---|---|
| [`tatsunokuchi1271/`](tatsunokuchi1271/) | The luminous object of the Tatsunokuchi persecution (AD 1271 October 18, Julian): a supernova in Vela Jr. or a fireball? |

---

## tatsunokuchi1271 — the Tatsunokuchi luminous object (AD 1271)

In his letter *Shuju Onfurumai Gosho*, Nichiren records that his execution at
Tatsunokuchi (Kamakura, Japan) before dawn on 1271 October 18 (Julian) was
halted when an object "shining like the Moon" and "like a ball" crossed the
sky from the direction of Enoshima, from the south-east (*tatsumi*) to the
north-west (*inui*). These scripts reconstruct the sky over the execution
ground (Ryūkōji, 35.311770° N, 139.489411° E, 18 m) and test two hypotheses:

* a supernova producing the remnant **Vela Jr.** (RX J0852.0−4622), and
* a **fireball** (bright meteor), e.g. of the Taurid complex.

Paper: Hayakawa, H., Motizuki, Y., & Abe, S., *The luminous object of the
Tatsunokuchi persecution (AD 1271): a supernova in Vela Jr. or a fireball?*
(in preparation for PASJ).

### Main results

| Quantity | Value (JST = UT1 + 9 h, Julian calendar) |
|---|---|
| New Moon before the event | 1271-10-06 02:52 (→ lunar month IX day 1) |
| Moon at 02:00 | waxing gibbous, age 11.96 d, 93 % lit, −12.1 mag, Az 255°, El +19° |
| Moonset | 03:37 (03:42 with a fixed standard altitude of −0.833°) |
| Full Moon | 1271-10-20 06:04 |
| Astronomical dawn / sunrise | 04:32 / 05:57 |
| Vela Jr. rise / transit / set | 02:43 (Az 147°) / 05:57 (El 10.9°) / 09:10 |
| Direction of Enoshima | Az 211–212°, 1.6 km |
| Solar longitude (J2000) | 221.5° (≈ present-day November 3; near the S. Taurid maximum) |
| Planets above the horizon | Saturn only (rose 01:16; Az 97°, El 30° at 03:45) |
| Moon-bright fireball (m ≤ −12) at one site | ~0.6–1.3 per year of dark, clear nights (EN catalogue) |

### Scripts

| Script | Output | Content |
|---|---|---|
| `analysis.py` | `results/results.json` | calendar, rise/set/transit, twilight, coordinates of Vela Jr., half-hourly positions, meteor radiants, Enoshima bearing |
| `sensitivity.py` | stdout | ΔT sensitivity, moonset conventions, radiant-side criterion for SE→NW motion, supernova magnitudes and extinction |
| `moon_phase.py` | stdout | Moon phase, age, bright-limb angle and principal phases |
| `planets.py` | stdout | visibility of Jupiter, Saturn and the other planets |
| `apex.py` | stdout | positions of the apex and antihelion sporadic sources |
| `fireball_rate.py` | stdout | fireball frequency from the European Fireball Network catalogue (Borovička et al. 2022, A&A 667, A157; downloaded from CDS on first run) |
| `figs_en.py` | `figures/en/` | figures of the paper (English) |
| `figs_ja.py` | `figures/ja/` | the same figures with Japanese labels |

Text outputs of the scripts are stored in `results/`.

### Requirements

* Python ≥ 3.10 with `numpy`, `pyerfa`, `jplephem`, `matplotlib` (`pip install -r requirements.txt`)
* the computation core of **Astrarium** (the `astrarium` Python package; Abe 2026,
  Astrarium v2.1). It is **not included** in this repository; point
  `ASTRARIUM_SRC` to the directory that contains the `astrarium` package.
  The core uses ERFA (IAU SOFA) for IAU 2006/2000A precession–nutation,
  aberration, light deflection, parallax and refraction.
* a JPL ephemeris covering AD 1271 — **DE406** (or DE441). DE421 does not
  cover the thirteenth century. Download e.g. from
  <https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/> and either place it at
  `tatsunokuchi1271/data/de406.bsp` or set `ASTRARIUM_EPHEMERIS`.
* Japanese figure labels use the Hiragino Sans font (macOS).

### Usage

```bash
cd tatsunokuchi1271
export ASTRARIUM_SRC=/path/to/astrarium/src
export ASTRARIUM_EPHEMERIS=/path/to/de406.bsp
python3 analysis.py        # -> results/results.json
python3 sensitivity.py
python3 moon_phase.py
python3 planets.py
python3 apex.py
python3 fireball_rate.py
python3 figs_en.py         # -> figures/en/*.png
python3 figs_ja.py         # -> figures/ja/*.png
```

### Notes on the computation

* Dates before 1582 are Julian-calendar dates. Instants are passed to the
  code as Julian Dates (`Time.from_ut1_jd`), because ERFA's calendar
  routines assume the proleptic Gregorian calendar.
* ΔT = 552 s (Espenak & Meeus 2006). An error of ±300 s shifts all clock
  times by ∓5 min and does not change the conclusions.
* Meteor-shower radiants and activity ranges are present-day values
  (Jenniskens 1994, 2006); the evolution of the streams over 750 yr is not
  modelled.
* Apparent meteor path lengths in the all-sky figure follow
  ℓ ≃ ΔH sin D / (r sin El_R), with ΔH = 20 km and H = 90 km.
* Fireball rates assume the uniform coverage of 7 × 10⁵ km² quoted by
  Borovička et al. (2022) for fireballs brighter than −10 mag, over 2 yr,
  with 0.25 mag per airmass extinction and a fireball height of 80 km.

### Citation

If you use these scripts, please cite the paper above and

* Abe, S. 2026, Astrarium, version 2.1, planetarium application for iOS and macOS.
