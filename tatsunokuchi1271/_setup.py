"""Common setup: locate the Astrarium computation core and a long-span JPL ephemeris.

Environment variables
---------------------
ASTRARIUM_SRC        path to the `src/` directory that contains the `astrarium` package
                     (not needed if `astrarium` is already importable)
ASTRARIUM_EPHEMERIS  path to a JPL SPK kernel covering AD 1271 (DE406 or DE441);
                     defaults to data/de406.bsp or data/de441.bsp next to this file
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

if os.environ.get("ASTRARIUM_SRC"):
    sys.path.insert(0, os.environ["ASTRARIUM_SRC"])

if "ASTRARIUM_EPHEMERIS" not in os.environ:
    for name in ("de406.bsp", "de441.bsp"):
        p = os.path.join(HERE, "data", name)
        if os.path.exists(p):
            os.environ["ASTRARIUM_EPHEMERIS"] = p
            break
    else:
        sys.exit("Set ASTRARIUM_EPHEMERIS to a JPL ephemeris covering AD 1271 "
                 "(DE406 or DE441), or put de406.bsp in data/ (see README).")


def shower_file():
    """Meteor-shower list shipped with the Astrarium data directory."""
    from astrarium.paths import data_dir
    return os.path.join(data_dir(), "meteor_showers.json")


def path(*parts):
    """Path relative to this directory (outputs go to results/, figures/)."""
    p = os.path.join(HERE, *parts)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    return p
