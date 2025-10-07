import requests
import json
from astroquery.jplhorizons import Horizons
from service_helpers.helpers import get_datetimes, get_paths_to_kernels, au_in_km
import spiceypy as spice
import numpy as np
from astropy.time import Time

SUN_ID   = 10
EARTH_ID = 399
AU_in_km = au_in_km()

FIELDS = [
    "full_name", "H", "G", "density", "albedo",
    "e", "a", "q", "i", "om", "w", "ma", "tp", "per", "n", "ad", "moid_ld"]


def call_sentry(designation):
    """ Get data from Sentry """
    response = requests.get(
        "https://ssd-api.jpl.nasa.gov/sentry.api", 
        params={"des": designation})
    data = response.json()
    if "error" in data:
        return
    mass = float(data["summary"]["mass"])
    diameter = float(data["summary"]["diameter"])
    energy = float(data["summary"]["energy"])
    # print(json.dumps(data, indent = 4))
    return mass, diameter, energy

# print(call_sentry("2008 JL3"))

def call_horizons(designation):
    """ Get data from JPL Horizons for the asteroid """
    spice.furnsh(get_paths_to_kernels())
    cal_str, _, jd = get_datetimes()
    obj = Horizons(
        id=f"{designation}",
        id_type="smallbody",
        location="@sun",
        epochs=jd
    )
    x, y, z, vx, vy, vz = get_vector_elements(obj)[2:8]
    # epoch = get_orbital_elements(obj)[-2]
    spice.kclear()    
    return x, y, z, vx, vy, vz, cal_str


def get_orbital_elements(obj):
    els_dict = dict(obj.elements()[0]) # obj.elements() return astropy.table by default
    # print(els_dict)
    _, gm_sun_pre = spice.bodvcd(bodyid=SUN_ID, item="GM", maxn=1)
    gm_sun = gm_sun_pre[0]
    orbital_elements = [
        spice.convrt(els_dict["q"],"AU", "km"),
        els_dict["e"],
        np.radians(els_dict["incl"]),
        np.radians(els_dict["Omega"]), 
        np.radians(els_dict["w"]),
        np.radians(els_dict["M"]),
        #  els_dict["datetime_str"][5:],
        spice.utc2et(els_dict["datetime_str"]),
        gm_sun # for calculations in SPICE
     ]
    return orbital_elements


def get_vector_elements(obj):
    els = dict(obj.vectors()[0])
    print(f"/n from horizons: {els}/n")
    physical_elements = [
        spice.convrt(els["H"], "AU", "km"),
        spice.convrt(els["G"], "AU", "km"),
        spice.convrt(els["x"], "AU", "km"),
        spice.convrt(els["y"], "AU", "km"),
        spice.convrt(els["z"], "AU", "km"),
        spice.convrt(els["vx"], "AU", "km"),
        spice.convrt(els["vy"], "AU", "km"),
        spice.convrt(els["vz"], "AU", "km"),
        spice.convrt(els["range"], "AU", "km"),
        (els["range_rate"] * AU_in_km), # AU/day -> km/day
    ]
    return physical_elements

# print(call_horizons("2008 JL3"))


def call_smdb(designation):
    """ Get data from JPL SBDB for the asteroid: https://ssd-api.jpl.nasa.gov/doc/sbdb.html  """
    response = requests.get(
        "https://ssd-api.jpl.nasa.gov/sbdb.api", 
        params={"des": designation})
    return json.dumps(response.json(), indent = 4)
