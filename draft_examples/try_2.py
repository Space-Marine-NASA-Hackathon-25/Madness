import requests
import json
from pathlib import Path
import spiceypy as spice
import pandas as pd


# --- SPICE kernels ---
DIR_KERNELS = "kernels"
KERNELS_PATHS = ["spk/de432s.bsp", "lsk/naif0012.tls", "pck/gm_de440.tpc"]
# --- SPICE body IDs ---
SUN_ID   = 10
EARTH_ID = 399
# --- Astronomical units ---
AU = spice.convrt(x=1, inunit="AU", outunit="km")
LD = 384401.0


def main():
    # spice.furnsh(get_paths_to_kernels(DIR_KERNELS, KERNELS_PATHS))
    # gm_sun, gm_earth, soi_earth, soi_earth_ld = get_constants()
    # spice.kclear()
    data = query_cad() # dist_max=soi_earth
    print(data)


def query_cad(dist_max=None, date_min="now", date_max="+60"):
    cad_url = "https://ssd-api.jpl.nasa.gov/cad.api"
    params = {
        # "date-m"in": date_min, # default "now"
        # "date-max": date_max, # default "+60"
        # "dist-max: dist_max, # AU; default for PHA is 0.05 (nominal dist)
        # "min-dist-max": 0.05, # AU; dist with uncertainty 
        "body": "Earth",
        "pha": 1, # PHA-only
        "sort": "date",
        # "diameter": 1, # None in most cases
    }
    r = requests.get(cad_url, params=params)
    r.raise_for_status()
    data = r.json()
    fields = data.get("fields", [])
    rows = data.get("data", [])
    df = pd.DataFrame(rows, columns=fields)
    return df
    
    
def get_constants():
    """ Calculates constants for the Solar System and the Earth """
    # https://naif.jpl.nasa.gov/pub/naif/misc/toolkit_docs_N0067/C/cspice/bodvcd_c.html
    _, gm_sun_pre = spice.bodvcd(bodyid=SUN_ID, item="GM", maxn=1)
    gm_sun = gm_sun_pre[0]
    _, gm_earth_pre = spice.bodvcd(bodyid=EARTH_ID, item="GM", maxn=1)
    gm_earth = gm_earth_pre[0]
    soi_earth = AU * (gm_earth / gm_sun) ** ( 2 / 5)
    soi_earth_ld = soi_earth / LD
    return gm_sun, gm_earth, soi_earth, soi_earth_ld       
    

def get_paths_to_kernels(dir_kernels, paths):
    """ Returns absolute path for each kernel file """
    # https://naif.jpl.nasa.gov/pub/naif/generic_kernels/
    BASE = Path(__file__).resolve().parent
    KERNELS = BASE / dir_kernels
    paths_to_kernels = []
    for p in paths:
        paths_to_kernels.append(str(KERNELS / p))
    return paths_to_kernels
    
    
if __name__ == "__main__":
    main()