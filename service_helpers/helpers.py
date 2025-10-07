from datetime import datetime, timezone
import spiceypy as spice
from pathlib import Path


DIR_KERNELS = "kernels"
KERNELS_PATHS = ["spk/de432s.bsp", "lsk/naif0012.tls", "pck/gm_de440.tpc"]

def au_in_km():
    """ Return high-precision AU in km """
    for k in get_paths_to_kernels():
        spice.furnsh(k)
    AU_in_km = spice.convrt(x=1, inunit="AU", outunit="km")
    spice.kclear()
    return AU_in_km


def get_paths_to_kernels(dir_kernels=DIR_KERNELS, paths=KERNELS_PATHS):
    """ Returns absolute path for each kernel file """
    # https://naif.jpl.nasa.gov/pub/naif/generic_kernels/
    BASE = Path(__file__).resolve().parent.parent
    KERNELS = BASE / dir_kernels
    return [str(KERNELS / p) for p in paths]


def sanitize_input(query):
    """ Return clear string without special characters """
    pass


def json2dict(js_obj):
    """ Make dict from JSON """
    dict_obj = ...
    return dict_obj


def get_datetimes():
    """ Return datetime in different formats """
    now = datetime.now(timezone.utc)
    calendar_str = now.strftime("%Y-%m-%d %H:%M:%S")
    datetime_et = spice.str2et(calendar_str)
    julian_date = datetime_et / 86400.0 + 2451545.0
    return calendar_str, datetime_et, julian_date