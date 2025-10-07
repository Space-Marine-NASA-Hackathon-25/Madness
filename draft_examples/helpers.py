import requests
from astroquery.jplhorizons import Horizons
from astropy.time import Time
from astropy.table import Talbe
import astropy.units as u


SSD_API = "https://ssd-api.jpl.nasa.gov"


def query_sentry(des=None):
    """Get Sentry VIs for a designation."""
    r = requests.get(
        f"{SSD_API}/sentry.api",
        params={
            
        })
    r.raise_for_status()
    return r.json()


def query_cad(start_date, end_date, des=None):
    """Return close-approach data within a date range (for a designation)."""
    r = requests.get(
        f"{SSD_API}/cad.api",
        params={
            
        })
    r.raise_for_status()
    return r.json()


def query_sbdb(des, fields=[
    "full_name", "H", "G", "density", "albedo",
    "e", "a", "q", "i", "om", "w", "ma", "tp",
    "per", "n", "ad", "moid_ld"]):
    """Return JSON from SBDB for a designation."""
    r = requests.get(
        f"{SSD_API}/sbdb.api",
        params={
            "des": f"{des}",
            "fields": ",".join(fields),
            # "sb-group":"pha", 
            # "limit":"15"
        })
    r.raise_for_status()
    return r.json()


def query_horizons(des, epoch_iso):
    """Return JSON from Horizons for a designation."""
    obj = Horizons(
        id=f"{des}",
        id_type="smallbody",
        location="@sun",
        epoch=Time(epoch_iso).jd
    )
    return obj.vectors(), obj.elements()


def json2df(data):
    """Convertss JSON data to Pandas DF""" # astropy.table
    table = ...
    return table


def get_constants():
    pass
