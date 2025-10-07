import pandas as pd
import requests
import json # test
import numpy as np
import spiceypy as spice
from pathlib import Path
from datetime import datetime, timezone
from astroquery.jplhorizons import Horizons
import re

from defines import (AU, LD, GM_SUN, GM_EARTH, SOI, SOI_LD, SUN_ID, EARTH_ID, FIELDS, NUM_COLS, DIR_KERNELS, KERNELS_PATHS)

def main():
    # PHA data
    data = query_sbdb(FIELDS)
    # data_json = json.dumps(data, indent=4)
    # print(data_json) # test
    asteroids_df = make_df(data)
    print(asteroids_df) # test

    # SPICE for orbit elements and geometry calculations
    spice.furnsh(get_paths_to_kernels(DIR_KERNELS, KERNELS_PATHS))
    convert_values(asteroids_df)
    
    global GM_SUN, GM_EARTH, SOI, SOI_LD 
    GM_SUN, GM_EARTH, SOI, SOI_LD = get_soi()
    datetime_et, calendar_str, julian_date = get_datetimes()
    print(f"{datetime_et = }, {calendar_str = }, {julian_date = }") # test: only ET_TDB needed

    rows = asteroids_df.to_dict('records')
    print("\nRIGHT CODE OUTPUT:\n")
    pattern = re.compile(r'\(([^()]+)\)') # to explicitly match the designation inside SBDB for use in Horizons
    for r in rows: 
        # HORIZONS System - ephemeris data 
        match = pattern.findall(r["full_name"])[0]
        obj = obj = get_horizons(match, datetime_et)
        pha_el = get_orbital_elements(obj)
        pha_v = get_vector_elements(obj)
        print(f"\nElements for {match}: ")
        print(f"\n{pha_el = }\n{pha_v = }\n") 
        while False: # test purposes only !!!!!!!!
            calculate_crossings(pha_elements, datetime_et, SOI)
            
    # TEST CODE 1:
    print("\nTEST CODE OUTPUT:\n")
    obj = get_horizons("2022 SS2", datetime_et)
    pha_el = get_orbital_elements(obj)
    pha_v = get_vector_elements(obj)
    range_rate = pha_v[-1]
    print(f"\n{pha_el = }\n{pha_v = }\n") # test
    calculate_crossings(pha_el, datetime_et, SOI)
    
    spice.kclear()
    
    ########### TODO ###############
    ################ Simulation & visualization:
    # Cosmographia - data-based scripting for simulation 
    # Flask web-app for visuals
    ################ Impact scenarios:
    # ??? + Geological data to stimulate & visualize consequences
    ################ Consequences:
    # calculator: https://www.purdue.edu/impactearth/
    # visuals: surface involved (+ text explanation with numbers)
    ################ Mitigation strategies:
    # text-based output for the specified asteroid parameters (from previously written possible solutions based on each parameter)
    ################ Flask shell
    
    
def query_sbdb(fields):
    """ Get data from JPL SBDB for PHA """
    request_dict = {
        "fields": ",".join(fields),
        # "AND": ["spkid|EQ|54543542"], # for query a specific object (not working); for close approach 2025 QB21
        "sb-group":"pha", 
        "limit":"15" # testing purposes only!
    }
    response = requests.get('https://ssd-api.jpl.nasa.gov/sbdb_query.api', params=request_dict)
    return response.json()


def make_df(data):
    """ Makes Pandas dataframe from JSON """
    df = pd.DataFrame.from_dict(data["data"])
    df.columns = data["fields"]
    return df


def get_paths_to_kernels(dir_kernels, paths):
    """ Returns absolute path for each kernel file """
    # https://naif.jpl.nasa.gov/pub/naif/generic_kernels/
    BASE = Path(__file__).resolve().parent
    KERNELS = BASE / dir_kernels
    paths_to_kernels = []
    for p in paths:
        paths_to_kernels.append(str(KERNELS / p))
    return paths_to_kernels


def get_soi():
    """ Calculates constants for the Solar System and the Earth """
    # https://naif.jpl.nasa.gov/pub/naif/misc/toolkit_docs_N0067/C/cspice/bodvcd_c.html
    _, gm_sun_pre = spice.bodvcd(bodyid=SUN_ID, item="GM", maxn=1)
    gm_sun = gm_sun_pre[0]
    _, gm_earth_pre = spice.bodvcd(bodyid=EARTH_ID, item="GM", maxn=1)
    gm_earth = gm_earth_pre[0]
    soi_earth = AU * (gm_earth / gm_sun) ** ( 2 / 5)
    soi_earth_ld = soi_earth / LD
    return gm_sun, gm_earth, soi_earth, soi_earth_ld


def convert_values(df, num_cols=NUM_COLS):
    """ Converts several columns to the units for SPICE """
    if NUM_COLS: # avoid running if NUM_COLS is empty
        df[num_cols] = df[num_cols].apply(pd.to_numeric)
        df["q_km"] = df['q'].apply(lambda x: spice.convrt(float(x), 'AU', 'km'))
        df["inc_rad"] = np.radians(df["i"])
        df["om_rad"] = np.radians(df["om"])
        df["w_rad"] = np.radians(df["w"])
        df["ma_rad"] = np.radians(df["ma"])


def get_datetimes():
    now = datetime.now(timezone.utc)
    print(f"{now = }") # test
    calendar_str = now.strftime("%Y-%m-%d %H:%M:%S")
    datetime_et = spice.str2et(calendar_str)
    julian_date = datetime_et / 86400.0 + 2451545.0
    return datetime_et, calendar_str, julian_date


def get_horizons(pha_id, et_date):
    # ERRORS: Horizons only recognize designation (not spkid, name or full name)
    jd = 2451545.0 + et_date / 86400.0
    obj = Horizons(
        id=f"{pha_id}",
        id_type="smallbody",
        location="@sun",
        epochs=jd
    )
    return obj


def get_orbital_elements(obj):
    els_dict = dict(obj.elements()[0]) # obj.elements() return astropy.table by default
    # print(els_dict)
    orbital_elements = [
         spice.convrt(els_dict["q"],"AU", "km"),
         els_dict["e"],
         np.radians(els_dict["incl"]),
         np.radians(els_dict["Omega"]), 
         np.radians(els_dict["w"]),
         np.radians(els_dict["M"]),
         spice.utc2et(els_dict["datetime_str"]),
         GM_SUN # for calculations in SPICE
     ]
    return orbital_elements


def get_vector_elements(obj):
    els = dict(obj.vectors()[0])
    physical_elements = [
        els["H"],
        els["G"],
        els["x"],
        els["y"],
        els["z"],
        els["vx"],
        els["vy"],
        els["vz"],
        spice.convrt(els["range"], "AU", "km"),
        (els["range_rate"] * AU), # AU/day -> km/day
    ]
    return physical_elements
   

def calculate_crossing_soi(pha_id, dt_et):
    obj = get_horizons(pha_id, dt_et)
    pha_el = get_orbital_elements(obj)
    distance_pha_to_earth = get_distance_pha_to_earth(pha_el, dt_et)
    range_rate = get_vector_elements(obj)[-1]
    while range_rate < 0:
        obj = get_horizons(pha_id, dt_et)
        pha_el = get_orbital_elements(obj)
        distance_pha_to_earth = get_distance_pha_to_earth(pha_el, dt_et)        
        range_rate = get_vector_elements(obj)[-1]
        ...
        dt_et += 5
    # TODO obtain new range_rate 
    if range_rate > 0:
        pass
    elif range_rate == 0:
        hit_date = spice.et2datetime(dt_et)
        return hit_date
       
    else:
        return None
    

def get_distance_pha_to_earth(pha_el, dt_et):
    pass


def calculate_date_crossing(pha_el, dt_et, soi): # or already calculated vectors
    pass


def get_earth_state_vector(dt_et):
    earth_state_vector, _ = spice.spkgeo(targ=EARTH_ID, 
                                         et=dt_et, 
                                         ref="ECLIPJ2000", 
                                         obs=SUN_ID)
    return earth_state_vector

def get_pha_state_vector(pha_el, dt_et):
    return spice.conics(pha_el, dt_et)


def calculate_crossings(pha_el, dt_et, soi):
    # ATTENTION! Currently running correctly only if crossing SOI (means that > 2.4 LD is never finishes calculations, as well as < 2.4LD from the start)
    # ? query horizons for when negativa range_rate become positive -> min distance
    # pha_el.append(GM_SUN)
    soi_crossing_counter = 0
    bracket_min_values = []
    min_distance = np.inf
    min_et = None
    
    while soi_crossing_counter == 0:
        pha_state_vector = spice.conics(pha_el, dt_et)
        earth_state_vector, _ = spice.spkgeo(targ=EARTH_ID, 
                                             et=dt_et, 
                                             ref="ECLIPJ2000", 
                                             obs=SUN_ID)
        pha_wrt_earth_state_vector = pha_state_vector - earth_state_vector
        earth_pha_distance = spice.vnorm(pha_wrt_earth_state_vector[:3])
            
        if earth_pha_distance <= soi:
            soi_crossing_counter += 1
            print(f"First SOI crossing at UTC: {spice.et2datetime(dt_et)}")
            print(f"dist in LD: {earth_pha_distance / LD}")
        dt_et += 20
        
    pha_el_wrt_earth = spice.oscelt(state=pha_wrt_earth_state_vector,
                                    et=dt_et,
                                    mu=GM_EARTH)
    print(f"\tPerigee of PHA w.r.t. the Earth in km: {pha_el_wrt_earth[0]}")
    print(f"\tEccentricity of PHA w.r.t. the Earth: {pha_el_wrt_earth[1]}")
    
    while soi_crossing_counter == 1:
        pha_wrt_earth_state_vector = spice.conics(pha_el_wrt_earth, dt_et)
        earth_pha_distance = spice.vnorm(pha_wrt_earth_state_vector[:3])
        
        # BROKEN CODE: doesn't work at obtaining minimum (needs debugging)
        bracket_min_values.append((earth_pha_distance, dt_et))
        if len(bracket_min_values) == 3:
            t1, d1 = bracket_min_values[0]
            t2, d2 = bracket_min_values[1]
            t3, d3 = bracket_min_values[2]
            if d2 < d1 and d2 < d3:
                print(f"{bracket_min_values = }")
                t_min = quadratic_approximation(t1, d1, t2, d2, t3, d3)
                pha_wrt_earth_state_vector_at_min = spice.conics(pha_el_wrt_earth, t_min)
                d_min = spice.vnorm(pha_wrt_earth_state_vector_at_min[:3])
                print("********** Closest approach:", spice.et2utc(t_min, "ISOC", 3), "distance (km):", d_min)
                break
            bracket_min_values.pop(0)
        
        if earth_pha_distance < min_distance:
            min_distance = earth_pha_distance
            min_et = dt_et
        
        if earth_pha_distance >= soi:
            soi_crossing_counter += 1
            print(f"Second SOI crossing at UTC: {spice.et2datetime(dt_et)}")
            print(f"dist in LD: {earth_pha_distance / LD}")
        dt_et += 1

    print(f"\nClosest approach distance: {min_distance} km.\nAt time: {spice.et2utc(min_et, "ISOC", 3)}")


def quadratic_approximation(*args):
    t1, d1, t2, d2, t3, d3 = args
    a = (d1 - 2*d2 + d3) / ((t1 - t2)*(t1 - t3))
    b = (d3 - d1)/(t3 - t1) - a*(t1 + t3)
    t_min = -b/(2*a)
    print(f"{t_min = }")
    return t_min


if __name__ == "__main__":
    main()