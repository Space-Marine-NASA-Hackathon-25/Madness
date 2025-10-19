import numpy as np
import spiceypy as spice
from scipy.optimize import brentq  # for root-finding

def quadratic_approximation(t1, d1, t2, d2, t3, d3):
    """
    Fit a parabola through three points (t, d) and return the time of minimum distance.
    Numerically stable for large ephemeris times.
    """
    ts = np.array([t1, t2, t3])
    ds = np.array([d1, d2, d3])
    ts_shift = ts - ts.mean()
    a, b, _ = np.polyfit(ts_shift, ds, 2)
    t_min_local = -b / (2 * a)
    return t_min_local + ts.mean()

def calculate_crossings(pha_el, dt_et_start, soi, max_iter=10**7):
    """
    Calculate first and second SOI crossings and closest approach of a PHA to Earth.
    
    Parameters:
        pha_el: heliocentric orbital elements of PHA
        dt_et_start: starting ephemeris time
        soi: Earth's sphere of influence radius (km)
        max_iter: maximum number of iterations for scanning (safety)
    """
    LD = 384400  # lunar distance in km
    GM_EARTH = 398600.4418  # km^3/s^2
    EARTH_ID = 399
    SUN_ID = 10

    dt_et = dt_et_start
    soi_crossings = []
    min_distance = np.inf
    min_et = None
    bracket = []

    def distance_from_earth(t, pha_elements):
        pha_vec = spice.conics(pha_elements, t)
        earth_vec, _ = spice.spkgeo(EARTH_ID, t, "ECLIPJ2000", SUN_ID)
        pha_wrt_earth = pha_vec - earth_vec
        return spice.vnorm(pha_wrt_earth[:3])

    # 1️⃣ Scan for first SOI crossing
    step = 20  # initial coarse step
    for i in range(max_iter):
        dist = distance_from_earth(dt_et, pha_el)
        if dist <= soi:
            # Refine first crossing using root-finding
            t_entry = brentq(lambda t: distance_from_earth(t, pha_el) - soi,
                             dt_et - step, dt_et)
            soi_crossings.append(t_entry)
            print(f"First SOI crossing at UTC: {spice.et2datetime(t_entry)}")
            print(f"Distance in LD: {distance_from_earth(t_entry, pha_el)/LD}")
            break
        dt_et += step
    else:
        raise RuntimeError("First SOI crossing not found")

    # 2️⃣ Compute PHA elements w.r.t Earth at entry
    pha_vec_entry = spice.conics(pha_el, soi_crossings[0])
    pha_el_wrt_earth = spice.oscelt(pha_vec_entry, soi_crossings[0], GM_EARTH)
    print(f"\tPerigee w.r.t Earth (km): {pha_el_wrt_earth[0]}")
    print(f"\tEccentricity w.r.t Earth: {pha_el_wrt_earth[1]}")

    # 3️⃣ Scan for closest approach inside SOI
    dt_et = soi_crossings[0]
    step = 1  # small step for min distance
    for i in range(max_iter):
        pha_vec = spice.conics(pha_el_wrt_earth, dt_et)
        earth_dist = spice.vnorm(pha_vec[:3])

        # Store 3-point bracket for quadratic fit
        bracket.append((dt_et, earth_dist))
        if len(bracket) == 3:
            (t1, d1), (t2, d2), (t3, d3) = bracket
            if d2 < d1 and d2 < d3:
                t_min = quadratic_approximation(t1, d1, t2, d2, t3, d3)
                pha_vec_min = spice.conics(pha_el_wrt_earth, t_min)
                d_min = spice.vnorm(pha_vec_min[:3])
                min_distance = d_min
                min_et = t_min
                print("********** Closest approach **********")
                print(f"UTC: {spice.et2utc(t_min, 'ISOC', 3)}, Distance: {d_min} km")
                # Increase step after closest approach to speed up exit scanning
                dt_et = t_min + 10
                step = 10
                break
            bracket.pop(0)

        if earth_dist < min_distance:
            min_distance = earth_dist
            min_et = dt_et

        dt_et += step

    # 4️⃣ Scan for second SOI crossing (exit)
    for i in range(max_iter):
        dist = distance_from_earth(dt_et, pha_el_wrt_earth)
        if dist >= soi:
            # Refine using root-finding
            t_exit = brentq(lambda t: distance_from_earth(t, pha_el_wrt_earth) - soi,
                            dt_et - step, dt_et)
            soi_crossings.append(t_exit)
            print(f"Second SOI crossing at UTC: {spice.et2datetime(t_exit)}")
            print(f"Distance in LD: {distance_from_earth(t_exit, pha_el_wrt_earth)/LD}")
            break
        dt_et += step

    print(f"\nClosest approach distance: {min_distance} km at UTC: {spice.et2utc(min_et, 'ISOC', 3)}")
    return soi_crossings, min_et, min_distance