from flask import Blueprint, render_template, request, redirect, url_for
from services import get_phys_meteor, sim_orbit_n_conseq
# from service_helpers.api_calls import call_sentry


bp = Blueprint("routes", __name__, url_prefix="")

@bp.route("/")
def main_page():
    """ Base route """
    return render_template("index.html")


@bp.route("/search", methods=["GET", "POST"])
def search_asteroid_daesignation():
    """ Search APIs for info + parameters? POST return text data """
    if request.method == "POST":
        des = request.form.get("des")
        threat, _, _, _ = get_phys_meteor(des)
        if threat:
            orbit, r_vec, mass, entry_angle, v_mag, kinetic_energy, lat, lon = sim_orbit_n_conseq(des)
            return render_template("result.html", des=des, orbit=orbit, r_vec=r_vec, mass=mass, entry_angle=entry_angle, v_mag=v_mag, kinetic_energy=kinetic_energy, lat=lat, lon=lon)
        else:
            return render_template("no_threat.html") 
    return render_template("search.html")
    
    
@bp.route("/result")
def result_page(threat=False):
    """ Provide current most-like trajectory for PHA """
    if threat:
        return render_template("result.html")
    return render_template("no_threat.html")


@bp.route("/simulation")
def simulation():
    """ Run simulation of hitting Earth """
    return render_template("simulation.html")


@bp.route("/mitigation", methods=["GET", "POST"])
def mitigation():
    """ Choose method of mitigation PHA """
    if request.method == "POST":
        name = request.form.get("name")
        if name not in ["general", "deflection", "tractor", "nuclear"]:
            return render_template("apology.html")
        return render_template(f"{name}.html")
    return render_template("mitigation.html")