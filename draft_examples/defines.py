import spiceypy as spice

# --- Astronomical units ---
AU = spice.convrt(x=1, inunit="AU", outunit="km")
LD = 384401.0

# --- SPICE body IDs ---
SUN_ID   = 10
EARTH_ID = 399

# --- Placeholders for constant result using SPICE ---
GM_SUN   = None
GM_EARTH = None
SOI      = None
SOI_LD   = None

# --- For querying SBDB ---
FIELDS = [
        # "spkid", 
        "full_name", 
        # "name", "pdes", "class", 
        # "orbit_id", 
        # "epoch",
        # "sb_used", 
        "H", "G", "density", "albedo", 
]
NUM_COLS = ["e", "a", "q", "i", "om", "w", "ma", "tp", "per", "n", "ad", "moid_ld"]
FIELDS.extend(NUM_COLS)

# --- SPICE kernels ---
DIR_KERNELS = "kernels"
KERNELS_PATHS = ["spk/de432s.bsp", "lsk/naif0012.tls", "pck/gm_de440.tpc"]