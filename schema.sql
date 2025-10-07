CREATE TABLE IF NOT EXISTS asteroids (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    designation TEXT NOT NULL,
    epoch_jd REAL,           -- Julian Date of the orbital elements
    q REAL,                  -- Perihelion distance [AU]
    e REAL,                  -- Eccentricity
    i REAL,                  -- Inclination [deg]
    om REAL,                 -- Longitude of ascending node [deg]
    w REAL,                  -- Argument of perihelion [deg]
    M REAL,                  -- Mean anomaly [deg]
    rot REAL,                -- Rotation period [hours]
    albedo REAL,             -- Geometric albedo [0-1]
    density REAL,            -- Density [g/cm³]
    H REAL,                  -- Absolute magnitude [mag]
    G REAL,                  -- Slope parameter
    x REAL,                  -- Heliocentric X [AU]
    y REAL,                  -- Heliocentric Y [AU]
    z REAL,                  -- Heliocentric Z [AU]
    vx REAL,                 -- Velocity X [AU/day]
    vy REAL,                 -- Velocity Y [AU/day]
    vz REAL,                 -- Velocity Z [AU/day]
    timestamp TEXT NOT NULL  -- retrieved at
);

