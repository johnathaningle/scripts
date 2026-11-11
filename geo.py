import math

# WGS-84 ellipsoid constants
WGS84_A = 6378137.0                     # semi-major axis (m)
WGS84_F = 1 / 298.257223563             # flattening
WGS84_B = WGS84_A * (1 - WGS84_F)       # semi-minor axis (m)
WGS84_E2 = WGS84_F * (2 - WGS84_F)      # 1st eccentricity squared
WGS84_EP2 = (WGS84_A**2 - WGS84_B**2) / WGS84_B**2   # 2nd eccentricity squared

FT_TO_M = 0.3048                        # exact, by definition of the foot

# Geoid undulation N: height of the EGM2008 geoid above the WGS-84
# ellipsoid, at the input location below (GeographicLib GeoidEval,
# interpolated RMS error ~1 mm; EGM96 gives -36.155 m, EGM84 -35.999 m).
# Relation:  h = N + H   (HAE = MSL + N).  N is negative here because the
# geoid surface sits ~36 m BELOW the ellipsoid on the US mid-Atlantic
# coast.  This value is location-specific: recompute it (GeoidEval, or
# pyproj with the EGM2008 grid) if the input coordinates change.
GEOID_UNDULATION_M = -35.8850

indata = {
    "COURSE": "83",
    "LAT"   : "38.292583",
    "LNG"   : "-76.293300",
    "ALT_FT": "28000",
    "ECEF_X": None,   # MSL METERS
    "ECEF_Y": None,   # MSL METERS
    "ECEF_Z": None,   # MSL METERS
}


def geodetic_to_ecef(lat_deg, lon_deg, height_m):
    """Geodetic (deg) + height (m) -> WGS-84 ECEF x, y, z (m).

    height_m is the height above the WGS-84 ellipsoid.
    """
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    # prime vertical radius of curvature
    n = WGS84_A / math.sqrt(1 - WGS84_E2 * sin_lat**2)
    return (
        (n + height_m) * cos_lat * math.cos(lon),
        (n + height_m) * cos_lat * math.sin(lon),
        ((1 - WGS84_E2) * n + height_m) * sin_lat,
    )


def ecef_to_geodetic(x, y, z):
    """WGS-84 ECEF x, y, z (m) -> lat, lon (deg), height above ellipsoid (m)."""
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    # Bowring's closed-form solution for geodetic latitude
    theta = math.atan2(WGS84_A * z, WGS84_B * p)
    lat = math.atan2(
        z + WGS84_EP2 * WGS84_B * math.sin(theta) ** 3,
        p - WGS84_E2 * WGS84_A * math.cos(theta) ** 3,
    )
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    # height formula that stays stable near the poles
    h = p * cos_lat + z * sin_lat - WGS84_A * math.sqrt(
        1 - WGS84_E2 * sin_lat**2
    )
    return math.degrees(lat), math.degrees(lon), h


# Convert altitude from feet to meters
alt_meters = float(indata["ALT_FT"]) * FT_TO_M
print(f"Alt in feet {indata['ALT_FT']} to {alt_meters:.3f} meters")

lat = float(indata["LAT"])
lon = float(indata["LNG"])

# fill in the placeholder ECEF values: the aircraft position as
# lat/lon/MSL-altitude expressed in ECEF meters (real feed values, if
# present, pass straight through)
if indata["ECEF_X"] is None:
    indata["ECEF_X"], indata["ECEF_Y"], indata["ECEF_Z"] = geodetic_to_ecef(
        lat, lon, alt_meters
    )

# Convert ECEF Mean Sea Level in Meters to Height Above Ellipsoid
x_msl, y_msl, z_msl = (float(indata[k]) for k in ("ECEF_X", "ECEF_Y", "ECEF_Z"))

# 1) ECEF -> geodetic. ECEF->geodetic always returns an ELLIPSOIDAL
#    height, but these coordinates were built from an MSL altitude, so
#    what comes back is really H: height above mean sea level.
glat, glon, h_msl = ecef_to_geodetic(x_msl, y_msl, z_msl)

# 2) MSL -> HAE: apply the geoid undulation, h = N + H. With N ~ -36 m
#    here the aircraft sits lower w.r.t. the ellipsoid than its MSL
#    altitude alone suggests.
h_hae = h_msl + GEOID_UNDULATION_M

# 3) back to ECEF, now referenced to the ellipsoid (HAE)
x_hae, y_hae, z_hae = geodetic_to_ecef(glat, glon, h_hae)

print(f"ECEF MSL (M) {x_msl:.3f} {y_msl:.3f} {z_msl:.3f} "
      f"to ECEF HAE (M) {x_hae:.3f} {y_hae:.3f} {z_hae:.3f}")
print(f"  {h_msl:.3f} m MSL -> {h_hae:.3f} m HAE "
      f"(geoid undulation N = {GEOID_UNDULATION_M:+.4f} m)")
