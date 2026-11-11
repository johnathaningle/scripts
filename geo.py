import math
import struct

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
    "COURSE"       : "83",
    "SPEED_KTS"    : "200",
    "LAT"          : "38.292583",
    "LNG"          : "-76.293300",
    "ALT_FT"       : "28000",
    "ECEF_X"       : None,   # MSL METERS
    "ECEF_Y"       : None,   # MSL METERS
    "ECEF_Z"       : None,   # MSL METERS
    "ENU_VEL_E_MPS": None,   # east velocity (m/s)
    "ENU_VEL_N_MPS": None,   # north velocity (m/s)
    "ENU_VEL_U_MPS": None,   # up velocity (m/s)
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
    """WGS-84 ECEF x, y, z (m) -> lat, lon (deg), height above ellipsoid."""
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


# --- speed units ---
# 1 knot = 1852 m / 3600 s, 1 mph = 1609.344 m / 3600 s (both exact)
KT_TO_MPS = 1852 / 3600
MPH_TO_MPS = 1609.344 / 3600
SPEED_TO_MPS = {"mps": 1.0, "kts": KT_TO_MPS, "mph": MPH_TO_MPS}


def convert_speed(value, unit_from, unit_to):
    """Convert a scalar speed between 'mps', 'kts', and 'mph'."""
    return value * SPEED_TO_MPS[unit_from] / SPEED_TO_MPS[unit_to]


# --- velocity in the local ENU frame ---
def course_speed_to_enu(course_deg, speed_mps, vertical_mps=0.0):
    """Course (deg clockwise from true north) + speed -> ENU velocity (m/s)."""
    c = math.radians(course_deg)
    return (
        speed_mps * math.sin(c),   # east
        speed_mps * math.cos(c),   # north
        vertical_mps,              # up
    )


def enu_to_course_speed(enu_mps):
    """ENU velocity (m/s) -> course (deg), ground speed (m/s), vertical."""
    e, n, u = enu_mps
    course = math.degrees(math.atan2(e, n)) % 360.0
    return course, math.hypot(e, n), u


# --- ENU <-> DIS world velocity ---
# DIS (IEEE 1278.1) linear velocity is an X/Y/Z vector, in m/s, in the
# world coordinate system: right-handed geocentric Cartesian, WGS-84 --
# i.e. the ECEF frame. Each component is a 32-bit float on the wire.
def enu_to_ecef_velocity(enu_mps, lat_deg, lon_deg):
    """Rotate local ENU velocity into the DIS world (ECEF) frame, m/s."""
    e, n, u = enu_mps
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    return (
        -sin_lon * e - sin_lat * cos_lon * n + cos_lat * cos_lon * u,
        cos_lon * e - sin_lat * sin_lon * n + cos_lat * sin_lon * u,
        cos_lat * n + sin_lat * u,
    )


def ecef_velocity_to_enu(ecef_mps, lat_deg, lon_deg):
    """Rotate DIS world (ECEF) velocity into the local ENU frame, m/s.

    Inverse of enu_to_ecef_velocity: the rotation is orthonormal, so
    this applies its transpose.
    """
    x, y, z = ecef_mps
    lat, lon = math.radians(lat_deg), math.radians(lon_deg)
    sin_lat, cos_lat = math.sin(lat), math.cos(lat)
    sin_lon, cos_lon = math.sin(lon), math.cos(lon)
    return (
        -sin_lon * x + cos_lon * y,
        -sin_lat * cos_lon * x - sin_lat * sin_lon * y + cos_lat * z,
        cos_lat * cos_lon * x + cos_lat * sin_lon * y + sin_lat * z,
    )


def _f32(x):
    """Round to IEEE 754 binary32, as carried in a DIS PDU."""
    return struct.unpack("f", struct.pack("f", x))[0]


def enu_to_dis_velocity(enu_mps, lat_deg, lon_deg):
    """ENU velocity (m/s) -> DIS Linear Velocity record (X, Y, Z, m/s)."""
    ecef = enu_to_ecef_velocity(enu_mps, lat_deg, lon_deg)
    return tuple(_f32(v) for v in ecef)


def dis_velocity_to_enu(dis_mps, lat_deg, lon_deg):
    """DIS Linear Velocity record (X, Y, Z, m/s) -> ENU velocity (m/s)."""
    return ecef_velocity_to_enu(dis_mps, lat_deg, lon_deg)


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
ECEF_KEYS = ("ECEF_X", "ECEF_Y", "ECEF_Z")
x_msl, y_msl, z_msl = (float(indata[k]) for k in ECEF_KEYS)

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

# --- velocity: ENU <-> raw speed <-> DIS ---

# fill in the placeholder ENU velocity: COURSE 83 at 200 kts, level
# flight (feed-supplied values pass straight through)
enu_keys = ("ENU_VEL_E_MPS", "ENU_VEL_N_MPS", "ENU_VEL_U_MPS")
if indata["ENU_VEL_E_MPS"] is None:
    speed_mps = convert_speed(float(indata["SPEED_KTS"]), "kts", "mps")
    enu = course_speed_to_enu(float(indata["COURSE"]), speed_mps)
    for key, val in zip(enu_keys, enu):
        indata[key] = val

enu_vel = tuple(float(indata[k]) for k in enu_keys)
e, n, u = enu_vel
print(f"ENU velocity (m/s): E {e:.4f}  N {n:.4f}  U {u:.4f}")

course, gs_mps, _ = enu_to_course_speed(enu_vel)
print(f"  course {course:.2f} deg  speed "
      f"{convert_speed(gs_mps, 'mps', 'kts'):.2f} kts  "
      f"{convert_speed(gs_mps, 'mps', 'mph'):.2f} mph  "
      f"{gs_mps:.2f} m/s")

dis_vel = enu_to_dis_velocity(enu_vel, lat, lon)
dx, dy, dz = dis_vel
print(f"DIS velocity (m/s): X {dx:.4f}  Y {dy:.4f}  Z {dz:.4f}")

# DIS -> ENU -> raw speed: round trip through the DIS record
back_enu = dis_velocity_to_enu(dis_vel, lat, lon)
course_rt, gs_rt_mps, _ = enu_to_course_speed(back_enu)
print(f"DIS -> ENU: course {course_rt:.2f} deg  speed "
      f"{convert_speed(gs_rt_mps, 'mps', 'kts'):.2f} kts  "
      f"{convert_speed(gs_rt_mps, 'mps', 'mph'):.2f} mph  "
      f"{gs_rt_mps:.2f} m/s")
