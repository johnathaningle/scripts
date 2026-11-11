# Geodesy & DIS (IEEE 1278.1) conversion reference
# Frames: geodetic (deg) / ECEF (m) / local ENU,NED (m/s). DIS floats are f32.
import math
import struct

# --- WGS-84 ellipsoid ---
WGS84_A = 6378137.0                     # semi-major axis (m)
WGS84_F = 1 / 298.257223563             # flattening
WGS84_B = WGS84_A * (1 - WGS84_F)       # semi-minor axis (m)
WGS84_E2 = WGS84_F * (2 - WGS84_F)      # 1st eccentricity squared
WGS84_EP2 = (WGS84_A**2 - WGS84_B**2) / WGS84_B**2   # 2nd eccentricity squared

FT_TO_M = 0.3048                        # exact, by definition of the foot

# --- inputs ---
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


# --- geodetic <-> ECEF ---
def geodetic_to_ecef(lat_deg, lon_deg, height_m):
    """Geodetic (deg) + ellipsoidal height (m) -> WGS-84 ECEF x, y, z (m)."""
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


# --- course/speed <-> ENU velocity ---
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


# --- ENU <-> ECEF vector rotation ---
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
    """Rotate DIS world (ECEF) velocity into local ENU, m/s.

    Inverse of enu_to_ecef_velocity (orthonormal transpose).
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


# --- DIS linear velocity record ---
# ECEF X/Y/Z vector in m/s, one 32-bit float per component.
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


# --- DIS orientation: psi/theta/phi <-> heading/pitch/roll ---
def euler_columns(psi, theta, phi):
    """Nose, right-wing, belly axes of Rz(psi) Ry(theta) Rx(phi).

    Angles in radians; axes are unit vectors in the frame the angles
    reference (NED for heading/pitch/roll, ECEF for DIS orientation).
    """
    cp, sp = math.cos(psi), math.sin(psi)
    ct, st = math.cos(theta), math.sin(theta)
    cr, sr = math.cos(phi), math.sin(phi)
    return (
        (cp * ct, sp * ct, -st),
        (cp * st * sr - sp * cr, sp * st * sr + cp * cr, ct * sr),
        (cp * st * cr + sp * sr, sp * st * cr - cp * sr, ct * cr),
    )


def ned_to_ecef(ned, lat_deg, lon_deg):
    """Rotate a local NED vector into ECEF (pure rotation; units preserved)."""
    n, e, d = ned
    return enu_to_ecef_velocity((e, n, -d), lat_deg, lon_deg)


def ecef_to_ned(ecef, lat_deg, lon_deg):
    """Rotate an ECEF vector into the local NED frame (units preserved)."""
    e, n, u = ecef_velocity_to_enu(ecef, lat_deg, lon_deg)
    return n, e, -u


def heading_to_dis_orientation(lat_deg, lon_deg, heading_deg,
                                pitch_deg=0.0, roll_deg=0.0):
    """Local heading/pitch/roll (deg) -> DIS psi, theta, phi (rad, f32).

    Pitch positive nose-up, roll positive right-wing-down.
    """
    axes = euler_columns(
        math.radians(heading_deg),
        math.radians(pitch_deg),
        math.radians(roll_deg),
    )
    nose, right, belly = (ned_to_ecef(a, lat_deg, lon_deg) for a in axes)
    # psi and theta come off the first column, phi off the last row
    psi = math.atan2(nose[1], nose[0])
    theta = -math.asin(max(-1.0, min(1.0, nose[2])))
    phi = math.atan2(right[2], belly[2])
    return tuple(_f32(a) for a in (psi, theta, phi))


def dis_orientation_to_heading(lat_deg, lon_deg, psi, theta, phi):
    """DIS psi, theta, phi (rad) -> heading, pitch, roll (deg).

    Inverse of heading_to_dis_orientation (gimbal-locked at pitch +-90).
    """
    nose, right, belly = (
        ecef_to_ned(a, lat_deg, lon_deg)
        for a in euler_columns(psi, theta, phi)
    )
    heading = math.degrees(math.atan2(nose[1], nose[0])) % 360.0
    pitch = math.degrees(-math.asin(max(-1.0, min(1.0, nose[2]))))
    roll = math.degrees(math.atan2(right[2], belly[2]))
    return heading, pitch, roll


# ======================================================================
# worked example on the inputs above
# ======================================================================

# altitude feet -> meters
alt_meters = float(indata["ALT_FT"]) * FT_TO_M
print(f"Alt in feet {indata['ALT_FT']} to {alt_meters:.3f} meters")

lat = float(indata["LAT"])
lon = float(indata["LNG"])

# Fill placeholder ECEF from lat/lon/MSL altitude (feed values pass through).
if indata["ECEF_X"] is None:
    indata["ECEF_X"], indata["ECEF_Y"], indata["ECEF_Z"] = geodetic_to_ecef(
        lat, lon, alt_meters
    )

x, y, z = (float(indata[k]) for k in ("ECEF_X", "ECEF_Y", "ECEF_Z"))
print(f"ECEF (M) {x:.3f} {y:.3f} {z:.3f}")

# --- velocity: ENU <-> raw speed <-> DIS ---

# Fill placeholder ENU velocity (feed values pass through).
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

# Round trip back through the DIS record.
back_enu = dis_velocity_to_enu(dis_vel, lat, lon)
course_rt, gs_rt_mps, _ = enu_to_course_speed(back_enu)
print(f"DIS -> ENU: course {course_rt:.2f} deg  speed "
      f"{convert_speed(gs_rt_mps, 'mps', 'kts'):.2f} kts  "
      f"{convert_speed(gs_rt_mps, 'mps', 'mph'):.2f} mph  "
      f"{gs_rt_mps:.2f} m/s")

# --- DIS orientation: heading/pitch/roll <-> Entity State PDU angles ---
psi, theta, phi = heading_to_dis_orientation(
    lat, lon, float(indata["COURSE"])
)
print(f"DIS orientation (rad): psi {psi:.6f}  theta {theta:.6f}  "
      f"phi {phi:.6f}")
print(f"  in degrees: psi {math.degrees(psi):8.3f}  "
      f"theta {math.degrees(theta):8.3f}  phi {math.degrees(phi):8.3f}")

# Level flight still gives a large phi: the angles are geocentric,
# not referenced to the local horizon.
h_rt, p_rt, r_rt = dis_orientation_to_heading(lat, lon, psi, theta, phi)
print(f"DIS -> heading: {h_rt:.2f} deg  pitch {p_rt:.3f} deg  "
      f"roll {r_rt:.3f} deg")
