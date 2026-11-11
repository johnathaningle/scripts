from datetime import datetime, timedelta, timezone
from decimal import Decimal

dt_data = {
    "DAY" : "17",
    "MONTH": "7",
    "YEAR" : "2025",
    "HOUR" : "12",
    "MINUTE" : "8",
    "SECONDS": "2.499999999999993"
}

#convert the above into a datetime object
seconds = Decimal(dt_data["SECONDS"])
dt = datetime(
    year=int(dt_data["YEAR"]),
    month=int(dt_data["MONTH"]),
    day=int(dt_data["DAY"]),
    hour=int(dt_data["HOUR"]),
    minute=int(dt_data["MINUTE"]),
    second=int(seconds),
    microsecond=int(seconds % 1 * 1_000_000),
)

# --- serialize dt to string, a few different ways ---

# str(): default rendering, space-separated, microseconds included
s_default = str(dt)  # '2025-07-17 12:08:02.499999'

# isoformat(): ISO 8601, the standard for data exchange
s_iso = dt.isoformat()  # '2025-07-17T12:08:02.499999'

# RFC 3339 style with explicit UTC offset: ADS-B time is UTC, so
# attach the timezone before serializing (don't just append 'Z' to a
# naive dt). Result: '2025-07-17T12:08:02.499999+00:00'
s_rfc3339 = dt.replace(tzinfo=timezone.utc).isoformat()

# isoformat() with fixed precision via timespec (truncates, doesn't
# round): '2025-07-17T12:08:02.499' and '2025-07-17T12:08:02'
s_iso_ms = dt.isoformat(timespec="milliseconds")
s_iso_s  = dt.isoformat(timespec="seconds")

# strftime(): fully custom formats
# '2025-07-17 12:08:02.499999' and '17 Jul 2025 12:08:02 UTC'
s_full   = dt.strftime("%Y-%m-%d %H:%M:%S.%f")
s_human  = dt.strftime("%d %b %Y %H:%M:%S UTC")
# '20250717120802' (sortable, no separators)
s_compact = dt.strftime("%Y%m%d%H%M%S")

# Unix epoch seconds: what the ADS-B feeds themselves use.
# NB: .timestamp() treats a naive datetime as LOCAL time, so attach
# UTC first or the number will be wrong by your UTC offset.
# Value: 1752754082.499999 (float, so the last digits are approximate)
s_epoch = dt.replace(tzinfo=timezone.utc).timestamp()
s_epoch_ms = f"{s_epoch:.3f}"  # '1752754082.500' (ms precision, rounded)

for label, s in [
    ("str()       ", s_default),
    ("isoformat() ", s_iso),
    ("rfc3339     ", s_rfc3339),
    ("iso ms      ", s_iso_ms),
    ("iso s       ", s_iso_s),
    ("strftime()  ", s_full),
    ("human       ", s_human),
    ("compact     ", s_compact),
    ("epoch s     ", s_epoch),
    ("epoch ms    ", s_epoch_ms),
]:
    print(f"{label}: {s}")

# --- Unix epoch -> datetime, the inverse of s_epoch above ---
# Feeds carry UTC Unix epoch time: readsb/tar1090 uses epoch seconds
# (floats, milliseconds in the fraction); ADSB Exchange v2 REST uses
# integer epoch milliseconds. Parsing through Decimal keeps the
# fractional digits exact instead of smearing them into float noise.
UNIX_EPOCH = datetime(1970, 1, 1, tzinfo=timezone.utc)


def epoch_to_utc(epoch_s):
    """Unix epoch seconds (int/float/str) -> aware UTC datetime."""
    d = Decimal(str(epoch_s))
    sec = int(d)             # int() truncates toward zero, so negative
    frac = d - sec           # epochs need one flooring correction
    if frac < 0:
        sec, frac = sec - 1, frac + 1
    dt_utc = UNIX_EPOCH + timedelta(seconds=sec)
    return dt_utc.replace(microsecond=int(frac * 1_000_000))


def epoch_ms_to_utc(epoch_ms):
    """Integer epoch milliseconds -> aware UTC datetime."""
    sec, ms = divmod(int(epoch_ms), 1000)
    dt_utc = UNIX_EPOCH + timedelta(seconds=sec)
    return dt_utc.replace(microsecond=ms * 1000)


def utc_to_epoch_ms(dt_utc):
    """Aware UTC datetime -> integer epoch milliseconds (exact)."""
    delta = dt_utc - UNIX_EPOCH
    secs = delta.days * 86400 + delta.seconds
    return secs * 1000 + delta.microseconds // 1000


# --- DIS PDU header timestamps ---
# Every DIS PDU header carries a 32-bit timestamp: the least
# significant bit is the mode flag (1 = absolute, i.e. clock-synced to
# UTC; 0 = relative, simulation-defined), and the upper 31 bits count
# time since the top of the hour, rolling over every hour. The hour is
# divided into 2^31 ticks, so one tick is 3600/2^31 s, about 1.68 us.
# Decoding needs the hour itself from context.
DIS_TICK_SECONDS = 3600 / 2**31


def utc_to_dis_timestamp(dt_utc, absolute=True):
    """UTC datetime -> 32-bit DIS PDU header timestamp."""
    hour_start = dt_utc.replace(minute=0, second=0, microsecond=0)
    ticks = int((dt_utc - hour_start).total_seconds() / DIS_TICK_SECONDS)
    return (ticks << 1) | int(absolute)


def dis_timestamp_to_seconds(ts):
    """DIS timestamp -> (seconds past the top of the hour, is_absolute)."""
    return (ts >> 1) * DIS_TICK_SECONDS, bool(ts & 1)


# --- the extra conversions, demonstrated on the same dt ---
print()
dt_utc = dt.replace(tzinfo=timezone.utc)

# exact both ways: Decimal preserves the sub-millisecond digits
print(f"epoch s  -> dt: {epoch_to_utc(s_epoch)}")
print(f"epoch ms -> dt: {epoch_ms_to_utc(utc_to_epoch_ms(dt_utc))}")

# DIS PDU header timestamp for the feed timestamp
ts = utc_to_dis_timestamp(dt_utc)
sec_hour, is_abs = dis_timestamp_to_seconds(ts)
print(f"DIS timestamp: 0x{ts:08X} = {sec_hour:.4f} s past the hour, "
      f"{'absolute' if is_abs else 'relative'}")