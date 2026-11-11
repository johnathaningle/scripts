from datetime import datetime, timezone
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
s_default = str(dt) # '2025-07-17 12:08:02.499999'

# isoformat(): ISO 8601, the standard for data exchange
s_iso = dt.isoformat()                                  # '2025-07-17T12:08:02.499999'

# RFC 3339 style with explicit UTC offset: ADS-B time is UTC, so attach
# the timezone before serializing (don't just append 'Z' to a naive dt)
s_rfc3339 = dt.replace(tzinfo=timezone.utc).isoformat() # '2025-07-17T12:08:02.499999+00:00'

# isoformat() with fixed precision via timespec (truncates, doesn't round)
s_iso_ms = dt.isoformat(timespec="milliseconds")        # '2025-07-17T12:08:02.499'
s_iso_s  = dt.isoformat(timespec="seconds")             # '2025-07-17T12:08:02'

# strftime(): fully custom formats
s_full   = dt.strftime("%Y-%m-%d %H:%M:%S.%f")          # '2025-07-17 12:08:02.499999'
s_human  = dt.strftime("%d %b %Y %H:%M:%S UTC")         # '17 Jul 2025 12:08:02 UTC'
s_compact = dt.strftime("%Y%m%d%H%M%S")                 # '20250717120802' (sortable, no separators)

# Unix epoch seconds: what the ADS-B feeds themselves use.
# NB: .timestamp() treats a naive datetime as LOCAL time, so attach UTC
# first or the number will be wrong by your UTC offset.
s_epoch = dt.replace(tzinfo=timezone.utc).timestamp()  # 1752754082.499999 (float)
s_epoch_ms = f"{s_epoch:.3f}"                           # '1752754082.500' (ms precision, rounded)

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
