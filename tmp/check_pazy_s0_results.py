import csv
import math
import statistics
import sys
from bisect import bisect_right
from pathlib import Path


path = Path(sys.argv[1])
step = int(sys.argv[2]) if len(sys.argv) > 2 else 1
speeds = [1e-3, *range(step, 61, step)]
expected_count = 3 * len(speeds)
with path.open(newline="") as stream:
    rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(stream)]

print(f"{len(rows)}/{expected_count} rows, {sum(row['uvlm_converged'] for row in rows):.0f} converged")
assert len(rows) == expected_count
for row in rows:
    assert row["uvlm_converged"] == 1.0
    assert all(math.isfinite(value) for value in row.values())
    assert math.isclose(
        row["dynamic_pressure_pa"], 0.5 * row["density_kg_m3"] * row["airspeed_m_s"] ** 2
    )

for quantity, tolerance in (
    ("geometry_residual", 1e-5),
    ("load_residual", 1e-4),
    ("wake_force_residual", 1e-4),
    ("wake_circulation_residual", 1e-4),
):
    maximum = max(row[f"uvlm_{quantity}"] for row in rows)
    assert maximum <= tolerance
    print(f"Peak accepted {quantity}: {maximum:.6g}")

data = Path("test/referenceData/Pazy")
mesh = "nominal" if step == 1 else "refined"
with (data / f"S0_inertia_matrices_{mesh}mesh.csv").open(newline="") as stream:
    positions = [(float(row[1]), float(row[2])) for row in list(csv.reader(stream))[1:]]
semispan = sum(math.dist(left, right) for left, right in zip(positions, positions[1:]))


def rms_error(curve, x_name, y_name, xs, ys):
    x_curve = [row[x_name] for row in curve]
    errors = []
    for x, y in zip(xs, ys):
        left = bisect_right(x_curve, x) - 1
        if left < 0 or left >= len(curve) - 1:
            continue
        fraction = (x - x_curve[left]) / (x_curve[left + 1] - x_curve[left])
        interpolated = (1 - fraction) * curve[left][y_name] + fraction * curve[left + 1][y_name]
        errors.append((interpolated - y) ** 2)
    return math.sqrt(statistics.mean(errors))


for pitch in (3, 5, 7):
    curve = [row for row in rows if row["pitch_deg"] == pitch]
    assert [row["airspeed_m_s"] for row in curve] == speeds
    endpoint = curve[-1]
    for field in ("oop_pct", "twist_deg"):
        print(f"{pitch} deg at 60 m/s {field}: {endpoint['native_' + field]:.4f} -> {endpoint['uvlm_' + field]:.4f}")
    comparisons = [("oop_pct", "dynamic_pressure_pa", f"tip_u3Vsq_rootPitch{pitch}_ExpTest1.txt", 100 / semispan)]
    if pitch != 3:
        comparisons.append(("twist_deg", "airspeed_m_s", f"tip_thetaVsU_rootPitch{pitch}_Exp.txt", 1))
    for field, x_name, filename, scale in comparisons:
        xs, ys = [list(map(float, line.split())) for line in (data / filename).read_text().splitlines() if line.strip()]
        ys = [value * scale for value in ys]
        errors = [rms_error(curve, x_name, prefix + "_" + field, xs, ys) for prefix in ("native", "uvlm")]
        print(f"{pitch} deg RMS vs experimental {field}: {errors[0]:.4f} -> {errors[1]:.4f}")
