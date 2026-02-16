#!/usr/bin/env python3
"""Visualize top strange attractors, split by method."""

import json
import os
import glob
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
GALLERY = os.path.join(ROOT, "gallery")

TOP_N = 15


def load_metadata():
    """Load just metadata (no trajectory) for ranking."""
    out = []
    for path in glob.glob(
        os.path.join(RESULTS, "*.json")
    ):
        with open(path) as f:
            data = json.load(f)
        # Strip trajectory to save memory.
        entry = {
            "id": data["id"],
            "spectrum": data["spectrum"],
            "ky_dim": data["ky_dim"],
            "method": data.get("method", ""),
            "_path": path,
        }
        out.append(entry)
    return out


def load_full(entry):
    """Load full entry including trajectory."""
    with open(entry["_path"]) as f:
        return json.load(f)


def score(a):
    """Favor high dimension, moderate lambda."""
    lam = a["spectrum"][0]
    dim = a["ky_dim"]
    # Penalize extreme lambda — sweet spot 0.05-1.0
    if lam < 0.01:
        return -1
    lam_score = math.exp(-((math.log(lam) - math.log(0.3)) ** 2) / 2.0)
    return dim * 3.0 + lam_score * 2.0


def meta_valid(a):
    """Quick validity check on metadata."""
    s = a["spectrum"]
    return (
        all(v == v for v in s)
        and a["ky_dim"] == a["ky_dim"]
        and a["spectrum"][0] < 100
    )


def render(a, out_path, label=""):
    traj = a["trajectory"]
    xs = [p[0] for p in traj]
    ys = [p[1] for p in traj]
    zs = [p[2] for p in traj]
    n = len(traj)

    fig = plt.figure(
        figsize=(8, 8), facecolor="black",
    )
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("black")

    ax.scatter(
        xs, ys, zs,
        c=range(n), cmap="inferno",
        s=0.1, alpha=0.6,
    )

    method = a.get("method", "?")
    ax.set_title(
        f"[{method}] "
        f"\u03bb\u2081={a['spectrum'][0]:.3f}"
        f"  dim={a['ky_dim']:.2f}"
        f"  id={a['id']:016x}",
        color="white", fontsize=9,
    )
    for spine in [ax.xaxis, ax.yaxis, ax.zaxis]:
        spine.label.set_color("white")
        spine.set_tick_params(colors="white")
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    fig.savefig(
        out_path, dpi=100,
        facecolor="black",
        bbox_inches="tight",
    )
    plt.close(fig)


def gen_html(sections):
    """Generate index.html with method sections."""
    body = ""
    for title, items in sections:
        cards = ""
        for a in items:
            fname = f"{a['id']:016x}.jpg"
            cards += (
                f'<div class="card">'
                f'<img src="{fname}">'
                f'<p>\u03bb\u2081='
                f'{a["spectrum"][0]:.3f}'
                f' dim={a["ky_dim"]:.2f}</p>'
                f'</div>\n'
            )
        body += (
            f'<h2>{title}</h2>\n'
            f'<div class="grid">\n{cards}</div>\n'
        )

    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Strange Attractor Gallery</title>
<style>
body {{
  background: #000; color: #ccc;
  font-family: 'JetBrains Mono', monospace;
  padding: 20px;
}}
h1 {{ text-align: center; color: #fff; }}
h2 {{
  color: #888; border-bottom: 1px solid #333;
  padding-bottom: 8px; margin-top: 40px;
}}
.grid {{
  display: grid;
  grid-template-columns: repeat(5, 1fr);
  gap: 12px; padding: 10px 0;
}}
.card img {{
  width: 100%; border-radius: 4px;
}}
.card p {{
  font-size: 11px; color: #666;
  margin: 4px 0; text-align: center;
}}
</style></head><body>
<h1>Strange Attractor Gallery</h1>
{body}
</body></html>"""
    path = os.path.join(GALLERY, "index.html")
    with open(path, "w") as f:
        f.write(html)


def print_table(label, items):
    print(f"\n--- {label} ---")
    print(
        f"{'#':>3}  {'ID':>16}"
        f"  {'lam1':>8}  {'dim':>6}"
        f"  {'score':>7}"
    )
    for i, a in enumerate(items, 1):
        print(
            f"{i:3d}  {a['id']:016x}"
            f"  {a['spectrum'][0]:8.4f}"
            f"  {a['ky_dim']:6.3f}"
            f"  {score(a):7.3f}"
        )


def main():
    os.makedirs(GALLERY, exist_ok=True)

    print("Loading metadata...")
    all_m = load_metadata()
    print(f"Loaded {len(all_m)}")

    valid = [a for a in all_m if meta_valid(a)]
    print(f"Valid: {len(valid)}")

    # Split by method
    random = [
        a for a in valid
        if a.get("method", "") != "evolve"
    ]
    evolve = [
        a for a in valid
        if a.get("method", "") == "evolve"
    ]

    random.sort(key=score, reverse=True)
    evolve.sort(key=score, reverse=True)

    top_r = random[:TOP_N]
    top_e = evolve[:TOP_N]

    print_table("Random search", top_r)
    print_table("Evolutionary search", top_e)

    # Render one at a time to save memory.
    all_top = top_r + top_e
    for i, meta in enumerate(all_top, 1):
        fname = f"{meta['id']:016x}.jpg"
        out = os.path.join(GALLERY, fname)
        print(
            f"[{i:2d}/{len(all_top)}] "
            f"Rendering {fname}..."
        )
        full = load_full(meta)
        # Skip if trajectory has bad values.
        bad = False
        for p in full["trajectory"]:
            for v in p:
                if v is None or v != v:
                    bad = True
                    break
                if abs(v) > 1e6:
                    bad = True
                    break
            if bad:
                break
        if bad:
            print(f"  skipping (bad trajectory)")
            continue
        render(full, out)

    gen_html([
        (
            f"Random Search (top {TOP_N})",
            top_r,
        ),
        (
            f"Evolutionary Search (top {TOP_N})",
            top_e,
        ),
    ])
    print(f"\nDone. {GALLERY}/")


if __name__ == "__main__":
    main()
