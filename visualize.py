#!/usr/bin/env python3
"""Visualize top 20 strange attractors by interestingness."""

import json
import os
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

ROOT = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(ROOT, "results")
GALLERY = os.path.join(ROOT, "gallery")


def load_attractors():
    """Load all JSON files from results/."""
    attractors = []
    for path in glob.glob(os.path.join(RESULTS, "*.json")):
        with open(path) as f:
            attractors.append(json.load(f))
    return attractors


def score(a):
    """Interestingness score."""
    lam1 = a["spectrum"][0]
    return a["ky_dim"] * 2 + min(lam1, 10)


def traj_valid(a):
    """Check trajectory has no NaN/None."""
    for p in a["trajectory"]:
        for v in p:
            if v is None or v != v:
                return False
            if abs(v) > 1e6:
                return False
    return True


def render(a, out_path):
    """Render 3D scatter plot of trajectory."""
    traj = a["trajectory"]
    xs = [p[0] for p in traj]
    ys = [p[1] for p in traj]
    zs = [p[2] for p in traj]
    n = len(traj)

    fig = plt.figure(
        figsize=(8, 8),
        facecolor="black",
    )
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("black")

    sc = ax.scatter(
        xs, ys, zs,
        c=range(n),
        cmap="inferno",
        s=0.5,
        alpha=0.8,
    )

    ax.set_title(
        f"\u03bb\u2081={a['spectrum'][0]:.3f}"
        f"  dim={a['ky_dim']:.2f}"
        f"  id={a['id']:016x}",
        color="white",
        fontsize=10,
    )
    for spine in [ax.xaxis, ax.yaxis, ax.zaxis]:
        spine.label.set_color("white")
        spine.set_tick_params(colors="white")
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    fig.savefig(
        out_path,
        dpi=100,
        facecolor="black",
        bbox_inches="tight",
    )
    plt.close(fig)


def gen_html(top20):
    """Generate index.html grid."""
    rows = []
    for a in top20:
        fname = f"{a['id']:016x}.png"
        rows.append(
            f'  <div class="card">'
            f'<img src="{fname}"></div>'
        )
    html = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Strange Attractor Gallery</title>
<style>
body {{
  background: #111; color: #eee;
  font-family: monospace;
}}
.grid {{
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: 10px; padding: 10px;
}}
.card img {{
  width: 100%; border-radius: 4px;
}}
h1 {{ text-align: center; }}
</style></head><body>
<h1>Top 20 Strange Attractors</h1>
<div class="grid">
{chr(10).join(rows)}
</div></body></html>"""
    path = os.path.join(GALLERY, "index.html")
    with open(path, "w") as f:
        f.write(html)


def main():
    os.makedirs(GALLERY, exist_ok=True)

    print("Loading attractors...")
    attractors = load_attractors()
    print(f"Loaded {len(attractors)} attractors.")

    valid = [a for a in attractors if traj_valid(a)]
    print(f"Valid trajectories: {len(valid)}")

    ranked = sorted(
        valid, key=score, reverse=True
    )
    top20 = ranked[:20]

    # Summary table
    print(
        f"{'Rank':>4}  {'ID':>16}"
        f"  {'λ1':>8}  {'dim':>6}  {'score':>7}"
    )
    print("-" * 50)
    for i, a in enumerate(top20, 1):
        print(
            f"{i:4d}  {a['id']:016x}"
            f"  {a['spectrum'][0]:8.4f}"
            f"  {a['ky_dim']:6.3f}"
            f"  {score(a):7.3f}"
        )

    # Render each
    for i, a in enumerate(top20, 1):
        fname = f"{a['id']:016x}.png"
        out = os.path.join(GALLERY, fname)
        print(f"[{i:2d}/20] Rendering {fname}...")
        render(a, out)

    gen_html(top20)
    print(f"\nDone. Gallery at {GALLERY}/")


if __name__ == "__main__":
    main()
