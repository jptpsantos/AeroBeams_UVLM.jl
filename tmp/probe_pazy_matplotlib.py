import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from pathlib import Path

print("Matplotlib", matplotlib.__version__)
print("Serif font:", fm.findfont("Times New Roman"))
with matplotlib.rc_context({"font.family": "Times New Roman", "mathtext.fontset": "cm", "text.usetex": False}):
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.plot([30, 60], [15, 35], color="#D95319", linestyle="--", linewidth=0.55)
    ax.set_xlabel(r"Airspeed, $U_\infty$ [m/s]")
    ax.set_ylabel(r"Tip OOP displacement, $u_3$ [% semispan]")
    fig.savefig(Path(__file__).with_suffix(".png"), dpi=150)
