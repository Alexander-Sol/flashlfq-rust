import numpy as np, csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.environ["CLAUDE_JOB_DIR"] + r"\tmp"

def load(path):
    score=[]; inten=[]; z=[]; niso=[]
    with open(path, newline="") as f:
        r=csv.DictReader(f, delimiter="\t")
        for row in r:
            score.append(float(row["Detector Score"]))
            inten.append(float(row["Summed Intensity"]))
            z.append(int(row["Charge"]))
            niso.append(int(row["Num Isotopes"]))
    return (np.array(score), np.array(inten), np.array(z), np.array(niso))

ts,ti,tz,tn = load(D+r"\target.detected.tsv")
ds,di,dz,dn = load(D+r"\decoy.detected.tsv")

print(f"target n={len(ts):,}   decoy n={len(ds):,}")
def q(a,name):
    qs=np.percentile(a,[10,25,50,75,90,99])
    print(f"{name:10s} p10={qs[0]:.3g} p25={qs[1]:.3g} p50={qs[2]:.3g} p75={qs[3]:.3g} p90={qs[4]:.3g} p99={qs[5]:.3g}")

print("\n-- raw Detector Score --")
q(ts,"target"); q(ds,"decoy")

# shape/RT fit proxy: score / summed intensity  (~ intensity-weighted mean of w*g, in (0,1])
tf = ts/np.maximum(ti,1e-9)
df = ds/np.maximum(di,1e-9)
print("\n-- fit proxy = score / summed intensity --")
q(tf,"target"); q(df,"decoy")

# target-decoy style: features >= threshold on the fit proxy (scale-free, unlike raw score)
print("\n-- count >= fit-proxy threshold (raw counts; decoy pass is coverage-capped) --")
print(f"{'thr':>6} {'target':>9} {'decoy':>9} {'decoy/target':>13}")
for t in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]:
    nt=int((tf>=t).sum()); nd=int((df>=t).sum())
    print(f"{t:6.2f} {nt:9,} {nd:9,} {nd/max(nt,1):13.3f}")

# ---- figure ----
fig,ax=plt.subplots(1,2,figsize=(13,5))

# panel 1: raw score, log-x density
lo=min(ts[ts>0].min(), ds[ds>0].min()); hi=max(ts.max(), ds.max())
bins=np.logspace(np.log10(lo), np.log10(hi), 60)
ax[0].hist(ts, bins=bins, density=True, alpha=.55, label=f"target (n={len(ts):,})", color="#2f7bd6")
ax[0].hist(ds, bins=bins, density=True, alpha=.55, label=f"decoy (n={len(ds):,})", color="#d64550")
ax[0].set_xscale("log"); ax[0].set_xlabel("Detector Score (matched-filter Σ w·g·I)")
ax[0].set_ylabel("density"); ax[0].set_title("Raw detector score (scales with intensity)")
ax[0].legend()

# panel 2: fit proxy, linear density
bins2=np.linspace(0, max(tf.max(),df.max())*1.02, 60)
ax[1].hist(tf, bins=bins2, density=True, alpha=.55, label="target", color="#2f7bd6")
ax[1].hist(df, bins=bins2, density=True, alpha=.55, label="decoy", color="#d64550")
ax[1].set_xlabel("fit proxy = score / summed intensity")
ax[1].set_ylabel("density"); ax[1].set_title("Shape/RT fit quality (intensity-normalized)")
ax[1].legend()

fig.suptitle("Target (averagine) vs Decoy (chlorinated averagine) — CA/Lumos, cov90, detect-only", fontsize=12)
fig.tight_layout()
out=D+r"\decoy_hist.png"
fig.savefig(out, dpi=130)
print("\nwrote", out)
