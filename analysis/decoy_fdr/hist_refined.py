import numpy as np, csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.environ["CLAUDE_JOB_DIR"] + r"\tmp"

def load_refined(path):
    s=[]
    with open(path, newline="") as f:
        r=csv.DictReader(f, delimiter="\t")
        for row in r:
            try: s.append(float(row["Decon Score"]))
            except: pass
    return np.array(s)

t = load_refined(D+r"\target_full.refined.tsv")
d = load_refined(D+r"\decoy_full.refined.tsv")
print(f"target refined n={len(t):,}   decoy refined n={len(d):,}")

def q(a,name):
    p=np.percentile(a,[10,25,50,75,90,99])
    print(f"{name:8s} p10={p[0]:.3f} p25={p[1]:.3f} p50={p[2]:.3f} p75={p[3]:.3f} p90={p[4]:.3f} p99={p[5]:.3f}")
print("\n-- Decon Score (envelope-fit cosine) --")
q(t,"target"); q(d,"decoy")

print("\n-- target-decoy at Decon-Score threshold --")
print(f"{'thr':>6} {'target>=':>9} {'decoy>=':>9} {'ratio d/t':>10}")
for thr in [0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.99]:
    nt=int((t>=thr).sum()); nd=int((d>=thr).sum())
    print(f"{thr:6.2f} {nt:9,} {nd:9,} {nd/max(nt,1):10.3f}")

fig,ax=plt.subplots(figsize=(8,5))
bins=np.linspace(min(t.min(),d.min()), max(t.max(),d.max()), 60)
ax.hist(t, bins=bins, density=True, alpha=.55, color="#2f7bd6", label=f"target (n={len(t):,})")
ax.hist(d, bins=bins, density=True, alpha=.55, color="#d64550", label=f"decoy (n={len(d):,})")
ax.set_xlabel("Decon Score = envelope-fit cosine (union grid)")
ax.set_ylabel("density")
ax.set_title("Refined envelope-fit: target (averagine) vs decoy (chlorinated)\nCA/Lumos cov90")
ax.legend()
fig.tight_layout()
out=D+r"\decoy_hist_refined.png"; fig.savefig(out,dpi=130); print("\nwrote",out)
