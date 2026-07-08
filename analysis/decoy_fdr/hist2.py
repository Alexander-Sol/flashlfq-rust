import numpy as np, csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = os.environ["CLAUDE_JOB_DIR"] + r"\tmp"

def load(path):
    sc=[];it=[];ni=[]
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sc.append(float(row["Detector Score"])); it.append(float(row["Summed Intensity"]))
            ni.append(int(row["Num Isotopes"]))
    return np.array(sc),np.array(it),np.array(ni)

sets = {
 "target (uniform)": D+r"\target.detected.tsv",
 "off-lattice x0.5": D+r"\off05.detected.tsv",
 "off-lattice x0.9": D+r"\off09.detected.tsv",
 "mixed-charge":     D+r"\mixed.detected.tsv",
}
data={k:load(v) for k,v in sets.items()}

print(f"{'set':18} {'n':>8} {'fit p50':>8} {'fit p90':>8} {'iso p50':>8} {'iso>=4 %':>9} {'iso==2 %':>9}")
for k,(sc,it,ni) in data.items():
    fit=sc/np.maximum(it,1e-9)
    print(f"{k:18} {len(sc):8,} {np.median(fit):8.3f} {np.percentile(fit,90):8.3f} "
          f"{int(np.median(ni)):8d} {100*np.mean(ni>=4):8.1f}% {100*np.mean(ni==2):8.1f}%")

# Num-isotopes distribution table
print("\nNum-isotopes distribution (% of features):")
maxk=8
print("iso  " + " ".join(f"{k[:10]:>11}" for k in data))
for j in range(2,maxk+1):
    lab = f">={j}" if j==maxk else str(j)
    row=[]
    for k,(sc,it,ni) in data.items():
        pct = 100*np.mean(ni>=j) if j==maxk else 100*np.mean(ni==j)
        row.append(f"{pct:10.1f}%")
    print(f"{lab:4} "+" ".join(row))

# figure: fit proxy + num isotopes
fig,ax=plt.subplots(1,2,figsize=(13,5))
colors={"target (uniform)":"#2f7bd6","off-lattice x0.5":"#e08a1e","off-lattice x0.9":"#d64550","mixed-charge":"#7a4fbf"}
bins=np.linspace(0,1.2,50)
for k,(sc,it,ni) in data.items():
    fit=sc/np.maximum(it,1e-9)
    ax[0].hist(fit,bins=bins,density=True,histtype="step",lw=2,color=colors[k],label=k)
ax[0].set_xlabel("fit proxy = detector score / summed intensity"); ax[0].set_ylabel("density")
ax[0].set_title("Detection fit proxy"); ax[0].legend(fontsize=9)

isob=np.arange(2,10)-0.5
for k,(sc,it,ni) in data.items():
    ax[1].hist(ni,bins=isob,density=True,histtype="step",lw=2,color=colors[k],label=k)
ax[1].set_xlabel("num isotopes observed"); ax[1].set_ylabel("fraction")
ax[1].set_title("Teeth found per feature"); ax[1].legend(fontsize=9)
fig.suptitle("Positional decoys vs target — CA/Lumos cov90, detect-only (averagine weights)")
fig.tight_layout()
out=D+r"\decoy_positional.png"; fig.savefig(out,dpi=130); print("\nwrote",out)
