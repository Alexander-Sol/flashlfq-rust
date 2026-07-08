import numpy as np, csv, os
D = os.environ["CLAUDE_JOB_DIR"] + r"\tmp"
METRICS = ["n_teeth","contig_run","ppm_spread","coelution","rt_gauss","cosine"]
# direction: +1 => higher is more target-like; -1 => lower is more target-like
DIRN = {"n_teeth":1,"contig_run":1,"ppm_spread":-1,"coelution":1,"rt_gauss":1,"cosine":1}

def load(p):
    d={m:[] for m in METRICS}
    with open(p,newline="") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            for m in METRICS:
                v=row[m]
                d[m].append(np.nan if v=="NA" else float(v))
    return {m:np.array(v) for m,v in d.items()}

T = load(D+r"\mt_target.tsv")
sets = {"rotated":load(D+r"\mt_rotated.tsv"),
        "off x0.9":load(D+r"\mt_off09.tsv"),
        "mixed":load(D+r"\mt_mixed.tsv"),
        "rt-unif":load(D+r"\mt_rtuni.tsv"),
        "rt-inv":load(D+r"\mt_rtinv.tsv")}

def rankdata(a):
    order=np.argsort(a,kind="mergesort"); ranks=np.empty(len(a)); ranks[order]=np.arange(1,len(a)+1)
    # average ties
    a_sorted=a[order]; i=0
    while i<len(a):
        j=i
        while j+1<len(a) and a_sorted[j+1]==a_sorted[i]: j+=1
        if j>i:
            avg=(i+1+j+1)/2.0
            ranks[order[i:j+1]]=avg
        i=j+1
    return ranks

def auc(pos,neg,dirn):
    pos=pos[~np.isnan(pos)]; neg=neg[~np.isnan(neg)]
    if len(pos)==0 or len(neg)==0: return np.nan,len(pos),len(neg)
    if dirn<0: pos,neg=-pos,-neg
    allv=np.concatenate([pos,neg]); r=rankdata(allv)
    a=(r[:len(pos)].sum()-len(pos)*(len(pos)+1)/2.0)/(len(pos)*len(neg))
    return a,len(pos),len(neg)

print("Per-metric AUC (target vs decoy; 0.5 = no separation, 1.0 = perfect). dir shown.\n")
hdr=f"{'metric':11}{'dir':>4}"+"".join(f"{k:>12}" for k in sets)
print(hdr); print("-"*len(hdr))
for m in METRICS:
    row=f"{m:11}{'hi' if DIRN[m]>0 else 'lo':>4}"
    for k,S in sets.items():
        a,_,_=auc(T[m],S[m],DIRN[m])
        row+=f"{a:12.3f}"
    print(row)

print("\nDefined-fraction (features with the metric non-NA):")
hdr=f"{'metric':11}"+f"{'target':>10}"+"".join(f"{k:>10}" for k in sets)
print(hdr)
for m in ["ppm_spread","coelution","rt_gauss"]:
    row=f"{m:11}{np.mean(~np.isnan(T[m])):10.2f}"
    for k,S in sets.items():
        row+=f"{np.mean(~np.isnan(S[m])):10.2f}"
    print(row)

print("\nMedians (target vs decoys):")
hdr=f"{'metric':11}{'target':>10}"+"".join(f"{k:>10}" for k in sets)
print(hdr)
for m in METRICS:
    row=f"{m:11}{np.nanmedian(T[m]):10.3f}"
    for k,S in sets.items():
        row+=f"{np.nanmedian(S[m]):10.3f}"
    print(row)

# orthogonality: correlation matrix among metrics on target (defined rows)
print("\nMetric correlation matrix (target features, pairwise-complete):")
print("           "+"".join(f"{m[:8]:>9}" for m in METRICS))
M=np.column_stack([T[m] for m in METRICS])
for i,mi in enumerate(METRICS):
    row=f"{mi:11}"
    for j in range(len(METRICS)):
        mask=~np.isnan(M[:,i])&~np.isnan(M[:,j])
        if mask.sum()>10:
            c=np.corrcoef(M[mask,i],M[mask,j])[0,1]
            row+=f"{c:9.2f}"
        else: row+="      NA "
    print(row)
