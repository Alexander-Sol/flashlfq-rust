import numpy as np, csv, os
D = os.environ["CLAUDE_JOB_DIR"] + r"\tmp"
FEATS = ["cosine","n_teeth","contig_run","ppm_spread","coelution","rt_gauss"]

def load(p):
    rows=[]
    with open(p,newline="") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            d={}
            for m in FEATS:
                v=r[m]; d[m]=np.nan if v=="NA" else float(v)
            rows.append(d)
    return rows

def matrix(rows):
    X=np.empty((len(rows),len(FEATS)))
    for i,r in enumerate(rows):
        for j,m in enumerate(FEATS):
            v=r[m]
            if np.isnan(v):
                v = 12.0 if m=="ppm_spread" else (0.0 if m in("cosine",) else -1.0)
            X[i,j]=v
    return X

T   = matrix(load(D+r"\mt_target.tsv"))
OFF = matrix(load(D+r"\mt_off09.tsv"))
ROT = matrix(load(D+r"\mt_rotated.tsv"))
MIX = matrix(load(D+r"\mt_mixed.tsv"))

# deterministic even/odd split (no RNG)
def split(X): return X[0::2], X[1::2]
T_tr,T_te   = split(T)
OFF_tr,OFF_te = split(OFF)

# standardize on training pool
Xtr_pool = np.vstack([T_tr, OFF_tr])
mu, sd = Xtr_pool.mean(0), Xtr_pool.std(0)+1e-9
z = lambda X:(X-mu)/sd

# balance decoy to target size for training (deterministic stride subsample)
def sub(X,n):
    if len(X)<=n: return X
    idx=np.linspace(0,len(X)-1,n).astype(int); return X[idx]
Xpos=z(T_tr); Xneg=z(sub(OFF_tr,len(T_tr)))
X=np.vstack([Xpos,Xneg]); y=np.concatenate([np.ones(len(Xpos)),np.zeros(len(Xneg))])
Xb=np.hstack([X,np.ones((len(X),1))])  # bias

# logistic regression via gradient descent (L2)
w=np.zeros(Xb.shape[1]); lr=0.5; lam=1e-3
for _ in range(2000):
    p=1/(1+np.exp(-Xb@w)); g=Xb.T@(p-y)/len(y)+lam*np.r_[w[:-1],0]; w-=lr*g
def score(X): return z(X)@w[:-1]+w[-1]

def rankdata(a):
    o=np.argsort(a,kind="mergesort"); r=np.empty(len(a)); r[o]=np.arange(1,len(a)+1)
    s=a[o]; i=0
    while i<len(a):
        j=i
        while j+1<len(a) and s[j+1]==s[i]: j+=1
        if j>i: r[o[i:j+1]]=(i+1+j+1)/2.0
        i=j+1
    return r
def auc(pos,neg):
    allv=np.concatenate([pos,neg]); r=rankdata(allv)
    return (r[:len(pos)].sum()-len(pos)*(len(pos)+1)/2)/(len(pos)*len(neg))

print("Held-out AUC (target_test vs decoy).  cosine-only vs full classifier:\n")
print(f"{'decoy':12}{'n_decoy':>9}{'cosine':>9}{'classifier':>12}{'lift':>8}")
cos_i = FEATS.index("cosine")
def cos_auc(pos,neg): return auc(pos[:,cos_i], neg[:,cos_i])
for name,Dset in [("off x0.9*",OFF_te),("rotated",ROT),("mixed",MIX)]:
    ca=cos_auc(T_te,Dset); la=auc(score(T_te),score(Dset))
    print(f"{name:12}{len(Dset):9,}{ca:9.3f}{la:12.3f}{la-ca:+8.3f}")
print("  (* off x0.9 test split is in-distribution; rotated/mixed are cross-decoy held-out)")

print("\nStandardized logistic weights (feature importance):")
for f,c in sorted(zip(FEATS,w[:-1]),key=lambda t:-abs(t[1])):
    print(f"  {f:12}{c:+7.3f}")

# --- q-values on held-out target_test vs held-out decoy (mixed), classifier vs cosine ---
def qvals(pos_s, neg_s, Nt, Nd):
    # scale decoy trials to target trial count; FDR(t)=(D>=t * Nt/Nd)/(T>=t); monotonic q
    s=Nt/Nd
    order=np.argsort(-np.concatenate([pos_s,neg_s]))
    lab=np.concatenate([np.ones(len(pos_s)),np.zeros(len(neg_s))])[order]
    T=0;Dc=0;q=[]
    fdrs=[]
    for l in lab:
        if l==1:T+=1
        else:Dc+=1
        fdrs.append((Dc*s)/max(T,1))
    # monotonize from the bottom
    fdrs=np.array(fdrs); qv=np.minimum.accumulate(fdrs[::-1])[::-1]
    # count targets passing q thresholds
    tmask=lab==1
    return qv, tmask
def targets_at_q(pos_s,neg_s,Nt,Nd,thr):
    qv,tmask=qvals(pos_s,neg_s,Nt,Nd)
    return int(((qv<=thr)&tmask).sum())

Nt,Nd=len(T),len(MIX)
print(f"\nTarget features passing q-cutoff (held-out null = mixed, scaled Nt/Nd={Nt/Nd:.2f}):")
print(f"{'cutoff':>8}{'cosine-only':>13}{'classifier':>12}")
for thr in [0.01,0.05,0.10]:
    c=targets_at_q(T_te[:,cos_i],MIX[:,cos_i],Nt,Nd,thr)
    l=targets_at_q(score(T_te),score(MIX),Nt,Nd,thr)
    print(f"{thr:8.2f}{c:13,}{l:12,}")
print(f"(target_test n={len(T_te):,}; counts are the held-out half, ~2x for full set)")
