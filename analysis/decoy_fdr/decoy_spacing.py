import numpy as np

ISO = {
 'H':[(1.0078250319,0.999885),(2.0141017779,0.000115)],
 'C':[(12.0,0.9893),(13.0033548378,0.0107)],
 'N':[(14.0030740052,0.99632),(15.0001088984,0.00368)],
 'O':[(15.9949146221,0.99757),(16.9991315,0.00038),(17.9991604,0.00205)],
 'S':[(31.97207069,0.9493),(32.9714585,0.0076),(33.96786683,0.0429),(35.96708088,0.0002)],
 'P':[(30.97376151,1.0)],
 'Cl':[(34.96885271,0.7578),(36.9659026,0.2422)],
 'Br':[(78.9183376,0.5069),(80.9162906,0.4931)],
}

def conv(dist, elem, n, prune=1e-9, mtol=1e-4):
    single = ISO[elem]
    for _ in range(n):
        nd={}
        for m,p in dist.items():
            for em,ep in single:
                mm=m+em; pp=p*ep
                if pp<prune: continue
                key=round(mm/mtol)
                nd[key]=nd.get(key,0.0)+pp
        dist={}
        for key,pp in nd.items():
            if pp<prune: continue
            dist[key*mtol]=pp
    return dist

def pattern(formula):
    dist={0.0:1.0}
    for elem,n in formula.items():
        if n>0: dist=conv(dist,elem,n)
    return dist

def peaks(dist, base_frac=0.01):
    mono=min(dist)
    groups={}
    for m,p in dist.items():
        idx=int(round(m-mono))
        groups.setdefault(idx,[0.0,0.0])
        groups[idx][0]+=m*p
        groups[idx][1]+=p
    pk=[]
    for idx in sorted(groups):
        wm,tp=groups[idx]
        pk.append((idx, wm/tp, tp))
    mx=max(p[2] for p in pk)
    pk=[(i,c,t/mx) for i,c,t in pk if t/mx>=base_frac]
    return mono,pk

def eff_spacing(pk):
    idx=np.array([p[0] for p in pk],float)
    cen=np.array([p[1] for p in pk],float)
    w  =np.array([p[2] for p in pk],float)
    W=np.diag(w); X=np.vstack([np.ones_like(idx),idx]).T
    beta=np.linalg.solve(X.T@W@X, X.T@W@cen)
    slope=beta[1]
    gaps=[pk[i+1][1]-pk[i][1] for i in range(len(pk)-1) if pk[i+1][0]-pk[i][0]==1]
    naive=np.mean(gaps) if gaps else float('nan')
    return slope,naive

def averagine(mass, extra=None):
    u=mass/111.1
    f={'C':round(4.9384*u),'H':round(7.7583*u),'N':round(1.3577*u),
       'O':round(1.4773*u),'S':round(0.0417*u)}
    if extra:
        for e,n in extra.items(): f[e]=f.get(e,0)+n
    return f

C13=1.0033548378
cases=[
  ('peptide averagine ~1500', averagine(1500)),
  ('  + 6 Cl',                 averagine(1500,{'Cl':6})),
  ('  + 4 Br',                 averagine(1500,{'Br':4})),
  ('  + 10 Cl (extreme)',      averagine(1500,{'Cl':10})),
  ('  + 6 Cl + 4 Br',          averagine(1500,{'Cl':6,'Br':4})),
]
print(f"{'case':24s} {'mono':>9s} {'slopeD':>9s} {'naiveD':>9s}  dev/step   +4 tooth offset vs 1.00336 lattice")
print("-"*104)
for name,f in cases:
    mono,pk=peaks(pattern(f))
    slope,naive=eff_spacing(pk)
    dev=slope-C13
    off4_da=(mono+4*slope)-(mono+4*C13)
    mz=mono/1.0
    off4_ppm=off4_da/mz*1e6
    off4_z2_th=off4_da/2
    print(f"{name:24s} {mono:9.3f} {slope:9.5f} {naive:9.5f} {dev*1000:+7.2f} mDa  "
          f"{off4_da*1000:+6.1f}mDa={off4_ppm:+5.1f}ppm(z1)/{off4_z2_th*1000:+5.1f}mTh(z2)")
    shape=" ".join(f"{i}:{t:.2f}" for i,c,t in pk[:8])
    print(f"{'':24s} shape: {shape}")
