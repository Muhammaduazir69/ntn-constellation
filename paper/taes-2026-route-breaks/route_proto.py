import numpy as np, isl_core as C, sys
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
WE = 7.2921159e-5
def gs_ecef(lat, lon):
    la, lo = np.radians(lat), np.radians(lon); return C.RE*np.array([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)])
def plusgrid(sh):
    P,S=sh['P'],sh['S']; E=[]
    for p in range(P):
        for s in range(S):
            a=p*S+s; E.append((a,p*S+(s+1)%S))
            if sh['star'] and p==P-1: continue
            q=(p+1)%P; s2=s
            if not sh['star'] and p==P-1: s2=(s+sh['F'])%S   # seam phase wrap of a delta shell (approx)
            E.append((a,q*S+s2))
    return np.array(E)
def run(sh, gsA, gsB, mask=25, dt=1.0, hours=1.6, latgate=None):
    n=sh['P']*sh['S']; E=plusgrid(sh); w=2*np.pi/np.sqrt(sh['r']**3/C.MU)
    gA,gB=gs_ecef(*gsA),gs_ecef(*gsB); prev=None; forced=opt=0; steps=int(hours*3600/dt); nchg=0
    for k in range(steps):
        t=k*dt; X=C.pos(sh,w*t); th=-WE*t   # rotate sats into ECEF
        Xe=X.copy(); Xe[:,0]=np.cos(th)*X[:,0]-np.sin(th)*X[:,1]; Xe[:,1]=np.sin(th)*X[:,0]+np.cos(th)*X[:,1]
        d=np.linalg.norm(Xe[E[:,0]]-Xe[E[:,1]],axis=1); ok=d<=5000e3
        if latgate: lat=np.degrees(np.arcsin(Xe[:,2]/sh['r'])); ok&=(abs(lat[E[:,0]])<latgate)&(abs(lat[E[:,1]])<latgate)
        rows=list(E[ok,0]); cols=list(E[ok,1]); wts=list(d[ok])
        vis={}
        for gi,g in ((n,gA),(n+1,gB)):
            v=Xe-g; rr=np.linalg.norm(v,axis=1); el=np.degrees(np.arcsin(v@g/np.linalg.norm(g)/rr))
            idx=np.where(el>=mask)[0]; vis[gi]=set(idx); rows+= [gi]*len(idx); cols+=list(idx); wts+=list(rr[idx])
        G=csr_matrix((wts,(rows,cols)),shape=(n+2,n+2)); D,Pr=dijkstra(G,directed=False,indices=n,return_predecessors=True)
        path=[]; x=n+1
        while x!=-9999 and x!=n: path.append(x); x=Pr[x]
        path=tuple(path[::-1]); edges=set(zip(path[:-1],path[1:]))
        if prev is not None and path!=prev[0]:
            nchg+=1
            # was old path still valid?
            okset=set(zip(rows,cols)); okset|={(b,a) for a,b in okset}
            pe=list(zip((n,)+prev[0][:-1],prev[0]))
            if all(e in okset for e in pe): opt+=1
            else: forced+=1
        prev=(path,D[n+1])
    print(f"changes {nchg} forced {forced} optional {opt} per hour {nchg/hours:.0f}")
if __name__=='__main__':
    sh=C.walker(72,22,39,53,550e3)
    run(sh,(51.5,-0.13),(40.7,-74.0),hours=float(sys.argv[1]) if len(sys.argv)>1 else 0.5)
