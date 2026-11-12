"""Prototype: shortest-delay vs lifetime-aware ground-to-ground routing on a +Grid Walker shell."""
import numpy as np, isl_core as C, sys
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from route_proto import plusgrid, gs_ecef
WE = 7.2921159e-5

def ephem(sh, T, dt):
    w = 2*np.pi/np.sqrt(sh['r']**3/C.MU); out=[]
    for k in range(int(T/dt)):
        t=k*dt; X=C.pos(sh,w*t); th=WE*t
        Xe=X.copy(); Xe[:,0]=np.cos(th)*X[:,0]+np.sin(th)*X[:,1]; Xe[:,1]=-np.sin(th)*X[:,0]+np.cos(th)*X[:,1]
        out.append(Xe)
    return np.array(out)

def gsl(Xs, g, mask):
    v=Xs-g; rr=np.linalg.norm(v,axis=2); el=np.degrees(np.arcsin((v@g)/np.linalg.norm(g)/rr))
    return el>=mask, rr

def residual(vis, dt):
    R=np.zeros(vis.shape); run=np.zeros(vis.shape[1])
    for k in range(vis.shape[0]-1,-1,-1):
        run=np.where(vis[k],run+dt,0.0); R[k]=run
    return R

def route(n, E, dI, visA, rA, visB, rB, prune=None):
    m = np.ones(len(E),bool)
    a=np.where(visA)[0]; b=np.where(visB)[0]
    if prune is not None:
        LA,LB,L=prune; a=a[LA[a]>=L]; b=b[LB[b]>=L]
    rows=np.r_[E[:,0],[n]*len(a),[n+1]*len(b)]; cols=np.r_[E[:,1],a,b]; w=np.r_[dI,rA[a],rB[b]]
    G=csr_matrix((w,(rows,cols)),shape=(n+2,n+2)); D,Pr=dijkstra(G,directed=False,indices=n,return_predecessors=True)
    if not np.isfinite(D[n+1]): return np.inf,None
    path=[]; x=n+1
    while x!=n: path.append(int(x)); x=Pr[x]
    return D[n+1], tuple(path[::-1])

def pdelay(path, Xk, gA, gB, n):
    pts=[gA]+[Xk[s] for s in path[:-1]]+[gB]
    return sum(np.linalg.norm(pts[j+1]-pts[j]) for j in range(len(pts)-1))

def sim(sh, A, B, hours=1.0, dt=1.0, mask=25, eps=0.05, look=900, policy='life'):
    n=sh['P']*sh['S']; E=plusgrid(sh); T=hours*3600
    X=ephem(sh,T+look,dt); gA,gB=gs_ecef(*A),gs_ecef(*B)
    vA,rA=gsl(X,gA,mask); vB,rB=gsl(X,gB,mask); LA,LB=residual(vA,dt),residual(vB,dt)
    path=None; exp=0; chg=forced=comps=0; st=[]
    for k in range(int(T/dt)):
        Xk=X[k]; dI=np.linalg.norm(Xk[E[:,0]]-Xk[E[:,1]],axis=1)
        dstar,pstar=route(n,E,dI,vA[k],rA[k],vB[k],rB[k])
        valid = path is not None and vA[k][path[0]] and vB[k][path[-2]]
        need = (not valid) if policy in ('life','sticky') else True
        if policy=='life' and valid and k*dt>=exp: need=True
        if need:
            comps+=1
            if policy=='life':
                cand=np.unique(np.r_[LA[k][vA[k]],LB[k][vB[k]]]); lo,hi=0,len(cand)-1; best=pstar
                while lo<=hi:
                    mid=(lo+hi)//2; d,p=route(n,E,dI,vA[k],rA[k],vB[k],rB[k],prune=(LA[k],LB[k],cand[mid]))
                    if p is not None and d<=(1+eps)*dstar: best=p; lo=mid+1
                    else: hi=mid-1
                newp=best; exp=k*dt+min(LA[k][newp[0]],LB[k][newp[-2]])
            else: newp=pstar
            if path is not None and newp!=path:
                chg+=1; forced+= (not valid)
            path=newp
        st.append(pdelay(path,Xk,gA,gB,n)/dstar-1)
    print(f"{policy:6s} eps={eps:.2f}: changes/h {chg/hours:6.0f}  forced/h {forced/hours:5.0f}  computations/h {comps/hours:6.0f}  mean stretch {100*np.mean(st):.2f}%  p99 {100*np.percentile(st,99):.2f}%")

if __name__=='__main__':
    sh=C.walker(72,22,39,53,550e3); h=float(sys.argv[1]) if len(sys.argv)>1 else 0.5
    for pol,e in (('short',0),('sticky',0),('life',0.02),('life',0.05),('life',0.10)):
        sim(sh,(51.5,-0.13),(40.7,-74.0),hours=h,policy=pol,eps=e)
