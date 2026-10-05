import sympy as sp
from fractions import Fraction as F

c,b,q,r = sp.symbols('c b q r')
a=(q*c+r)/2
def C2(x): return x*(x-1)/2
def C3(x): return x*(x-1)*(x-2)/6
al_c=2*q+r-4; al_cm1=3*q+r-4
x2,x1,y1,y2=sp.symbols('x2 x1 y1 y2'); sizes=[c-2,c-1,c+1,c+2]; w=[x2,x1,y1,y2]
eqs=[sum(w)-b, sum(C3(sizes[j])*w[j] for j in range(4))-2*C3(a)]
eqs.append(C2(c-2)*((c-2)-2-al_cm1)/((c-3)-al_cm1)*x2 + C2(c-1)*x1+C2(c+1)*y1+C2(c+2)*y2 - C2(a)*(2*a-4-al_cm1)/(c-3))
eqs.append(C2(c-2)*((c-2)-2-al_c)/((c-2)-al_c)*x2 + C2(c-1)*((c-1)-2-al_c)/((c-2)-al_c)*x1 + C2(c+1)*y1+C2(c+2)*y2 - C2(a)*(2*a-4-al_c)/(c-2))
sol=sp.solve(eqs,[x2,x1,y1,y2],dict=True)[0]
wsol=[sol[w[j]] for j in range(4)]
Fexpr=sp.together(sum(sizes[j]*wsol[j] for j in range(4)))
ye,ys,ycm1,yc=sp.symbols('ye ys ycm1 yc')
def Srow(i,lt_cm1,lt_c):
    dcm1=(C2(i)*((i-2-al_cm1)/((c-3)-al_cm1)) if lt_cm1 else C2(i))
    dc  =(C2(i)*((i-2-al_c)/((c-2)-al_c))   if lt_c   else C2(i))
    return ye+ys*C3(i)+ycm1*dcm1+yc*dc
te=[Srow(c-2,1,1)-(c-2),Srow(c-1,0,1)-(c-1),Srow(c+1,0,0)-(c+1),Srow(c+2,0,0)-(c+2)]
dsol=sp.solve(te,[ye,ys,ycm1,yc],dict=True)[0]
# band endpoints (two distinct b*): from n_{c-2} pair and n_{c+1} pair
bA=sp.simplify(sp.solve(sp.Eq(wsol[0],0),b)[0])   # n_{c-2}=0
bB=sp.simplify(sp.solve(sp.Eq(wsol[2],0),b)[0])   # n_{c+1}=0

def C2i(i): return F(i*(i-1),2)
def C3i(i): return F(i*(i-1)*(i-2),6)
def es(A,B):
    isz=list(range(2,A+1)); nv=len(isz)
    rows=[([F(1)]*nv,F(B)),([C3i(i) for i in isz],2*C3i(A))]
    for k in range(3,A+1):
        kk=k-2; al=(2*(A-2))%kk
        rows.append(([(C2i(i)*F(i-2-al,kk-al) if i<k else C2i(i)) for i in isz],C2i(A)*F(2*(A-2)-al,kk)))
    m=len(rows); o=[F(i) for i in isz]; ncol=nv+m
    T=[rows[rr][0][:]+[F(1) if j==rr else F(0) for j in range(m)]+[rows[rr][1]] for rr in range(m)]
    cost=[o[j] for j in range(nv)]+[F(0)]*m; basis=[nv+rr for rr in range(m)]
    while True:
        e=next((j for j in range(ncol) if cost[j]>0),-1)
        if e==-1: break
        leave=-1; best=None
        for rr in range(m):
            if T[rr][e]>0:
                rt=T[rr][ncol]/T[rr][e]
                if best is None or rt<best or (rt==best and basis[rr]<basis[leave]): best=rt; leave=rr
        pv=T[leave][e]; T[leave]=[x/pv for x in T[leave]]
        for rr in range(m):
            if rr!=leave and T[rr][e]!=0:
                fc=T[rr][e]; T[rr]=[T[rr][cc]-fc*T[leave][cc] for cc in range(ncol+1)]
        if cost[e]!=0:
            fc=cost[e]; cost=[cost[cc]-fc*T[leave][cc] for cc in range(ncol)]
        basis[leave]=e
    supp={isz[basis[rr]]:T[rr][ncol] for rr in range(m) if basis[rr]<nv and T[rr][ncol]!=0}
    val=sum(o[basis[rr]]*T[rr][ncol] for rr in range(m) if basis[rr]<nv)
    return val,supp

print("Full certificate check over a grid (exact arithmetic):")
allgood=True
for (qv,cv,rv) in [(5,18,2),(5,20,2),(5,22,2),(5,19,1),(5,21,3),(6,20,0),(6,22,2),(7,23,1),(7,25,3),(8,28,2)]:
    if (qv*cv+rv)%2: continue
    av=(qv*cv+rv)//2
    sub={q:qv,c:cv,r:rv}
    Bl=bA.subs(sub); Bh=bB.subs(sub)
    lo=int(sp.ceiling(sp.Min(Bl,Bh))); hi=int(sp.floor(sp.Max(Bl,Bh)))
    # y>=0 (rational, at an interior b)
    bm=(lo+hi)//2; subb={q:qv,c:cv,r:rv,b:bm}
    yvals=[sp.Rational(dsol[s].subs(subb)) for s in (ye,ys,ycm1,yc)]
    ynn=all(v>=0 for v in yvals)
    # Es=F across the whole band. The solver's returned support is diagnostic:
    # at degenerate optima it may choose a different optimal basis.
    band_ok=True; shape_ok=True
    for bb in range(lo,hi+1):
        val,supp=es(av,bb)
        Fv=sp.Rational(Fexpr.subs({q:qv,c:cv,r:rv,b:bb}))
        if F(Fv.p,Fv.q)!=val: band_ok=False
        if sorted(k-cv for k in supp)!=[-2,-1,1,2] and lo<bb<hi: shape_ok=False
    # dual feasibility S(i)>=i for ALL i (exact), find window (binding non-support sizes)
    feas=True; window=[]
    for i in range(2,av+1):
        dcm1=(C2i(i)*F(i-2-(3*qv+rv-4),(cv-3)-(3*qv+rv-4)) if i<cv-1 else C2i(i))
        dc  =(C2i(i)*F(i-2-(2*qv+rv-4),(cv-2)-(2*qv+rv-4)) if i<cv   else C2i(i))
        Sv=sp.Rational(dsol[ye].subs(subb))+sp.Rational(dsol[ys].subs(subb))*C3i(i)+sp.Rational(dsol[ycm1].subs(subb))*sp.Rational(dcm1)+sp.Rational(dsol[yc].subs(subb))*sp.Rational(dc)-i
        if Sv<0: feas=False
        if 0<=Sv<=1 and i not in (cv-2,cv-1,cv+1,cv+2): window.append(i)
    ok = ynn and band_ok and feas
    allgood&=ok
    print(f"  q={qv} c={cv} rho={rv} a={av}: band=[{lo},{hi}] EsEqF={band_ok} solverShape={shape_ok} y>=0:{ynn} dualFeas:{feas} window={window}  {'OK' if ok else 'FAIL'}")

# The manuscript's sub-threshold example: the fixed-basis value is not optimal.
subthreshold_ok=True
for bb in range(28,44):
    val,_=es(30,bb)
    Fv=sp.Rational(Fexpr.subs({q:5,c:12,r:0,b:bb}))
    if val==F(Fv.p,Fv.q):
        subthreshold_ok=False
allgood&=subthreshold_ok
print(f"  sub-threshold q=5 a=30 b=28..43: Es!=F_q at every point: {subthreshold_ok}")

print("\nALL CERTIFICATE CHECKS PASS" if allgood else "\nSOME FAILED")
print("\nBand endpoints (exact):")
print("  b=0 for n_{c-2}:  bA =", bA)
print("  b=0 for n_{c+1}:  bB =", bB)

if not allgood:
    raise SystemExit("companion certificate verification failed")
