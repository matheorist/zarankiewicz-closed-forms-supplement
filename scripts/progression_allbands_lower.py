# =====================================================================
#  Lower bound Es>=F_r (hence identity Es=F_r) for ALL progression bands r=5..9.
#  Same machine per band: 4x4 primal, objective F_r, nonneg threshold j_r, all cuts.
#  Prints the factored slacks (sign-definite) + exact feasibility per band.
# =====================================================================
import sympy as sp
from fractions import Fraction as F
import math
j,a,n,b,q = sp.symbols('j alpha n beta q', positive=True)
def C2(x): return sp.Rational(1,2)*x*(x-1)
def C3(x): return sp.Rational(1,6)*x*(x-1)*(x-2)
def band(r): return (-(-r**3//4), -(-(r+1)**3//4)-1)
def Fr(r,m):
    M=r*(r+1)//2
    P=(2*r+1)*m**3+sp.Rational((3*r+1)*(3*r+2),2)*m**2+3*(2*r+1)*M*m+2*M**2
    Q=(6*M+1)*m**2+3*(2*r+1)*M*m+2*M**2
    return 2*(m+2)*(n*P+M**2*m*(m+4)*(m-1))/((m+4)*Q)
def thf(k,i,m):
    kk=k-2; al=(2*(m-2))%kk
    return F(i*(i-1),2)*F(i-2-al,kk-al) if i<k else F(i*(i-1),2)

for r in [5,6,7,8,9]:
    M=r*(r+1)//2; m=M*j; lo,hi=band(r)
    S=[r*j+1,r*j+2,(r+1)*j+1,(r+1)*j+2]
    k1=(2*(m-2)-a)/r; d1=k1-a; k2=k1+1; d2=k2-(a-r)
    t1=lambda x:C2(x)*(x-2-a)/d1; t2=lambda x:C2(x)*(x-2-(a-r))/d2
    A=sp.Matrix([[sp.Integer(1)]*4,[C3(s) for s in S],
                 [t1(S[0]),t1(S[1]),C2(S[2]),C2(S[3])],[t2(S[0]),t2(S[1]),C2(S[2]),C2(S[3])]])
    W=[sp.cancel(w) for w in A.solve(sp.Matrix([n,2*C3(m),r*C2(m),r*C2(m)]))]  # alpha-free
    obj_ok=sp.simplify(sum(S[i]*W[i] for i in range(4))-Fr(r,m))==0
    U0=sp.solve(sp.numer(sp.together(W[0])),n)[0]; V0=sp.solve(sp.numer(sp.together(W[2])),n)[0]
    jr=math.ceil(max([float(x) for x in sp.real_roots(sp.Poly(sp.numer(sp.together(V0-hi)),j))]+[0]))
    # slacks
    kb=(2*(m-2)-b)/r; kmb=kb-b
    qr=sp.simplify(r*C2(m)-(C2(S[0])*(S[0]-2-b)/kmb*W[0]+C2(S[1])*(S[1]-2-b)/kmb*W[1]+C2(S[2])*W[2]+C2(S[3])*W[3]))
    al=r*j-4; d1t=sp.simplify(r*C2(m)-(C2(S[0])*(S[0]-2-al)/(r*j-al)*W[0]+C2(S[1])*W[1]+C2(S[2])*W[2]+C2(S[3])*W[3]))
    sl_hi=sp.factor((r+1)*C2(m)-sum(C2(S[i])*W[i] for i in range(4)))
    al=(r+1)*j-4; sl_d2=sp.factor((r-1)*C2(m)-(C2(S[0])*(S[0]-2-al)/((r+1)*j-al)*W[0]+C2(S[1])*(S[1]-2-al)/((r+1)*j-al)*W[1]+C2(S[2])*(S[2]-2-al)/((r+1)*j-al)*W[2]+C2(S[3])*W[3]))
    kg=(2*(m-2)-b)/q; kmbg=kg-b
    BIG=sp.cancel(sp.numer(sp.together(q*C2(m)-sum(C2(S[i])*(S[i]-2-b)/kmbg*W[i] for i in range(4))))/(b*j*q))
    big_pos=all(
        all(
            c >= 0
            for c in sp.Poly(sp.expand(BIG.subs({n: nn, q: qq})), j).all_coeffs()
        )
        for nn in range(lo, hi + 1)
        for qq in range(1, r)
    )
    # exact feasibility using alpha-free weights
    def feas(jj,nn):
        mm=M*jj; Sn=[r*jj+1,r*jj+2,(r+1)*jj+1,(r+1)*jj+2]
        wv=[W[i].subs({j:jj,n:nn,a:0}) for i in range(4)]   # a-free => any a works
        wv=[F(int(sp.numer(v)),int(sp.denom(v))) for v in wv]
        if any(v<0 for v in wv): return False
        for k in range(3,mm+1):
            kc=k-2; al2=(2*(mm-2))%kc
            if sum(thf(k,Sn[i],mm)*wv[i] for i in range(4))>F(mm*(mm-1),2)*F(2*(mm-2)-al2,kc): return False
        return True
    fe=all(feas(jj,nn) for jj in range(jr,jr+12) for nn in range(lo,hi+1))
    print(f"--- r={r}: band [{lo},{hi}],  j_r={jr} ---")
    print(f"   obj=F_r:{obj_ok}  qr-tight:{qr==0}  D_(rj+2)-tight:{d1t==0}  q<=r-1 BIG>0:{big_pos}  exact-feasible(j_r..j_r+11):{fe}")
    print(f"   q>=r+1 slack = {sl_hi}")
    print(f"   D_((r+1)j+2) slack = {sl_d2}")
