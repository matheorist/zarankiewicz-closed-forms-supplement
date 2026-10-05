import sympy as sp
from fractions import Fraction as F
import math
j, a, i, n = sp.symbols('j alpha i n', positive=True)
def C2(x): return sp.Rational(1,2)*x*(x-1)
def C3(x): return sp.Rational(1,6)*x*(x-1)*(x-2)

def pos_large_j(expr, j0=15):
    """True if rational expr in j is > 0 for all integer j >= j0 (exact real_roots)."""
    e=sp.together(expr)
    def polypos(P):
        P=sp.Poly(sp.expand(P), j)
        if P.degree()<0 or P.LC()<=0: return False
        try:
            return all(rt < j0 for rt in P.real_roots())
        except Exception:
            return False
    return polypos(sp.numer(e)) and polypos(sp.denom(e))

def positive_coeff_rational(expr, gens):
    """Certify positivity on positive variables by positive coefficients."""
    num, den = sp.fraction(sp.cancel(expr))
    return (
        all(c > 0 for c in sp.Poly(sp.expand(num), *gens).coeffs())
        and all(c > 0 for c in sp.Poly(sp.expand(den), *gens).coeffs())
    )

def affine_positive_on_window(expr, left, right):
    """Certify an affine alpha-factor at both ends of [left,right]."""
    expr = sp.cancel(expr)
    if sp.degree(sp.numer(expr), a) > 1 or sp.degree(sp.denom(expr), a) > 0:
        return False
    return pos_large_j(expr.subs(a, left)) and pos_large_j(expr.subs(a, right))

def factored_positive_on_window(expr, left, right):
    """Certify a factored rational expression throughout an alpha-window."""
    num, den = sp.fraction(sp.factor(sp.cancel(expr)))
    coeff, factors = sp.factor_list(num)
    if not pos_large_j(coeff):
        return False
    for factor, _power in factors:
        if factor.has(a):
            if not affine_positive_on_window(factor, left, right):
                return False
        elif not pos_large_j(factor):
            return False
    return positive_coeff_rational(den, [j])

print("Symbolic upper-bound proof, all bands:")
symbolic_ok = True
for r in [5,6,7,8,9]:
    M=r*(r+1)//2; m=M*j
    k1m2=(2*(m-2)-a)/r; k2m2=k1m2+1
    den1=k1m2-a; den2=k2m2-(a-r)
    th1=lambda x: C2(x)*(x-2-a)/den1
    th2=lambda x: C2(x)*(x-2-(a-r))/den2
    sb=[r*j+1,r*j+2]; sa=[(r+1)*j+1,(r+1)*j+2]
    y0,ys,y1,y2=sp.symbols('y0 ys y1 y2')
    eqs=[sp.Eq(y0+ys*C3(s)+y1*th1(s)+y2*th2(s),s) for s in sb]+[sp.Eq(y0+ys*C3(s)+y1*C2(s)+y2*C2(s),s) for s in sa]
    sol=sp.solve(eqs,[y0,ys,y1,y2],dict=True)[0]
    Y0,Ys,Y1,Y2=sol[y0],sol[ys],sol[y1],sol[y2]
    # Multiplier window: y1 vanishes at B1 and y2 at A2.
    B1=[sp.simplify(x) for x in sp.solve(sp.numer(sp.together(Y1)),a) if sp.degree(sp.numer(sp.together(sp.simplify(x))),j)>=1][0]
    A2=[sp.simplify(x) for x in sp.solve(sp.numer(sp.together(Y2)),a) if sp.degree(sp.numer(sp.together(sp.simplify(x))),j)>=1][0]
    window_ok = sp.simplify(B1-A2-r) == 0
    # The adjacent cut keeps quotient r only when its remainder a-r is nonnegative.
    remainder_ok = pos_large_j(A2-r) and pos_large_j(2*(m-2)-(r+1)*B1)
    base_multipliers_ok = positive_coeff_rational(Y0,[j]) and positive_coeff_rational(Ys,[j])
    y1_window = sp.cancel(Y1/(B1-a))
    y2_window = sp.cancel(Y2/(a-A2))
    cut_multipliers_ok = (
        affine_positive_on_window(y1_window,A2,B1)
        and affine_positive_on_window(y2_window,A2,B1)
    )

    # The support-root quadratics are nonnegative at every integer edge size;
    # certify that the remaining factors have positive coefficients.
    g_hi=sp.cancel(Y0+Ys*C3(i)+(Y1+Y2)*C2(i)-i)
    g_lo=sp.cancel(Y0+Ys*C3(i)+Y1*th1(i)+Y2*th2(i)-i)
    hi_residual=sp.cancel(g_hi/((i-sa[0])*(i-sa[1])))
    lo_residual=sp.cancel(g_lo/((i-sb[0])*(i-sb[1])))
    residuals_ok = (
        positive_coeff_rational(hi_residual,[i,j])
        and positive_coeff_rational(lo_residual,[i,j])
    )

    # The only edge size not covered by the two generic slack formulas is i=k.
    k1=k1m2+2
    g_k=sp.factor(sp.simplify(Y0+Ys*C3(k1)+Y1*C2(k1)+Y2*th2(k1)-k1))
    gk_ok=factored_positive_on_window(g_k,A2,B1)

    Mr=M
    Pr=(2*r+1)*m**3+sp.Rational((3*r+1)*(3*r+2),2)*m**2+3*(2*r+1)*Mr*m+2*Mr**2
    Qr=(6*Mr+1)*m**2+3*(2*r+1)*Mr*m+2*Mr**2
    Fr=2*(m+2)*(n*Pr+Mr**2*m*(m+4)*(m-1))/((m+4)*Qr)
    dual=Y0*n+Ys*2*C3(m)+(Y1+Y2)*r*C2(m)
    objective_ok=sp.simplify(dual-Fr)==0

    band_ok=all((window_ok,remainder_ok,base_multipliers_ok,
                 cut_multipliers_ok,residuals_ok,gk_ok,objective_ok))
    symbolic_ok &= band_ok
    print(
        f"  r={r}: window={window_ok} remainder={remainder_ok} "
        f"base-y={base_multipliers_ok} cut-y={cut_multipliers_ok} "
        f"generic-slacks={residuals_ok} g_k={gk_ok} objective={objective_ok}"
    )

# ---- exact numeric feasibility for all bands, j in [5,120] ----
def Cf2(x): return F(x*(x-1),2)
def Cf3(x): return F(x*(x-1)*(x-2),6)
def theta(k,x,mm):
    kk=k-2; al=(2*(mm-2))%kk
    return Cf2(x)*F(x-2-al,kk-al) if x<k else Cf2(x)
def solve4(A,b):
    Mx=[A[r][:]+[b[r]] for r in range(4)]
    for c in range(4):
        p=next((r for r in range(c,4) if Mx[r][c]!=0),None)
        if p is None: return None
        Mx[c],Mx[p]=Mx[p],Mx[c]; pv=Mx[c][c]; Mx[c]=[x/pv for x in Mx[c]]
        for r in range(4):
            if r!=c and Mx[r][c]!=0:
                f=Mx[r][c]; Mx[r]=[Mx[r][t]-f*Mx[c][t] for t in range(5)]
    return [Mx[r][4] for r in range(4)]
print("\nExact feasibility (consecutive cut found) for j in [5,120]:")
finite_ok = True
for r in [5,6,7,8,9]:
    M=r*(r+1)//2; bad=[]
    for jj in range(5,121):
        mm=M*jj
        found=False
        for k in range(r*jj, (r+1)*jj+2):
            S=[r*jj+1,r*jj+2,(r+1)*jj+1,(r+1)*jj+2]
            Am=[[F(1),Cf3(s),theta(k,s,mm),theta(k+1,s,mm)] for s in S]
            y=solve4(Am,[F(s) for s in S])
            if y is None: continue
            if y[1]>=0 and y[2]>=0 and y[3]>=0 and min(y[0]+y[1]*Cf3(x)+y[2]*theta(k,x,mm)+y[3]*theta(k+1,x,mm)-x for x in range(2,mm+1))>=0:
                found=True; break
        if not found: bad.append(jj)
    finite_ok &= not bad
    print(f"  r={r}: {'ALL j feasible' if not bad else f'FAIL at {bad}'}")

if not (symbolic_ok and finite_ok):
    raise SystemExit("certificate verification failed")
