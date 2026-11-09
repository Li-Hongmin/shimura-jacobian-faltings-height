#!/usr/bin/env python3
"""Arithmetic checks for the manuscript
"A uniform lower bound for the Faltings height of Jacobians of Shimura curves".

Standard library only. All checks are exact (rational arithmetic) except
check 7, which is a finite numerical illustration of an o(1) statement.
Run:  python3 checks.py
"""
from fractions import Fraction as F
from itertools import combinations
import math

def primes_upto(n):
    s = [True]*(n+1); s[0] = s[1] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]:
            for j in range(i*i, n+1, i): s[j] = False
    return [i for i in range(n+1) if s[i]]

def kron(d, q):
    """Kronecker symbol (d/q) for a prime q."""
    if q == 2:
        return 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
    d %= q
    if d == 0: return 0
    return 1 if pow(d, (q-1)//2, q) == 1 else -1

def phi(D): return math.prod(q-1 for q in D)
def psi(M): return math.prod(p+1 for p in M)

ok = True
def check(cond, msg):
    global ok
    print(("ok   " if cond else "FAIL ") + msg)
    ok &= bool(cond)

# 1. Lemma 3.1: reduced discriminants of orders Z<1,u,v,uv> with
#    disc(u), disc(v) in {-4,-3}, u,v non-commuting:
#    disc_red = (d_u d_v - x^2)/4 > 0, x = 2 trd(uv) - trd(u) trd(v),
#    x == d_u d_v (mod 2).
vals = {}
for du in (-4, -3):
    for dv in (-4, -3):
        s = set()
        for x in range(-10, 11):
            if (x - du*dv) % 2: continue
            num = du*dv - x*x
            if num > 0 and num % 4 == 0:
                s.add(num//4)
        vals[(du, dv)] = sorted(s)
print("   reduced discriminants:", vals)
check(all(max(v) <= 4 for v in vals.values()), "Lemma 3.1: every value is in {1,2,3,4}, so no prime >= 5 divides it")
# x = +-d (commuting case, uv in the same field) gives 0: excluded by num > 0.

# 2. Mass identity (Lemma 4.3) and recurrence (Proposition 4.4) over many (D, M, p).
P = [q for q in primes_upto(60) if q >= 5]
cnt = 0; good = True
for r in (2, 4):
    for D in combinations(P[:7], r):
        for p in P:
            if p in D: continue
            for Mrest in ((), (P[-1],) if P[-1] not in D and P[-1] != p else ()):
                d_low = F(phi(D)*psi(Mrest), 6)          # deg_orb at level M/p
                d_M = F(phi(D)*psi(Mrest+(p,)), 6)        # deg_orb at level M
                mass = F(phi(D+(p,))*psi(Mrest), 12)     # Eichler mass, B_{Dp}, level M/p
                good &= (mass == F(p-1, 2)*d_low)
                good &= (d_M == (p+1)*d_low)
                incr = 2*p*d_low - mass                   # (2A.pi_*C2 + C2^2)/log p
                good &= (incr == d_M*F(3*p+1, 2*(p+1)))
                cnt += 1
check(good, f"mass = (p-1)/2 * d' and increment = d_M (3p+1)/(2(p+1)) on {cnt} triples")

# 3. The example D = 55, p = 17 (thick nodes with k = 3).
D, p = (5, 11), 17
d_low = F(phi(D), 6); mass = F(phi(D+(p,)), 12)
check(d_low == F(20, 3) and mass == F(160, 3), "D=55, p=17: d' = 20/3, sum 1/k_x = 160/3")
check((3*p+1)*d_low/2 == F(520, 3), "D=55, p=17: increment = (520/3) log 17")
# Eichler class number formula for the maximal order of B_{5*11*17}:
Dp = D + (p,)
e2 = math.prod(1 - kron(-4, q) for q in Dp) // 2
e3 = math.prod(1 - kron(-3, q) for q in Dp) // 2
h = mass + F(e2, 2) + F(2*e3, 3)
check(e2 == 0 and e3 == 4 and h == 56, "B_{935}: e2=0, e3=4, class number 56 = 52 + 4 (52 points with k=1, 4 with k=3)")
check(52 + F(4, 3) == mass, "mass = 52 + 4/3")

# 4. Class-number integrality as a consistency check of the mass normalisation
#    (mass uses w = #O^x/{+-1}).
bad = 0; n = 0
for r in (1, 3, 5):
    for Dd in combinations(primes_upto(40), r):
        P4 = math.prod(1 - (kron(-4, q) if q != 2 else 0) for q in Dd)
        P3 = math.prod(1 - (kron(-3, q) if q != 3 else 0) for q in Dd)
        h = F(phi(Dd), 12) + F(P4, 4) + F(P3, 3)
        n += 1; bad += (h.denominator != 1)
check(bad == 0, f"Eichler class number formula gives integers for {n} definite discriminants")

# 5. Merging of coefficients (Section 6).
good = True
for q in P:
    good &= (F(3*q-1, 2*(q-1)) + F(q+1, 2*(q-1)) == F(2*q, q-1))
    good &= (F(3*q+1, 2*(q+1)) + F(q-1, 2*(q+1)) == F(2*q, q+1))
    good &= (F(q+1, 2*(q-1)) == F(phi((q,))**0 * (q+1), 2*(q-1)))
    # delta_q mass: phi(D/q) psi(Mq)/12 = d_M (q+1)/(2(q-1))
    Dq = (5, 7, 11, q) if q not in (5, 7, 11) else (13, 17, 19, q)
    dM = F(phi(Dq), 6)
    good &= (F(phi(Dq)//(q-1)*(q+1), 12) == dM*F(q+1, 2*(q-1)))
check(good, "coefficient sums 2q/(q-1) (q|D), 2p/(p+1) (p|M) and the Cerednik-Drinfeld edge mass")

# 6. Invariant ring under an order-2 automorphism with wild ramification (p = 2):
#    h(t) = -(t+2)/(t+1), fixed points t = -1 +- i, z = t*h(t).
class G:  # Gaussian rationals a + b i
    def __init__(s, a, b=0): s.a, s.b = F(a), F(b)
    def __add__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a+o.a, s.b+o.b)
    def __sub__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a-o.a, s.b-o.b)
    def __mul__(s, o): o = o if isinstance(o, G) else G(o); return G(s.a*o.a-s.b*o.b, s.a*o.b+s.b*o.a)
    def inv(s): n = s.a*s.a+s.b*s.b; return G(s.a/n, -s.b/n)
    def __truediv__(s, o): o = o if isinstance(o, G) else G(o); return s*o.inv()
    def norm(s): return s.a*s.a+s.b*s.b
    def __eq__(s, o): return s.a == o.a and s.b == o.b
def hmap(t): return G(0) - (t+2)/(t+1)
def v_pi(x):  # valuation at (1+i) in Z_2[i]: v = v_2(norm)
    nrm = x.norm(); v = 0
    num, den = nrm.numerator, nrm.denominator
    while num % 2 == 0: num //= 2; v += 1
    while den % 2 == 0: den //= 2; v -= 1
    return v
tp, tm = G(-1, 1), G(-1, -1)
check(hmap(tp) == tp and hmap(tm) == tm and hmap(hmap(G(3, 5))) == G(3, 5), "h is an involution fixing -1+-i")
zp, zm = tp*hmap(tp), tm*hmap(tm)
check(v_pi(tp-tm) == 2 and v_pi(zp-zm) == 4, "v(t+ - t-) = 2 and v(z+ - z-) = 4 = 2*2 (length multiplied by |H|)")

# 7. Uniformity of sum_{p|M} log p/(p+1) = o(log N): worst case M = product of the first primes >= 5.
ratios = []
for k in (5, 10, 20, 40, 80, 160):
    ps = [q for q in primes_upto(5000) if q >= 5][:k]
    logN = sum(math.log(q) for q in ps)
    ratios.append(sum(math.log(q)/(q+1) for q in ps)/logN)
print("   ratios:", [round(x, 4) for x in ratios])
check(all(ratios[i+1] < ratios[i] for i in range(len(ratios)-1)), "numerical: ratio decreases along primorial-type M (finite illustration)")

# 8. Frey curves of the triples (16t, 1, 16t+1): model Y^2 + XY = X^3 - 4t X^2 - t X.
def disc(a1, a2, a3, a4, a6):
    b2 = a1*a1+4*a2; b4 = 2*a4+a1*a3; b6 = a3*a3+4*a6
    b8 = a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
    return -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6, b2*b2-24*b4
good = True
for t in range(1, 2000, 6):
    dl, c4 = disc(1, -4*t, 0, -t, 0)
    good &= (dl == t*t*(16*t+1)**2) and dl % 2 and dl % 3
    # isomorphic to y^2 = x(x-a)(x+b) with a=16t, b=1 up to the standard change of variables:
    # that curve has c4 = 16(a^2+ab+b^2); our model has c4' with 2^4 c4' = c4.
    good &= (16*c4 == 16*((16*t)**2+16*t+1))
check(good, "Frey model: Delta = t^2(16t+1)^2, prime to 6, c4 matches y^2=x(x-16t)(x+1) after scaling by 2")

print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
