"""Task 5: Lemma 4.  V_G(r^2) = -(1/(pi r)) int_0^inf w^{-1/2} ln(1+r^2/w) G(w) dw  by mpmath quadrature."""
import mpmath as mp
mp.mp.dps = 30
def V(G, r):
    f = lambda w: w**(-mp.mpf(1)/2) * mp.log1p(r*r/w) * G(w)
    pts = [0, r*r/100, r*r, 100*r*r, 1, 10, 1000, mp.inf]
    pts = sorted(set(pts))
    return -mp.quad(f, pts) / (mp.pi * r)
def a_of(G):
    return -mp.quad(lambda w: w**(-mp.mpf(3)/2)*G(w), [0, 1e-6, 1e-3, 1, 10, 1000, mp.inf]) / mp.pi
print("constant: (1/pi) int_0^inf x^{-5/2}(x - ln(1+x)) dx =", mp.quad(lambda x: x**(-mp.mpf(5)/2)*(x-mp.log1p(x)), [0,1,mp.inf])/mp.pi, " (claim 2/3)")
s2 = mp.sqrt(2)
profiles = {
 "G_c = w/(1+w)^2 - 2sqrt2 w/(2+w)^2 (TB-4)": (lambda w: w/(1+w)**2 - 2*s2*w/(2+w)**2, 1 - 1/s2),
 "G = w e^-w - sqrt2 w e^-2w": (lambda w: w*mp.exp(-w) - s2*w*mp.exp(-2*w), 1 - s2),
 "G = w/(1+w)^3 - (8/3)... a=0 combo": None,
}
# third: w/(1+w)^3 has a = -(1/pi) B(1/2,5/2) = -(3/8); w/(b+w)^3 scales as b^{-5/2}; combine b=1 and b=4: coefficient 32
profiles["G = w/(1+w)^3 - 32 w/(4+w)^3"] = (lambda w: w/(1+w)**3 - 32*w/(4+w)**3, 1 - mp.mpf(32)/64)
del profiles["G = w/(1+w)^3 - (8/3)... a=0 combo"]
for name,(G,c1) in profiles.items():
    print(f"\n{name}:  a(G) = {mp.nstr(a_of(G),5)},  c1 = {mp.nstr(c1,8)},  (2/3)c1 = {mp.nstr(2*c1/3,10)}")
    for r in ("1e-1","1e-2","1e-3","1e-4","1e-5"):
        r=mp.mpf(r); v=V(G,r)
        print(f"   r={mp.nstr(r,2)}: V/r^2 = {mp.nstr(v/r**2,10)}   [V/r^2-(2/3)c1]/(r|ln r|) = {mp.nstr((v/r**2-2*c1/3)/(r*abs(mp.log(r))),6)}   /(r) = {mp.nstr((v/r**2-2*c1/3)/r,6)}")
print("\nGeneral G in X^0 with non-analytic edge behaviour G ~ w^gamma (5/8<gamma<1): claimed |V| = O(r^{5/4}), d_r V = O(r^{1/4})")
for gam in (mp.mpf("0.65"), mp.mpf("0.75"), mp.mpf("0.9")):
    # G = w^gam/(1+w)^2 - c w^gam/(2+w)^2 with a(G)=0 :  a(w^gam/(b+w)^2) ~ b^{gam-5/2}
    c = mp.mpf(2)**(mp.mpf(5)/2-gam)
    G = lambda w: w**gam/(1+w)**2 - c*w**gam/(2+w)**2
    print(f" gamma={mp.nstr(gam,3)}: a(G)={mp.nstr(a_of(G),3)}")
    prev=None
    for r in ("1e-1","1e-2","1e-3","1e-4"):
        r=mp.mpf(r); v=V(G,r)
        slope = "" if prev is None else f" local exponent = {mp.nstr(mp.log(abs(v/prev[1]))/mp.log(r/prev[0]),5)} (expect 2 gamma = {mp.nstr(2*gam,3)})"
        print(f"   r={mp.nstr(r,2)}: V = {mp.nstr(v,8)}  V/r^(5/4) = {mp.nstr(v/r**(mp.mpf(5)/4),6)}{slope}")
        prev=(r,v)
print("\nControl with a(G) != 0: G = w/(2+w)^2: V/r -> a ?")
G=lambda w: w/(2+w)**2
for r in ("1e-2","1e-4"):
    r=mp.mpf(r); print("   r=",r," V/r =",mp.nstr(V(G,r)/r,10)," a =",mp.nstr(a_of(G),10))
