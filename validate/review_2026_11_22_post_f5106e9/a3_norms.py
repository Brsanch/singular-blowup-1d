"""Review script (2026-11-22, post-f5106e9 review). Lemma A3 (Track A) norm-compatibility probe.
Model element F_K with ||F_K||_{X_rho} = e^{rho K} (Mellin data concentrated at |xi| ~ K).
dW(s') = psi(s') F_K. The doc's W-component weight is (d^3/s^2)^{1/2}||dW||_{L^2([0,s];X_rho)};
the G-component weight is (d/s)^{1/2}||dG(s)||_{X_rho}. The map dW -> dG = eps*int_0^s dW (no derivative
loss at all) must be bounded between these norms for the claimed contraction. Ratio = G-norm/(W-norm/gamma) / eps."""
import numpy as np
rho0, a, S = 1.0, 1.0/6, 2.0
s = np.linspace(1e-4, S, 4000)
for K in [10, 30, 100, 300]:  # K >= 1000 overflows float64
    psi = np.exp(a*s*K - rho0*K) * K**2 * s          # saturates the W-constraint near the cone edge
    I2 = np.concatenate([[0], np.cumsum(0.5*(psi[1:]**2+psi[:-1]**2)*np.diff(s))])
    I1 = np.concatenate([[0], np.cumsum(0.5*(psi[1:]+psi[:-1])*np.diff(s))])
    Wn = Gn = 0.0
    for j in range(len(s)):
        rmax = rho0 - a*s[j]
        rho = np.linspace(0, rmax, 3000, endpoint=False)
        d = rho0 - rho - a*s[j]
        Wn = max(Wn, np.max(d**1.5/s[j]*np.exp(rho*K))*np.sqrt(I2[j]))
        Gn = max(Gn, np.max(np.sqrt(d/s[j])*np.exp(rho*K))*I1[j])
    print(f"K={K:5d}  ||dG||_G / ||dW||_W (per eps*gamma) = {Gn/Wn:10.3f}   sqrt(K)={np.sqrt(K):.2f}")
