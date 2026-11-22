#!/usr/bin/env python3
"""Task 6: symmetric-part norm of A_b = N^{-1} b (sigma+d_z) N on H^k, k=0..3 (gate's matrix, vectorised),
for smooth / rough / narrow b and growing frequency cutoff.  Also the hierarchical weight sum c_m xi^{2m}."""
import numpy as np, math, wk
sigma = 0.25; L = 20.0
def sym_norms(bhat, N, keepfrac=0.4, weights=("k0", "k1", "k2", "k3", "hier")):
    dz = L / N
    xi = 2 * np.pi * np.fft.fftfreq(N, dz)
    Nm = wk.build_Nm(xi, sigma)
    idx = np.where(np.abs(xi) <= keepfrac * np.max(np.abs(xi)))[0]
    B = bhat[(idx[:, None] - idx[None, :]) % N]
    out = {}
    for wn in weights:
        if wn == "hier":
            W = np.sqrt(64 + 16 * xi**2 + 4 * xi**4 + xi**6)
        else:
            W = (1 + xi**2) ** (0.5 * int(wn[1]))
        col = (sigma + 1j * xi[idx]) * Nm[idx] / W[idx]
        row = W[idx] / Nm[idx]
        M = B * col[None, :] * row[:, None]
        S = 0.5 * (M + M.conj().T)
        out[wn] = float(np.max(np.abs(np.linalg.eigvalsh(S))))
    return out, float(np.max(np.abs(xi[idx])))
def bnorms(bhat, N):
    xi = 2 * np.pi * np.fft.fftfreq(N, L / N)
    return (math.sqrt(L * np.sum((1 + xi**2) ** 2 * np.abs(bhat) ** 2)), float(np.sum(np.abs(xi) * np.abs(bhat))),
            math.sqrt(L * np.sum((1 + xi**2) ** 3 * np.abs(bhat) ** 2)), math.sqrt(L * np.sum((1 + xi**2) ** 3.5 * np.abs(bhat) ** 2)))
def mk(name, N):
    z = -L / 2 + np.arange(N) * L / N
    xi = 2 * np.pi * np.fft.fftfreq(N, L / N)
    if name == "gate":
        return np.fft.fft(np.exp(-0.25 * z * z) * np.cos(z)) / N
    if name.startswith("rough"):          # bhat = <eta>^-(s+0.55): b in H^s but barely
        s = float(name[5:]); bh = (1 + xi**2) ** (-0.5 * (s + 0.55)); return bh / np.sum(bh)   # b(0)=1 normalisation
    if name.startswith("narrow"):
        d = float(name[6:]); return np.fft.fft(np.exp(-0.5 * ((z) / d) ** 2)) / N
print("== (i)/(iii) cutoff growth, fixed b.  columns: N, xi_cut, ||sym|| in H^0,H^1,H^2,H^3, hierarchical-weight ; ||b||_H2, ||xi bhat||_l1")
for name in ("gate", "rough2", "rough2.5", "rough3.5"):
    for N in (256, 512, 1024, 2048, 4096):
        bh = mk(name, N); o, xc = sym_norms(bh, N); bn = bnorms(bh, N)
        print(f"{name:9s} N={N:5d} xi_cut={xc:6.1f}: " + " ".join(f"{k}={o[k]:.4f}" for k in o) + f" | ||b||H2={bn[0]:.3f} ||xi bhat||l1={bn[1]:.3f} | ratio sym_k/||b||H2: " + " ".join(f"{o[k]/bn[0]:.4f}" for k in o))
print("== (ii) narrow Gaussian b=exp(-z^2/2d^2), N=2048 (xi_cut=128.7): does sym_k/||b||_H2 stay bounded as d->0 ?")
for d in (2.0, 1.0, 0.5, 0.25, 0.125, 0.0625):
    N = 2048; bh = mk(f"narrow{d}", N); o, xc = sym_norms(bh, N); bn = bnorms(bh, N)
    print(f"d={d:6.4f}: " + " ".join(f"{k}={o[k]:.4f}" for k in o) + f" | ||b||H2={bn[0]:.3f} H3={bn[2]:.3f} H3.5={bn[3]:.3f} | sym_k/||b||H2: " + " ".join(f"{o[k]/bn[0]:.4f}" for k in o) + f" | sym_k0/||b||H3.5={o['k0']/bn[3]:.4f} sym_k1/||b||H3={o['k1']/bn[2]:.4f}")
print("== full operator norm ||A_b|| (not symmetrised), gate b, H^3: " + " ".join(
    (lambda N: (lambda xi, Nm, idx: f"N={N}:{np.linalg.norm(((np.fft.fft(np.exp(-0.25*(-L/2+np.arange(N)*L/N)**2)*np.cos(-L/2+np.arange(N)*L/N))/N)[(idx[:,None]-idx[None,:])%N])*((sigma+1j*xi[idx])*Nm[idx]/(1+xi[idx]**2)**1.5)[None,:]*((1+xi[idx]**2)**1.5/Nm[idx])[:,None],2):.2f}")(
        2*np.pi*np.fft.fftfreq(N, L/N), wk.build_Nm(2*np.pi*np.fft.fftfreq(N, L/N), sigma), np.where(np.abs(2*np.pi*np.fft.fftfreq(N, L/N)) <= 0.4*np.pi*N/L)[0]))(N) for N in (256, 512, 1024)))
# FS-2: exact symbol infimum
xi = np.linspace(-200, 200, 4000001)
zc = (xi**2 - 1/16) + 0.5j * xi
c = [64.0, 16.0, 4.0, 1.0]
lhs = sum(c[m] * (zc**m).real for m in range(4)); rhs = sum(c[m] * (xi**2 + 1/16) ** m for m in range(4))
r = lhs / rhs; k = int(np.argmin(r))
print(f"== FS-2 exact Fourier symbol (phi = e^(-z/4) v): inf over all v of lhs/rhs = min_xi = {r[k]:.6f} at xi={xi[k]:+.4f}; unweighted c_m=1: {np.min(sum((zc**m).real for m in range(4))/sum((xi**2+1/16)**m for m in range(4))):.6f}")
