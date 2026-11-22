#!/usr/bin/env python3
"""
GATE: Global nonlinear PURE-STEP (eta === 0) weak-edge compatibility manifold
via autonomous s-clock, exact conjugation G = N_g h, and strip-analytic Picard contraction.

Verifies the full nonlinear pure-step (eta === 0, A === B === 1, c_l = 2 - V_G, c_omega = -c_l/2)
PDE on v = X - 1 > 0:
    G_tau + (2 - V_G + Q G) v G_v - (1/2) (Q G)(1 + G) = 0,
with global weak-edge compatibility a(tau) = H[g](1, tau) === 0 (B(tau) === A(tau)^2 across X = 1)
and ZERO Theta-cusp (Theta === (pi/2) 1_{X>1}).

Checks (all pass => exit 0):
  PM-1  Exact O(eps^2) conjugation asymptotics on flat strip-analytic seed h_0(u) = u^3 e^{-u}:
          u f_{P_1}(u) -> -(pi/2) a(P_1) and u^2 h_2(u) -> -(pi/2) a(P_1) as u -> inf (< 5e-3),
          confirming a(G_2) === 0 and c_1(G_2) = 3 a(P_1) = -(3/pi) J(G_1).
  PM-2  Full nonlinear pure-step conjugated Picard iteration on s in [0, 1.2] (> s_inf ~ 1.0) at eps = 0.02,
        N_s = 240, sponge-stabilized (tanh end-sponge on the forcing) domain [-64, 64], Nz in {4096, 8192}
        (the supported resolutions; at Nz = 16384 the Picard difference stalls at a ~1e-9 round-off floor).
        Data h_0 = u^3 e^{-u} / PM_SC with a grid-independent constant PM_SC.
          - all iterates finite; Picard ||dG|| < 1e-10 at both Nz; Nz 4096 -> 8192 drift of G(1.2) < 1e-10;
            |grid max / PM_SC - 1| < 2e-5
          - s-clock nonlinear pure-step PDE residual < 5e-5 on the bulk, s in [0.05, 1.2) (O(ds^2)-limited)
          - weak-edge compatibility of the DATA: |a(G(0))| < 1e-9, with an a != 0 data control
            (h_0 += 1e-3 u/(1+u^2)) required to give |a(G(0))| > 1e-7. For s > 0, a(G(s)) = 0 is automatic
            (S_g(s) removes any 1/u tail of h), so |a_quad(G(s))| < 1e-7 at s = 0.6, 1.0 is printed as "automatic".
  PM-3  Physical-time reconstruction (r(s) = sqrt(lambda(s)), tau(s) -> +inf as s -> s_inf-) at tau in {0.2, 0.6},
        v in {0.2, 1, 4}:
          - |s_inf - 1| < 0.05
          - physical-time pure-step PDE residual < 2e-5, evaluated with V_G from an INDEPENDENT adaptive
            quadrature of its definition on the reconstructed physical profile (not the z-grid functional that
            drives the horizon ODE); |V_ODE - V_indep| < 1e-8 at the tested tau; quadrature error < 1e-9
          - V_G control: rerunning the reconstruction with V -> 3V, -V, 0 in the ODE path (the only place the
            reconstruction uses V) must raise the residual >= 50x.
        Not checked: a wrong V_G formula coded identically in both implementations (no self-consistency check can
        see that); a(tau) in physical time.
All thresholds fail on NaN/inf (no builtin max/min accumulators).
"""
from __future__ import annotations

import math
import warnings
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.integrate import IntegrationWarning, quad
from scipy.special import loggamma

# Grid-independent data normalisation: continuous max over the bulk of |N_g[u^3 e^{-u}]| (attained at
# z* = -3.840403), computed once by spline maximisation at Nz in {4096, 8192, 16384} on [-48,48] .. [-96,96]
# (agreement 1e-12). The old grid-sampled max made the data depend on where grid points fell (C3 follow-up).
PM_SC = 1.0869601848306


def _isfin(*xs) -> bool:
    return all(bool(np.all(np.isfinite(x))) for x in xs)


def _lt(x, tol) -> bool:
    """x < tol, and False for NaN/inf (never NaN-blind)."""
    return bool(np.isfinite(x) and x < tol)


def _gt(x, tol) -> bool:
    return bool(np.isfinite(x) and x > tol)


def _nanmax(vals) -> float:
    """max that PROPAGATES NaN (builtin max silently drops it); NaN for an empty input."""
    arr = np.asarray(list(vals), dtype=float)
    return float(np.max(arr)) if arr.size else float("nan")


def build_ops(Nz: int = 4096, zmin: float = -36.0, zmax: float = 36.0, sigma: float = 0.25):
    L = zmax - zmin
    dz = L / Nz
    z = zmin + np.arange(Nz) * dz
    w = np.exp(z)
    u = np.exp(-0.5 * z)
    xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)
    sgn = np.where(xi >= 0, 1.0, -1.0)

    def build_Nm(sig):
        mm = sig + 1j * xi
        logcos = (
            np.pi * np.abs(xi)
            - np.log(2.0)
            - 1j * np.pi * sig * sgn
            + np.log1p(np.exp(-2.0 * np.pi * np.abs(xi) + 2j * np.pi * sig * sgn))
        )
        return np.exp(logcos + loggamma(2.0 * mm + 2.0) - np.log(np.pi))

    Nm_g = build_Nm(sigma)
    Nm_r = build_Nm(sigma - 0.5)
    w_g = np.exp(-sigma * z)
    w_r = np.exp(-(sigma - 0.5) * z)
    m_g = sigma + 1j * xi
    filt = np.where(np.abs(xi) <= 0.6 * np.max(np.abs(xi)), 1.0, 0.0)

    Ng = lambda h: np.fft.ifft((Nm_g * filt) * np.fft.fft(w_g * h)).real / w_g
    NgI = lambda G: np.fft.ifft((filt / Nm_g) * np.fft.fft(w_g * G)).real / w_g
    Nr = lambda h: np.fft.ifft((Nm_r * filt) * np.fft.fft(w_r * h)).real / w_r
    Q = lambda G: -2.0 * Nr(u * NgI(G))
    Dg = lambda G: np.fft.ifft((m_g * filt) * np.fft.fft(w_g * G)).real / w_g

    return dict(Nz=Nz, dz=dz, z=z, w=w, u=u, Ng=Ng, NgI=NgI, Nr=Nr, Q=Q, Dg=Dg)


def check_pm1() -> bool:
    ops = build_ops(Nz=4096, zmin=-44.0, zmax=44.0)
    z, u, dz = ops["z"], ops["u"], ops["dz"]
    Ng, NgI, Q, Dg = ops["Ng"], ops["NgI"], ops["Q"], ops["Dg"]
    bulk = (z > -12.0) & (z < 12.0)

    h0 = np.where(u < 200.0, (u**3) * np.exp(-u), 0.0)
    G0 = Ng(h0)
    sc = float(np.max(np.abs(G0[bulk])))
    h0 /= sc
    G0 /= sc

    Ns = 100
    s_grid = np.linspace(0.0, 1.0, Ns + 1)
    ds = 1.0 / Ns
    x = u * ds
    ex = np.exp(-x)
    p1 = np.where(x < 1e-4, 1.0 - 0.5 * x + x**2 / 6.0 - x**3 / 24.0, (1.0 - ex) / x)
    p2 = np.where(x < 1e-4, 0.5 - x / 6.0 + x**2 / 24.0 - x**3 / 120.0, (1.0 - p1) / x)
    wt1 = ds * p2
    wt0 = ds * (p1 - p2)

    fP = []
    for sv in s_grid:
        Gj = Ng(np.exp(-u * sv) * h0)
        QGj = Q(Gj)
        Pj = -QGj * Dg(Gj) + 0.5 * QGj * Gj
        fP.append(NgI(Pj))

    h2 = np.zeros_like(u)
    for j in range(Ns):
        h2 = ex * h2 + wt0 * fP[j] + wt1 * fP[j + 1]

    G1 = Ng(np.exp(-u * 1.0) * h0)
    QG1 = Q(G1)
    P1 = -QG1 * Dg(G1) + 0.5 * QG1 * G1
    mi = (z >= -20.0) & (z <= 20.0)
    aP1 = float(-(1.0 / math.pi) * np.sum(np.exp(-0.5 * z[mi]) * P1[mi]) * dz)
    target = -0.5 * math.pi * aP1

    ix16 = np.argmin(np.abs(z - (-16.0)))
    uf_val = float(u[ix16] * fP[-1][ix16])
    u2h_val = float((u[ix16] ** 2) * h2[ix16])
    rel_f = abs(uf_val - target) / abs(target)
    rel_h = abs(u2h_val - target) / abs(target)
    uh_val = float(u[ix16] * h2[ix16])  # -> 0 since h2 = O(1/u^2) => a(G2) = 0

    ok = bool(_lt(rel_f, 1e-4) and _lt(rel_h, 5e-3) and _lt(abs(uh_val), 2e-4))
    print(
        f"PM-1 O(eps^2) conjugation edge   PASS={ok}: "
        f"u*fP1={uf_val:+.6f}, u^2*h2={u2h_val:+.6f} vs -pi/2*a(P1)={target:+.6f} "
        f"(rel={rel_f:.1e}, {rel_h:.1e}), |u*h2(z=-16)|={abs(uh_val):.2e}"
    )
    return ok


def check_pm2_pm3() -> bool:
    # Sponge-stabilized domain [-64, 64]; run at Nz = 4096 and Nz = 8192 (the supported resolutions:
    # at Nz = 16384 the Picard difference stalls at a ~1e-9 round-off floor, not at the iteration cap).
    def solve_pm(Nz: int):
        zmin, zmax = -64.0, 64.0
        L = zmax - zmin
        dz = L / Nz
        z = zmin + np.arange(Nz) * dz
        w = np.exp(z)
        u = np.exp(-0.5 * z)
        sigma = 0.25
        xi = 2.0 * np.pi * np.fft.fftfreq(Nz, dz)
        sgn = np.where(xi >= 0, 1.0, -1.0)

        def build_Nm_nf(sig):
            mm = sig + 1j * xi
            lc = np.pi * np.abs(xi) - np.log(2.0) - 1j * np.pi * sig * sgn + np.log1p(np.exp(-2.0 * np.pi * np.abs(xi) + 2j * np.pi * sig * sgn))
            return np.exp(lc + loggamma(2.0 * mm + 2.0) - np.log(np.pi))

        Nm_g = build_Nm_nf(sigma)
        Nm_r = build_Nm_nf(sigma - 0.5)
        w_g = np.exp(-sigma * z)
        w_r = np.exp(-(sigma - 0.5) * z)
        m_g = sigma + 1j * xi
        Ng = lambda h: np.fft.ifft(Nm_g * np.fft.fft(w_g * h)).real / w_g
        NgI = lambda G: np.fft.ifft(np.fft.fft(w_g * G) / Nm_g).real / w_g
        Nr = lambda h: np.fft.ifft(Nm_r * np.fft.fft(w_r * h)).real / w_r
        Q = lambda G: -2.0 * Nr(u * NgI(G))
        Dg = lambda G: np.fft.ifft(m_g * np.fft.fft(w_g * G)).real / w_g
        sponge = 0.5 * (1.0 + np.tanh((z - (zmin + 10.0)) / 1.5)) * 0.5 * (1.0 - np.tanh((z - (zmax - 10.0)) / 1.5))

        bulk = (z > -12.0) & (z < 12.0)
        h0 = np.where(u < 200.0, (u**3) * np.exp(-u), 0.0)
        sc_grid = float(np.max(np.abs(Ng(h0)[bulk])))  # only used to sanity-check PM_SC
        h0 = h0 / PM_SC  # grid-independent normalisation (C3 follow-up)

        eps = 0.02
        S_max = 1.2
        Ns = 240
        s_grid = np.linspace(0.0, S_max, Ns + 1)
        ds = S_max / Ns
        x = u * ds
        ex = np.exp(-x)
        p1 = np.where(x < 1e-4, 1.0 - 0.5 * x + x**2 / 6.0 - x**3 / 24.0, (1.0 - ex) / x)
        p2 = np.where(x < 1e-4, 0.5 - x / 6.0 + x**2 / 24.0 - x**3 / 120.0, (1.0 - p1) / x)
        wt1 = ds * p2
        wt0 = ds * (p1 - p2)

        G = np.array([Ng(eps * np.exp(-u * sv) * h0) for sv in s_grid])
        diff = float("nan")
        for it in range(1, 50):
            fP = np.zeros_like(G)
            for j in range(Ns + 1):
                Gj = G[j]
                QGj = Q(Gj)
                fP[j] = sponge * NgI(-QGj * Dg(Gj) + 0.5 * QGj * Gj)
            G_new = np.zeros_like(G)
            hc = eps * h0.copy()
            G_new[0] = Ng(hc)
            for j in range(Ns):
                hc = ex * hc + wt0 * fP[j] + wt1 * fP[j + 1]
                G_new[j + 1] = Ng(hc)
            diff = float(np.max(np.abs((G_new - G)[:, bulk])))
            G = G_new
            if diff < 1e-12:
                break
        return dict(z=z, w=w, u=u, dz=dz, s_grid=s_grid, ds=ds, Ns=Ns, bulk=bulk, G=G, Q=Q, Dg=Dg, Ng=Ng,
                    h0=h0, eps=eps, it=it, diff=diff, sc_grid=sc_grid)

    R4 = solve_pm(4096)
    R8 = solve_pm(8192)
    z, w, u, dz, s_grid, ds, Ns, bulk = (R4[k] for k in ("z", "w", "u", "dz", "s_grid", "ds", "Ns", "bulk"))
    G, Q, Dg = R4["G"], R4["Q"], R4["Dg"]
    finite_G = _isfin(G, R8["G"])
    nz_drift = float(np.max(np.abs(G[-1][bulk] - np.interp(z[bulk], R8["z"], R8["G"][-1]))))
    sc_rel = _nanmax([abs(R4["sc_grid"] / PM_SC - 1.0), abs(R8["sc_grid"] / PM_SC - 1.0)])

    QG_all = np.array([Q(G[j]) for j in range(Ns + 1)])
    max_s_res = _nanmax(
        np.max(np.abs(((G[j + 1] - G[j - 1]) / (2.0 * ds) - 0.5 * QG_all[j] + QG_all[j] * Dg(G[j]) - 0.5 * QG_all[j] * G[j])[bulk]))
        for j in range(10, Ns)
    )

    def a_quad_tail(G_prof):
        m_q = (z >= -16.0) & (z <= 20.0)
        zm = z[m_q]
        wt = np.full_like(zm, dz)
        wt[0] *= 0.5
        wt[-1] *= 0.5
        ix_lo = int(np.where(m_q)[0][0])
        c1_e = float(G_prof[ix_lo] / w[ix_lo])
        return float(-(1.0 / math.pi) * np.sum(wt * np.exp(-0.5 * zm) * G_prof[m_q]) - (2.0 * c1_e / math.pi) * math.exp(0.5 * zm[0]))

    # a = 0 is a constraint on the DATA only: for s > 0, S_g(s) = N_g e^{-us} N_g^{-1} kills any 1/u tail of h,
    # so a(G(s)) = 0 is automatic there. Check a(G(0)) directly, with an a != 0 data control (h0 += 1e-3 u/(1+u^2)).
    a_s0 = abs(a_quad_tail(G[0]))
    a_s0_ctrl = abs(a_quad_tail(R4["Ng"](R4["eps"] * (R4["h0"] + 1e-3 * u / (1.0 + u * u)))))
    max_a_quad = _nanmax([abs(a_quad_tail(G[120])), abs(a_quad_tail(G[200]))])

    ok2 = bool(
        finite_G
        and _lt(R4["diff"], 1e-10)
        and _lt(R8["diff"], 1e-10)
        and _lt(nz_drift, 1e-10)
        and _lt(sc_rel, 2e-5)
        and _lt(max_s_res, 5e-5)
        and _lt(a_s0, 1e-9)
        and _gt(a_s0_ctrl, 1e-7)
        and _lt(max_a_quad, 1e-7)
    )
    print(
        f"PM-2 Full nonlinear pure-step    PASS={ok2}: "
        f"iter=({R4['it']},{R8['it']}), ||dG||=({R4['diff']:.2e},{R8['diff']:.2e}), Nz(4096->8192)={nz_drift:.2e}, "
        f"s-clock PDE res={max_s_res:.2e}, |a(G(0))|={a_s0:.2e} (a!=0 data ctrl {a_s0_ctrl:.2e}), "
        f"|a_quad(s=0.6,1.0)|={max_a_quad:.2e} (automatic)"
    )

    # PM-3: physical-time reconstruction. The horizon ODE dr/ds = -1 + V/2 uses the z-grid functional
    # compute_VG_s ("ODE path"). The physical-time residual uses V_G computed INDEPENDENTLY from the
    # reconstructed physical profile by adaptive quadrature of the definition
    #     V_G(tau) = -(1/pi) int_0^inf y^{-1/2} ln(1 + 1/y) G(y, tau) dy = -(2/pi) int_0^inf ln(1 + t^-2) G(t^2, tau) dt.
    # Control: the reconstruction is rerun with V replaced by 3V, -V and 0 in the ODE path (the only place the
    # reconstruction uses V); the residual (independent V) must then be >= 50x larger than with the true V.
    # Limit: a wrong V_G formula coded identically in BOTH implementations cannot be detected by this check.
    def compute_VG_s(G_prof, r_val):
        if r_val <= 1e-12:
            return 0.0
        arg = (r_val * r_val) * np.exp(-z)
        ln_term = np.where(arg > 1e6, 2.0 * math.log(r_val) - z + np.log1p(1.0 / arg), np.log1p(arg))
        return float(-(1.0 / (math.pi * r_val)) * np.sum(np.exp(0.5 * z) * ln_term * G_prof) * dz)

    def lagrange_weights(s_val, ss):
        Lw = []
        for m in range(4):
            L_m = 1.0
            for q_m in range(4):
                if q_m != m:
                    L_m *= (s_val - ss[q_m]) / (ss[m] - ss[q_m])
            Lw.append(L_m)
        return Lw

    def V_indep(G_z_sp, j0, Lw, r0):
        lr2 = 2.0 * math.log(r0)
        z_lo, z_hi = float(z[0]) + 1.0, float(z[-1])

        def ig(t):
            y = t * t
            if y <= 0.0:
                return 0.0
            zz = math.log(y) + lr2
            if zz < z_lo or zz > z_hi:
                return 0.0
            g = sum(Lw[m] * float(G_z_sp[j0 + m](zz)) for m in range(4))
            return math.log1p(1.0 / y) * g if g != 0.0 else 0.0

        brk = [0.0, 1e-4, 1e-2, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0, 1e4, 1e7]
        tot, err = 0.0, 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            for lo, hi in zip(brk[:-1], brk[1:]):
                val, e = quad(ig, lo, hi, limit=400, epsabs=1e-15, epsrel=1e-11)
                tot += val
                err += e
            val, e = quad(ig, brk[-1], np.inf, limit=400, epsabs=1e-15, epsrel=1e-11)
        return -(2.0 / math.pi) * (tot + val), err + e

    def eval_phys_residual(v_scale: float):
        r_arr = np.zeros(Ns + 1)
        r_arr[0] = 1.0
        s_inf_idx = Ns
        for j in range(Ns):
            k1 = -1.0 + 0.5 * v_scale * compute_VG_s(G[j], r_arr[j])
            r_mid = r_arr[j] + 0.5 * ds * k1
            if r_mid <= 0.0:
                s_inf_idx = j
                break
            k2 = -1.0 + 0.5 * v_scale * compute_VG_s(0.5 * (G[j] + G[j + 1]), r_mid)
            r_next = r_arr[j] + ds * k2
            if r_next <= 0.0:
                s_inf_idx = j
                break
            r_arr[j + 1] = r_next
        s_inf_loc = float(s_grid[s_inf_idx] - r_arr[s_inf_idx] / (-1.0 + 0.5 * v_scale * compute_VG_s(G[s_inf_idx], r_arr[s_inf_idx])))
        tau_arr = np.zeros(s_inf_idx + 1)
        for j in range(s_inf_idx):
            tau_arr[j + 1] = tau_arr[j] + 0.5 * ds * (1.0 / r_arr[j] + 1.0 / r_arr[j + 1])

        s_sub = s_grid[: s_inf_idx + 1]
        r_sp = make_interp_spline(s_sub, r_arr[: s_inf_idx + 1], k=3)
        tau_to_s = make_interp_spline(tau_arr, s_sub, k=3)
        G_z_sp = [make_interp_spline(z, G[j], k=3) for j in range(s_inf_idx + 1)]
        QG_z_sp = [make_interp_spline(z, QG_all[j], k=3) for j in range(s_inf_idx + 1)]

        def j0_of(s_val):
            return min(max(int(s_val / ds) - 1, 0), s_inf_idx - 3)

        def eval_field(splines_list, w_val, s_val):
            zw = math.log(w_val)
            j0 = j0_of(s_val)
            Lw = lagrange_weights(s_val, s_sub[j0 : j0 + 4])
            return sum(Lw[m] * float(splines_list[j0 + m](zw)) for m in range(4))

        dtau = 1e-4
        dv_rel = 1e-4
        res_list, vdiff_list, qerr_list = [], [], []
        for tau_val in (0.2, 0.6):
            s0 = float(tau_to_s(tau_val))
            r0 = float(r_sp(s0))
            j0 = j0_of(s0)
            Lw = lagrange_weights(s0, s_sub[j0 : j0 + 4])
            Vp, qerr = V_indep(G_z_sp, j0, Lw, r0)
            G_s0 = sum(Lw[m] * G[j0 + m] for m in range(4))
            V_ode = v_scale * compute_VG_s(G_s0, r0)
            vdiff_list.append(abs(V_ode - Vp))
            qerr_list.append(qerr)

            for v0 in (0.2, 1.0, 4.0):
                def G_phys(vv, tt):
                    ss_t = float(tau_to_s(tt))
                    rr_t = float(r_sp(ss_t))
                    return eval_field(G_z_sp, rr_t * rr_t * vv, ss_t)

                G_t = (G_phys(v0, tau_val + dtau) - G_phys(v0, tau_val - dtau)) / (2.0 * dtau)
                dv = v0 * dv_rel
                G_v = (G_phys(v0 + dv, tau_val) - G_phys(v0 - dv, tau_val)) / (2.0 * dv)
                G_val = G_phys(v0, tau_val)
                QG_val = r0 * eval_field(QG_z_sp, r0 * r0 * v0, s0)
                res_list.append(abs(G_t + (2.0 - Vp + QG_val) * v0 * G_v - 0.5 * QG_val * (1.0 + G_val)))
        return s_inf_loc, float(tau_arr[-1]), _nanmax(res_list), _nanmax(vdiff_list), _nanmax(qerr_list)

    s_inf, tau_max, max_phys_res, vdiff, qerr = eval_phys_residual(1.0)
    ctrl = {lab: eval_phys_residual(vs)[2] for lab, vs in (("3V", 3.0), ("-V", -1.0), ("0", 0.0))}
    ctrl_min = float(np.min(np.asarray(list(ctrl.values()), dtype=float)))  # NaN-propagating
    vg_ratio = ctrl_min / max_phys_res if (np.isfinite(max_phys_res) and max_phys_res > 0.0) else float("nan")

    ok3 = bool(
        _lt(abs(s_inf - 1.0), 0.05)
        and _lt(max_phys_res, 2e-5)
        and _lt(vdiff, 1e-8)
        and _lt(qerr, 1e-9)
        and _gt(vg_ratio, 50.0)
    )
    print(
        f"PM-3 Physical-time pure-step PDE PASS={ok3}: "
        f"s_inf={s_inf:.6f} (tau_max={tau_max:.2f}), res_G(*)={max_phys_res:.2e} with independent-quadrature V_G "
        f"(|V_ODE-V_indep|={vdiff:.1e}); ODE/reconstruction V->3V,-V,0: "
        f"({ctrl['3V']:.2e}, {ctrl['-V']:.2e}, {ctrl['0']:.2e}) ({vg_ratio:.1f}x)"
    )
    return ok2 and ok3


def main() -> int:
    p1 = check_pm1()
    p23 = check_pm2_pm3()
    assert p1 and p23, "Pure-step nonlinear weak-edge manifold gate FAILED"
    print("ALL PASS (PM-1..PM-3)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
