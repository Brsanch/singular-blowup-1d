"""Task 2: mutation tests of TB-1..TB-4.  Reads the gate source (read-only), patches a copy in memory, runs one check."""
import io, contextlib, sys, traceback

SRC = open("/Volumes/4TB SD/ClaudeCode/singular-blowup-1d/validate/edge_track_b_proof_lemmas.py").read()


def run(label, fn, reps):
    s = SRC
    for old, new in reps:
        assert old in s, f"pattern not found: {old!r}"
        s = s.replace(old, new)
    ns = {"__name__": "mut"}
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            exec(compile(s, "mut", "exec"), ns)
            ok = ns[fn]()
    except Exception as e:
        ok = f"EXC {e!r}"
    print(f"--- {label}\n    -> returns {ok}\n    {buf.getvalue().strip()}")


# ---------------- TB-1 ----------------
run("TB-1 baseline", "check_tb1", [])
run("TB-1 MUT Q -> wrong operator  Q = -2 Nr(u^2 NgI G)", "check_tb1",
    [("    Q = lambda G: -2.0 * Nr(u * NgI(G))\n    Dg = lambda G: np.fft.ifft((sigma", "    Q = lambda G: -2.0 * Nr(u * u * NgI(G))\n    Dg = lambda G: np.fft.ifft((sigma")])
run("TB-1 MUT Q -> unrelated function QG := exp(-z^2/50)*cos(z)", "check_tb1",
    [("    QG = Q(G)\n    P = -QG * Dg(G)", "    QG = np.exp(-z*z/50.0)*np.cos(z)\n    P = -QG * Dg(G)")])
run("TB-1 MUT bilinear sign: P = +QG*Dg(G) + 0.5 QG G (wrong sign of advection)", "check_tb1",
    [("P = -QG * Dg(G) + 0.5 * QG * G\n    DQG", "P = +QG * Dg(G) + 0.5 * QG * G\n    DQG")])
run("TB-1 MUT coefficient: P = -QG*Dg(G) + 0.25 QG G", "check_tb1",
    [("P = -QG * Dg(G) + 0.5 * QG * G\n    DQG", "P = -QG * Dg(G) + 0.25 * QG * G\n    DQG")])
run("TB-1 MUT data with a(G) != 0: G = w/(2+w)^2", "check_tb1",
    [("    G = Ng(h0)\n    QG = Q(G)\n    P = -QG", "    G = np.exp(z)/(2.0+np.exp(z))**2\n    QG = Q(G)\n    P = -QG")])

# ---------------- TB-2 ----------------
run("TB-2 baseline", "check_tb2", [])
run("TB-2 MUT exponent 0.25 -> 0.30", "check_tb2", [("(t ** 0.25) * norm2", "(t ** 0.30) * norm2")])
run("TB-2 MUT exponent 0.25 -> 0.27", "check_tb2", [("(t ** 0.25) * norm2", "(t ** 0.27) * norm2")])
run("TB-2 MUT exponent 0.25 -> 0.235", "check_tb2", [("(t ** 0.25) * norm2", "(t ** 0.235) * norm2")])
run("TB-2 MUT residue: phi_* -> 3*phi_* (residue -3)", "check_tb2",
    [("phi_star = 1.0 / (2.0 * np.cosh(0.5 * z))", "phi_star = 3.0 / (2.0 * np.cosh(0.5 * z))")])
run("TB-2 MUT phi_* -> 1/(u+1)  (different function, residue -1 but hat = pi/sin(2 pi m)-type)", "check_tb2",
    [("phi_star = 1.0 / (2.0 * np.cosh(0.5 * z))", "phi_star = 1.0 / (u + 1.0)")])
run("TB-2 MUT phi_* -> u/(u^2+1)^2 * u^2 ... i.e. 1/u-tail with coefficient 1 but otherwise unrelated: u^3/(u^4+1)", "check_tb2",
    [("phi_star = 1.0 / (2.0 * np.cosh(0.5 * z))", "phi_star = u**3 / (u**4 + 1.0)")])
run("TB-2 MUT line 5/8 -> 3/4 (true exponent -1/2)", "check_tb2", [("sig2 = 0.625  # 5/8", "sig2 = 0.75")])
run("TB-2 MUT line 5/8 -> 0.6 (true exponent -0.2)", "check_tb2", [("sig2 = 0.625  # 5/8", "sig2 = 0.6")])

# ---------------- TB-3 ----------------
run("TB-3 baseline", "check_tb3", [])
run("TB-3 MUT data with a(G) != 0 (pole NOT cancelled): G = w/(2+w)^2", "check_tb3",
    [("    G = np.fft.ifft(Nm1 * np.fft.fft(w1 * h0)).real / w1\n", "    G = np.exp(z)/(2.0+np.exp(z))**2\n")])
run("TB-3 MUT symbol: cot -> 3*cot", "check_tb3",
    [("QGhat1 = (cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]", "QGhat1 = 3*(cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]"),
     ("QGhat2 = (cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]", "QGhat2 = 3*(cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]")])
run("TB-3 MUT symbol: cot/(zeta+1) -> tan/(zeta+1)", "check_tb3",
    [("QGhat1 = (cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]", "QGhat1 = (1.0/cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]"),
     ("QGhat2 = (cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]", "QGhat2 = (1.0/cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]")])
run("TB-3 MUT symbol without smoothing: cot (no 1/(zeta+1))", "check_tb3",
    [("QGhat1 = (cot_c(z1) / (z1 + 1.0)) * Ghat1[keep]", "QGhat1 = (cot_c(z1)) * Ghat1[keep]"),
     ("QGhat2 = (cot_c(z2) / (z2 + 1.0)) * Ghat2[keep]", "QGhat2 = (cot_c(z2)) * Ghat2[keep]")])

# ---------------- TB-4 ----------------
run("TB-4 baseline", "check_tb4", [])
run("TB-4 MUT constant 2/3 -> 0.664", "check_tb4", [("target = (2.0 / 3.0) * c1_exact", "target = 0.664 * c1_exact")])
run("TB-4 MUT constant 2/3 -> 0.65", "check_tb4", [("target = (2.0 / 3.0) * c1_exact", "target = 0.65 * c1_exact")])
run("TB-4 MUT constant 2/3 -> 3/4", "check_tb4", [("target = (2.0 / 3.0) * c1_exact", "target = 0.75 * c1_exact")])
run("TB-4 MUT profile with a(G) != 0: G = w/(1+w)^2 - 2.0*w/(2+w)^2", "check_tb4",
    [("G_c = w / (1.0 + w) ** 2 - (2.0 * math.sqrt(2.0)) * w / (2.0 + w) ** 2", "G_c = w / (1.0 + w) ** 2 - 2.0 * w / (2.0 + w) ** 2")])
run("TB-4 MUT V_G sign", "check_tb4", [("return float(-(1.0 / (math.pi * r_val))", "return float(+(1.0 / (math.pi * r_val))")])
