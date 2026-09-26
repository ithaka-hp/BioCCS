"""
Figure 1 of the BioCCS policy commentary.
Vintage and atmospheric effect of carbon removed by a forest stand
harvested in 2027 after 100 years of growth and used for BioCCS with 90% capture.

Growth:  Chapman-Richards, M(t) = Mmax * (1 - exp(-k t))**p
IRF:     Jeltsch-Thoemmes & Joos (2020), five exponential terms plus constant
Effect:  R(t) = f_cap * sum_i dM_i * IRF(t - t_i)   (Equation 3, annual pulses)
"""
import numpy as np
import matplotlib.pyplot as plt

# ---------------- parameters ----------------
T_GROW   = 100          # years of growth
YR_PLANT = 1927
YR_HARV  = YR_PLANT + T_GROW    # 2027
YR_END   = YR_HARV + 100        # 2127, end of plotting window
M_HARV   = 1000.0       # t CO2 in the stand at harvest
F_CAP    = 0.90         # capture efficiency -> 900 t stored
K, P     = 0.03, 3.0    # Chapman-Richards rate constant and shape exponent

# IRF (Jeltsch-Thoemmes & Joos 2020)
A0  = 0.008
A   = np.array([0.044, 0.112, 0.224, 0.310, 0.297])
TAU = np.array([68521.0, 5312.0, 362.0, 47.0, 6.0])

def irf(t):
    t = np.asarray(t, dtype=float)
    out = A0 + np.sum(A * np.exp(-t[..., None] / TAU), axis=-1)
    return np.where(t < 0, 0.0, out)

def chapman_richards(t, mmax, k, p):
    return mmax * (1.0 - np.exp(-k * t)) ** p

def stand(k, p, t_grow=T_GROW, m_harv=M_HARV):
    """Annual removals dM_i (t CO2) for i = 1..t_grow, scaled so that M(t_grow) = m_harv."""
    mmax = m_harv / (1.0 - np.exp(-k * t_grow)) ** p
    ages = np.arange(0, t_grow + 1)
    M = chapman_richards(ages, mmax, k, p)
    dM = np.diff(M)                       # removal in year i (between age i-1 and i)
    t_removal = np.arange(1, t_grow + 1) - 0.5  # annual removal booked at mid-year
    return mmax, ages, M, dM, t_removal

def remaining(t_eval, dM, t_removal, f_cap=1.0):
    """Reduction of atmospheric CO2 at ages t_eval attributable to removals dM at ages t_removal."""
    t_eval = np.atleast_1d(t_eval).astype(float)
    lag = t_eval[:, None] - t_removal[None, :]
    return f_cap * np.sum(dM[None, :] * irf(lag), axis=1)

# ---------------- base case ----------------
mmax, ages, M, dM, t_rem = stand(K, P)
mean_age = np.sum(t_rem * dM) / np.sum(dM)
mean_year = YR_PLANT + mean_age

t_axis = np.arange(0, YR_END - YR_PLANT + 1)          # ages 0..200
years  = YR_PLANT + t_axis
R_stored = remaining(t_axis, dM, t_rem, F_CAP)        # effect of the removal underlying the 900 t
R_cert   = np.where(t_axis >= T_GROW, F_CAP * M_HARV * irf(t_axis - T_GROW), np.nan)  # 900 t removed in 2027

R_at_harv = remaining(np.array([T_GROW]), dM, t_rem, F_CAP)[0]
R_at_end  = remaining(np.array([YR_END - YR_PLANT]), dM, t_rem, F_CAP)[0]
excess_area = np.nansum((R_cert - R_stored)[t_axis >= T_GROW])   # t CO2 * yr

print(f"Mmax = {mmax:.1f} t CO2")
print(f"Mean age of removal = {mean_age:.1f} yr  -> mean year {mean_year:.0f}")
print(f"Remaining effect of stored share at harvest ({YR_HARV}) = {R_at_harv:.1f} t CO2 "
      f"= {100*R_at_harv/(F_CAP*M_HARV):.1f}% of 900 t")
print(f"Remaining effect of stored share in {YR_END} = {R_at_end:.1f} t CO2")
print(f"Certificate (900 t in {YR_HARV}) remaining in {YR_END} = {F_CAP*M_HARV*irf(np.array([100.0]))[0]:.1f} t CO2")
print(f"Excess claim, integrated {YR_HARV}-{YR_END} = {excess_area:.0f} t CO2 yr")

# ---------------- sensitivity ----------------
print("\nSensitivity (k, p): mean age before harvest, retained share at harvest")
res = []
for k in [0.02, 0.03, 0.04, 0.05]:
    for p in [2.5, 3.0, 3.5, 4.0]:
        _, _, _, dMs, trs = stand(k, p)
        ma = np.sum(trs * dMs) / np.sum(dMs)
        ra = remaining(np.array([T_GROW]), dMs, trs, 1.0)[0] / M_HARV
        res.append((k, p, T_GROW - ma, ra))
        print(f"  k={k:.2f} p={p:.1f}: {T_GROW-ma:5.1f} yr, {100*ra:4.1f}%")
lags = [r[2] for r in res]; shares = [r[3] for r in res]
print(f"Range: {min(lags):.0f}-{max(lags):.0f} yr before harvest; {100*min(shares):.0f}-{100*max(shares):.0f}%")

# ---------------- figure ----------------
GREEN, RED, GREY = "#2e7d32", "#c62828", "#555555"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

# (a) growth
ax1.bar(YR_PLANT + t_rem, dM, width=1.0, color=GREEN, alpha=0.45, label="Annual removal")
ax1.set_ylabel("Annual removal (t CO$_2$ yr$^{-1}$)")
ax1.set_xlabel("Year")
ax1b = ax1.twinx()
ax1b.plot(YR_PLANT + ages, M, color=GREEN, lw=2, label="Carbon stock")
ax1b.set_ylabel("Carbon stock (t CO$_2$)")
ax1b.set_ylim(0, 1100)
ax1.axvline(mean_year, color=GREY, ls=":", lw=1)
ax1.text(mean_year + 1, ax1.get_ylim()[1]*0.93, f"mean year of removal\n{mean_year:.0f}", fontsize=8, color=GREY)
ax1.axvline(YR_HARV, color=GREY, lw=0.8)
ax1.text(YR_HARV - 1, ax1.get_ylim()[1]*0.93, "harvest", fontsize=8, color=GREY, ha="right")
ax1.set_xlim(YR_PLANT, YR_HARV + 3)
h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax1b.get_legend_handles_labels()
ax1.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8, frameon=False)
ax1.set_title("(a) Growth of the stand", loc="left", fontsize=10)

# (b) atmospheric effect
ax2.plot(years, R_stored, color=GREEN, lw=2,
         label="Removal underlying the stored 900 t (1927-2027)")
ax2.plot(years, R_cert, color=RED, lw=2,
         label="900 t removed in 2027, as the certificate implies")
mask = t_axis >= T_GROW
ax2.fill_between(years[mask], R_stored[mask], R_cert[mask], color=RED, alpha=0.12, lw=0,
                 label="Excess claim of the certificate")
ax2.axvline(YR_HARV, color=GREY, lw=0.8)
ax2.axhline(F_CAP * M_HARV, color=GREY, ls=":", lw=1)
ax2.text(YR_PLANT + 2, F_CAP * M_HARV + 15, "900 t captured and stored", fontsize=8, color=GREY)
ax2.annotate(f"{R_at_harv:.0f} t", xy=(YR_HARV, R_at_harv), xytext=(YR_HARV - 30, R_at_harv + 120),
             fontsize=8, color=GREEN, arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.8))
ax2.set_xlim(YR_PLANT, YR_END)
ax2.set_ylim(0, 1200)
ax2.set_xlabel("Year")
ax2.set_ylabel("Reduction of atmospheric CO$_2$ (t)")
ax2.legend(loc="upper right", fontsize=8, frameon=False)
ax2.set_title("(b) Atmospheric effect", loc="left", fontsize=10)

fig.tight_layout()
fig.savefig("Figure_1.png", dpi=300, bbox_inches="tight")
fig.savefig("Figure_1.pdf", bbox_inches="tight")
np.savetxt("Figure_1_data.csv",
           np.column_stack([years, np.r_[np.nan, dM, np.full(len(years)-len(dM)-1, np.nan)],
                            np.r_[M, np.full(len(years)-len(M), np.nan)], R_stored, R_cert]),
           delimiter=",", header="year,annual_removal_tCO2,stock_tCO2,effect_stored_share_tCO2,effect_certificate_tCO2",
           comments="", fmt="%.3f")
