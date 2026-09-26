# BioCCS: vintage and atmospheric effect of the removal behind a certified unit

Calculation and figure for the policy commentary on the dating of removals
certified under the EU CRCF for BioCCS.

## What is computed

A forest stand grows for 100 years (1927-2027) along a Chapman-Richards function

    M(t) = Mmax (1 - exp(-k t))^p

scaled so that the stand holds 1,000 t CO2 at harvest, with k = 0.03 yr^-1 and p = 3
(illustrative values within the range fitted to Pinus stands). The annual increments
dM_i are the removals from the atmosphere, booked at mid-year.

At harvest, 90% of the carbon (900 t CO2) is captured and stored; 10% is lost as
capture losses.

The reduction of atmospheric CO2 attributable to the stored share at time t is
(Equation 3 of the commentary)

    R(t) = 0.9 * sum_i dM_i * IRF(t - t_i)

with the impulse response function of Jeltsch-Thoemmes and Joos (2020),

    IRF(t) = a0 + sum_{j=1..5} a_j exp(-t / tau_j)

a0 = 0.008; a = 0.044, 0.112, 0.224, 0.310, 0.297; tau = 68521, 5312, 362, 47, 6 yr.

This is compared with the same IRF applied to a single pulse of 900 t CO2 removed in
2027, which is what a certificate issued in 2027 implies.

## Results (base case)

| Quantity | Value |
|---|---|
| Mean age of removal | 49 yr after planting (mean year 1976, 51 yr before harvest) |
| Remaining effect of the stored share at harvest, 2027 | 438 t CO2 (49% of 900 t) |
| Remaining effect of the stored share in 2127 | 291 t CO2 |
| Remaining effect of a 900 t pulse removed in 2027, in 2127 | 332 t CO2 |
| Excess claim of the certificate, integrated 2027-2127 | about 11,300 t CO2 yr |

Sensitivity over k = 0.02-0.05 yr^-1 and p = 2.5-4: mean year of removal 36-67 yr
before harvest; retained effect at harvest 43-55% of the stored amount. The full
grid is printed by the script.

## Files

- `figure1_vintage.py` - calculation and figure
- `Figure_1.png`, `Figure_1.pdf` - the figure
- `Figure_1_data.csv` - annual removals, stock and both effect curves, 1927-2127
- `run_output.txt` - console output of the script, including the sensitivity grid

## Usage

    pip install numpy matplotlib
    python figure1_vintage.py

## References

Jeltsch-Thoemmes, A. & Joos, F. (2020). Modeling the evolution of pulse-like
perturbations in atmospheric carbon and carbon isotopes: the role of
weathering-sedimentation imbalances. Climate of the Past 16, 423-451.

Richards, F.J. (1959). A flexible growth function for empirical use.
Journal of Experimental Botany 10, 290-300.

Schmidt, H.-P. & Hagemann, N. (2026). Methane emissions offsetting with temporary
carbon sinks.
