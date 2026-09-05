# ASME Sec VIII Div 1 Pressure Vessel Optimization (Industrial Fabrication Sizing)

A rigorous computational mechanics framework for techno-economic optimization of horizontal cylindrical pressure vessels with hemispherical heads under ASME Boiler and Pressure Vessel Code (Section VIII, Division 1), fully calibrated with Iranian steel fabrication costs and commercial plate gauges.

## Problem Overview

Standard academic benchmarks frequently optimize vessel dimensions over continuous thickness variables, producing non-manufacturable values. This repository provides an end-to-end mixed-integer nonlinear programming (MINLP) solver that maps thickness variables directly onto commercial structural steel plates (ASTM A516-Gr70 / DIN EN 10029) available in the domestic industrial steel market.

## Mathematical Formulation

$$\min_{\mathbf{x}} f(\mathbf{x}) = C_1 x_1 x_3 x_4 + C_2 x_2 x_3^2 + C_3 x_1^2 x_4 + C_4 x_1^2 x_3 \quad [\text{IRR}]$$

Subject to:
* $g_1(\mathbf{x}) = 0.02994 x_3 + 2.0 - x_1 \le 0 \quad (\text{Shell Hoop Stress - ASME UG-27})$
* $g_2(\mathbf{x}) = 0.01475 x_3 + 2.0 - x_2 \le 0 \quad (\text{Head Membrane Stress - ASME UG-32})$
* $g_3(\mathbf{x}) = 2.5 \times 10^{10} - \left(\pi x_3^2 x_4 + \frac{4}{3}\pi x_3^3\right) \le 0 \quad (\text{Min Volume: } 25\,\text{m}^3)$
* $g_4(\mathbf{x}) = x_4 - 8000 \le 0 \quad (\text{Logistics Envelope})$

### Design Variables and Commercial Bounds

| Parameter | Symbol | Classification | Domain Bounds | Discrete Step / Resolution |
|---|---|---|---|---|
| Shell Thickness ($T_s$) | $x_1$ | Discrete | $[6, 60]\,\text{mm}$ | $\{6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60\}\,\text{mm}$ |
| Head Thickness ($T_h$) | $x_2$ | Discrete | $[6, 60]\,\text{mm}$ | $\{6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60\}\,\text{mm}$ |
| Inside Radius ($R$) | $x_3$ | Continuous | $[500, 3000]\,\text{mm}$ | Continuous |
| Tangent Length ($L$) | $x_4$ | Continuous | $[1000, 8000]\,\text{mm}$ | Continuous |

## Domestic Cost Calibration (Iranian Rial)

* **$C_1$ ($36.98\,\text{IRR/mm}^3$):** Rolled alloy steel plate procurement ($750,000\,\text{IRR/kg}$).
* **$C_2$ ($118.32\,\text{IRR/mm}^3$):** Plate material including head press spinning and flanging ($1,200,000\,\text{IRR/kg}$).
* **$C_3$ ($35.00\,\text{IRR/mm}^3$):** Longitudinal seam SAW wire, flux, and edge preparation.
* **$C_4$ ($263.89\,\text{IRR/mm}^3$):** Circumferential girth joint fit-up and radiographic quality multi-pass welding.

## Repository Execution

Clone the repository and run the standalone optimization routine in MATLAB:

```matlab
[best_sol, best_cost_irr, conv_history] = run_iran_asme_vessel_optimization();```
