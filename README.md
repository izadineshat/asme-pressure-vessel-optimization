# ASME Sec VIII Div 1 Pressure Vessel Structural Optimization (IRR-calibrated)

A robust computational mechanics benchmark for the optimal structural design of horizontal pressure vessels with hemispherical heads under ASME Boiler and Pressure Vessel Code (Section VIII, Division 1), fully calibrated with Iranian steel fabrication costs (IRR).

## Overview

The solver optimizes the vessel geometry and assigns commercially available structural steel plate gauges (ASTM A516-Gr70 / DIN EN 10029) to minimize overall manufacturing expenses, spanning plate procurement, rolling, head flanging, and radiographic-controlled submerged arc welding.

## Quick Start

Execute the driver script inside MATLAB:

```matlab
[best_sol, best_cost_irr, conv_history] = run_iran_asme_vessel_optimization();
```

Please replace placeholder files with the full implementation and add the ASME notes PDF under docs/ when ready.
