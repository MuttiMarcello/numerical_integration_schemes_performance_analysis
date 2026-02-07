# Numerical Methods Performance Analysis

Individual project implementing and analysing performance of diverse numerical methods for zero-finding and ODE integration.

## Overview

This repository contains a study developed as part of a university modeling and simulation project.

The work includes:
1. Implementation and performance assessment of zero-finding Newton's method
2. Implementation and performance assessment of RK2, RK4 methods for ODE initial-value-problem solution
3. Stability region analysis of RK2, RK4 methods
4. Accuracy region analysis of RK1, RK2, RK4 methods
5. Implementation and stability region analysis of BI2-theta methods
6. Implementation and performance assessment of RK4, IEX4 methods for ODE initial-value-problem solution
7. Implementation and performance assessment of AB3, AM3, ABM3, BDF3 methods for ODE initial-value-problem solution

## Results and Validation

Key results:
- RK methods struggle to capture stiff dynamics and are outperformed by IEX4, which can handle larger step sizes.
- BI-theta methods, theta>0.5, handle stiff dynamics well with larger step sizes. BI-theta methods, theta<0.5, are not A-stable, meaning they require smaller step sizes for stiff dynamics, leading to increased computational cost.
- AM3 outperforms AB3 for stiff systems, given the same step size. ABM3 combines the two methods, but an increase in performance is not guaranteed. BDF3 outperforms AM3, AB3, ABM3 for stiff systems.

Representative outputs:
-RK2, RK2, BI-theta stability regions
-RK1, RK2, RK4 accuracy regions
-AB3, AM3, ABM3, BDF3 performance comparison

Representative figures are available in 'results/' (PNG format).
See 'results/results.txt' for figure-by-figure explanations.
The full methodology and results are documented in 'docs/report.pdf'.

## Repository structure

- 'src/' - MATLAB implementations of each study
- 'docs/' - Project report
- 'results/' - Key result figures (PNG)
- 'figures/' - Source figures (EPS)

## Development notes

This repository was uploaded after project completion.
Commit history does not reflect the original development timeline.