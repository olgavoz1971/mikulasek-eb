# Features of the light curves in data/

| star | A1 | A2 | fwhm1 | fwhm2 | pts/fwhm | expo/fwhm | W90/W50 | W10/W50 | Gamma est. | coverage | contact | maxima diff (mmag) | drift prim. (s) | drift sec. (s) | fold scatter/noise | trend amp. (mag) | red/white | corr. time (d) | gap steps (mmag) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AA And S16 2-min | 0.536 | 0.142 | 0.138 | 0.173 | 87.7 | 0.0107 | 0.513 | 2.11 | 0.592 | 0.745 | False | 1.8 (13 sigma) | -73 +- 77 | 143 +- 60 | 1.27 | 0.00504 | 14.9 | 0.281 | 1.3+-0.2 |
| AA And S16 30-min | 0.551 | 0.145 | 0.139 | 0.174 | 5.89 | 0.16 | 0.501 | 2.07 | 0.63 | 0.735 | False | 2.0 (3.8 sigma) | - | 7 +- 254 | 2.24 | 0.0055 | 8.77 | 0.281 | 2.1+-0.4 |
| XY Boo S50 10-min | 0.317 | 0.293 | 0.196 | 0.197 | 8.13 | 0.0956 | 0.362 | 1.89 | 0.713 | 0.833 | True | -9.7 (22 sigma) | -51 +- 12 | 9 +- 12 | 3.69 | 0.00437 | 16.6 | 0.111 | 0.0+-0.4; -1.2+-0.6; 0.1+-0.8 |
| EK Com S49 | 0.572 | 0.519 | 0.171 | 0.166 | 5.06 | 0.152 | 0.533 | 2.01 | 0.67 | 0.784 | True | 9.1 (18 sigma) | -7 +- 50 | 30 +- 27 | 2.61 | 0.0102 | 21.5 | 0.08 | 0.8+-0.4 |
| V477 Lyr S26 | 0.606 | 0.0715 | 0.0501 | 0.0199 | 16.1 | 0.0587 | 0.373 | 7.96 | 0.3 | 0.522 | False | 3.0 (1.9 sigma) | 25 +- 31 | - | 1.02 | 0.0316 | 7.44 | 1.51 | 9.8+-2.0 |
| AV Del S81 | 0.966 | 0.274 | 0.0773 | 0.102 | 90.8 | 0.00777 | 0.391 | 1.91 | 0.713 | 0.66 | False | 14.2 (37 sigma) | -26 +- 44 | 126 +- 56 | 2.04 | 0.0025 | 40.8 | 1.16 | -0.7+-0.1; 2.5+-0.2 |
| AR Aur S19 | 0.637 | 0.549 | 0.0249 | 0.0258 | 69.3 | 0.0135 | 0.317 | 1.97 | 0.63 | 0.108 | False | 0.6 (13 sigma) | -47 +- 42 | 67 +- 41 | 1.97 | 0.000626 | 50.9 | - | 0.6+-0.0 |

Notes (features that could not be measured):

- **AA And S16 2-min** (hash 1251a90d0246): gamma_estimate: poor match to profile 1 (distance 0.283)
- **AA And S16 30-min** (hash f7b28f984d12): gamma_estimate: poor match to profile 1 (distance 0.258); drift_primary: centre not measurable in both halves
- **XY Boo S50 10-min** (hash 94f918703c76): none
- **EK Com S49** (hash 7dfb8321907c): gamma_estimate: poor match to profile 1 (distance 0.267)
- **V477 Lyr S26** (hash acb7f68e0f3d): secondary_offset: error not available: bottom parabola is not convex towards fainter magnitudes; gamma_estimate: poor match to profile 1 (distance 4.803); drift_secondary: centre not measurable in both halves
- **AV Del S81** (hash 826ac61f4ff2): none
- **AR Aur S19** (hash 89109a8ca9dd): correlation_time_days: autocorrelation stays above 0.5 over half the longest segment
