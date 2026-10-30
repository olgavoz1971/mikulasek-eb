# Features of the light curves in data/

| star | A1 | A2 | fwhm1 | fwhm2 | pts/fwhm | expo/fwhm | W90/W50 | W10/W50 | Gamma est. | coverage | contact | maxima diff (mmag) | drift prim. (s) | drift sec. (s) | fold scatter/noise | trend order |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AA And S16 2-min | 0.533 | 0.139 | 0.138 | 0.171 | 87.5 | 0.0108 | 0.514 | 2.07 | 0.63 | 0.927 | True | 1.8 (13 sigma) | -73 +- 77 | 107 +- 61 | 1.27 | 3 |
| AA And S16 30-min | 0.55 | 0.142 | 0.139 | 0.173 | 5.89 | 0.16 | 0.501 | 2.06 | 0.63 | 0.934 | True | 2.0 (3.8 sigma) | - | 7 +- 254 | 2.24 | 3 |
| XY Boo S50 10-min | 0.31 | 0.287 | 0.193 | 0.195 | 7.99 | 0.0971 | 0.365 | 1.87 | 0.735 | 1 | True | -9.7 (22 sigma) | -51 +- 12 | 9 +- 12 | 3.69 | 3 |
| EK Com S49 | 0.567 | 0.511 | 0.17 | 0.163 | 5.03 | 0.153 | 0.535 | 1.98 | 0.691 | 0.999 | True | 9.1 (18 sigma) | -7 +- 50 | 32 +- 26 | 2.61 | 3 |
| V477 Lyr S26 | 0.6 | - | 0.0486 | - | 15.6 | 0.0606 | 0.384 | 7.87 | 0.3 | - | - | 3.0 (1.9 sigma) | 25 +- 31 | - | 1.02 | 3 |
| AV Del S81 | 0.963 | 0.269 | 0.0772 | 0.1 | 90.6 | 0.00778 | 0.391 | 1.9 | 0.713 | 0.532 | False | 14.2 (37 sigma) | -26 +- 44 | 126 +- 56 | 2.04 | 3 |
| AR Aur S19 | 0.637 | 0.548 | 0.0249 | 0.0257 | 69.3 | 0.0135 | 0.317 | 1.97 | 0.63 | 0.152 | False | 0.6 (13 sigma) | -47 +- 42 | 67 +- 41 | 1.97 | 3 |

Notes (features that could not be measured):

- **AA And S16 2-min** (hash 6dad59351ae5): gamma_estimate: poor match to profile 1 (distance 0.273)
- **AA And S16 30-min** (hash bbe89e1146ba): gamma_estimate: poor match to profile 1 (distance 0.255); harmonic_amplitudes: too few bins outside the eclipse windows; drift_primary: centre not measurable in both halves
- **XY Boo S50 10-min** (hash 31e37dca84f5): ooe_amplitude_over_depth: eclipse windows cover the whole cycle; harmonic_amplitudes: too few bins outside the eclipse windows
- **EK Com S49** (hash 1680b925a91a): gamma_estimate: poor match to profile 1 (distance 0.261); ooe_amplitude_over_depth: eclipse windows cover the whole cycle; harmonic_amplitudes: too few bins outside the eclipse windows
- **V477 Lyr S26** (hash 6f05aec8f774): A2: no minimum deeper than the reference level; gamma_estimate: poor match to profile 1 (distance 4.717); coverage: eclipse widths not measurable; parity_secondary: width not measurable; drift_secondary: centre not measurable in both halves
- **AV Del S81** (hash 65301a737224): trend_amplitude: segments too short
- **AR Aur S19** (hash faccdde69a36): trend_amplitude: segments too short
