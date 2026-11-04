# Noise levels for TESS-like synthetic light curves: Tmag 11 at 120, 600 and 1800 s

Notes prepared 3 November 2026 · Setup for degrading PHOEBE light curves before comparing Mikulášek and WSL/WSAP timings

## Choice

Use a single magnitude, Tmag 11. It is far enough from the bright-star systematic floor that the methods are still separated by noise, and it is close to AA And at about Tmag 10.5, judging from its measured noise.

## Noise to add to the fluxes (Tmag 11, typical star, σ of 1 hour = 438 ppm)

| Cadence | Scaling factor √(3600 s / Δt) | σ per point, relative flux | Percent | Equivalent in mag |
| --- | --- | --- | --- | --- |
| 120 s | 5.48 | 2400 ppm | 0.24% | 2.60 mmag |
| 600 s | 2.45 | 1070 ppm | 0.11% | 1.16 mmag |
| 1800 s | 1.41 | 620 ppm | 0.06% | 0.67 mmag |

For reference, at Tmag 10 (σ of 1 hour = 241 ppm) the values are 1320 ppm, 590 ppm and 340 ppm, which is 1.43, 0.64 and 0.37 mmag.

**Source of the noise model.** The 1-hour noise of TESS 2-min targets is fitted as σ(1 hr) = a + b·10^(0.2(T−10)) + c·10^(0.4(T−10)) in the TESS exoplanet yield study (arXiv:2202.03656). For the median star the constants are a = 50.2, b = 97.4 and c = 92.9 ppm. The cadence scaling is the white-noise rule σ(point) = σ(1 hr)·√(3600 s / Δt), and the conversion to magnitudes is a factor 1.0857.

## How to apply it

- Draw Gaussian noise with the σ above on the normalised flux. Do this after the exposure integration and the instrumental trend, then convert to magnitudes.
- In code: `sigma_rel = 438e-6 * np.sqrt(3600 / cadence_s)`, then `flux_noisy = flux * (1 + rng.normal(0, sigma_rel, len(flux)))`.

## What to expect

The white-noise timing precision for a canonical eclipse (depth 0.30 mag, FWHM 0.07 d) comes out at about 7–8 s at all three cadences. The cadences therefore differ mainly in shape resolution and exposure smearing, not in timing S/N. That keeps the comparison between the three cadences clean.

## Caveats

- These are the median values of the published noise model for 2-min targets. FFI photometry at 600 s and 1800 s may run somewhat noisier than the 1/√Δt scaling predicts, depending on the aperture. Check this against one of your own FFI light curves with the out-of-eclipse method: the scatter of successive differences divided by √2.
- At a given Tmag the real scatter is about a factor of 1.5–2 either side of the median.
