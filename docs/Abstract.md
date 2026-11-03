# The Mikulášek Method for Timing Eclipsing Binaries: Implementation and Comparison with Wall-supported Local Approximations

Automation is essential for efficient eclipse timing in the era of space telescopes and large photometric surveys. However, automation also requires careful selection of appropriate approximation methods to ensure reliable measurements and avoid systematic errors.

We present a public Python implementation (GitHub link) of the Mikulášek method for phenomenological light-curve modelling and eclipse timing (Mikulášek 2015). We compare its timing results with local approximations of eclipse minima by the Wall-Supported Asymptotic Parabola (WSAP) and Wall-Supported Line (WSL) methods (Andronov 1996).

One of the features of our implementation is a semi-automatic choice of model terms. It fits a set of candidates, ranks them using the Bayesian information criterion (BIC), and flags the results for the user’s review when no candidate is clearly preferred or when the data raise other warnings. This provides an automated alternative to the manual model-selection procedure described in the original Mikulášek (2015) paper. We present preliminary tests of this feature on synthetic TESS-like light curves.

The implementation also includes several practical extensions intended to simplify its application. These include reusing a model fitted to one TESS sector for other sectors of the same star, an optional second profile component for minima with complex shapes and accounting for the smoothing of the light curve by long exposures.

We compare the timing results obtained with the Mikulášek model and the WSAP/WSL methods using synthetic light curves generated with PHOEBE and degraded to resemble TESS observations by adding noise and instrumental trends. The test set includes a spotted system, for which we compare the results with both the geometric time of conjunction and the photometric minimum.

We also apply the methods to real TESS light curves of several eclipsing binaries. These examples demonstrate the practical application of the methods to real observations and provide preliminary O–C diagrams with their associated uncertainties. We discuss the cases where each method performs better.