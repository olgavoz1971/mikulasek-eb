# 

1. # The Mikulášek Model for Eclipsing-Binary Timing: Implementation and Comparison with WSAP/WSL  

2. # Mikulášek and WSAP/WSL eclipse timing compared: an open implementation, synthetic and TESS data 

3. # The Mikulášek Model for Timing Eclipsing Binaries: Implementation and Comparison

Several well-established methods are available for measuring times of minima of eclipsing binaries, including the Kwee–van Woerden method, the wall-supported asymptotic-parabola (WSAP) and flat-bottom (WSL) methods of Andronov \[refs\] , and the phenomenological model of Mikulášek (2015). We present a public Python implementation (GitHub link) of the Mikulášek model and its application to timing of minima.

One of the features of the implementation is a semi-automatic model-selection procedure. It fits a set of candidate models, ranks them using the Bayesian information criterion (BIC), and flags the result for the user's review when no candidate is clearly preferred or when the data raise other warnings. This provides an automated alternative to the manual model-selection procedure described in the original paper. We present preliminary tests of this feature on synthetic TESS-like light curves.

The implementation also includes several practical extensions and features intended to simplify its application to survey data. These include reusing a model fitted to a high-cadence light curve as a template for fitting lower-cadence observations of the same system, an optional second profile component to describe distorted minima, and accounting for the finite exposure time of long-cadence observations.

We compare the timing results obtained with the Mikulášek model and the WSAP/WSL methods using synthetic light curves generated with PHOEBE and degraded to resemble TESS observations by adding noise and instrumental trends. The test set includes a spotted system, for which we compare the results with both the geometric time of conjunction and the photometric minimum.

We also apply the methods to real TESS light curves of NSV 807 and \[two other stars\] from multiple sectors. These examples demonstrate the practical application of the methods to real data and provide preliminary O–C diagrams with their associated uncertainties. We discuss cases where each method performs better.