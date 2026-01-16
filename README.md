---
title: MRFM_BrownianFit 1.1.1
date: 08/5/2025
---

The brownian fit package has been built to enable the extraction of cantilever parameters
from the position fluctuation power spectrum for MRFM cantilevers. The current version
houses three modules: brownian_fit, LVprocessing, and MCMC. The brownian_fit module has functions
defined for fitting and plotting brownian motion data, as well as functions for . The LVprocessing module has
functions defined to convert input data from LabView into the format that the brownian_fit module
expects. The MCMC module defines the functions needed to run additional bayesian inference such as
sampling the log likelihood surface using emcee methods to find credible intervals.

This package was built to run with python 3.9.12. The author recommends installing in a
new virtual environment to protect any existing dependency installations.

 Versions of LabView running 32-bit python will not be able to run methods in the MCMC module. The dependencies for this module are optional and can be installed with the package by adding the optional key word 'MCMC" (see installation below).

## Installation
To install the package for general use, run the following code in the root directory of the package folder:

'''
pip install .
'''

To install the package with the optional dependencies for the MCMC module, run the following code in the root directory of the package folder:

'''
pip install .[MCMC]
'''

## Testing
TBD

## brownian_fit

To use the brownian_fit class, import the class from the module as shown below

'''
from MRFM_BrownianFit.brownian_fit import brownian_fit
'''

Call the class with a 5-element tuple:
(N_avgs:int, temp:float, x:list, y:list, file)
The class assumes that x is given in units of Hz and y is given in units of nm^2/Hz. To
avoid numerical errors, the data is converted into kHz and pm^2/Hz respectively.

### User Functions:

The plot_peak function will just plot the cantilever peak in a semilogy plot. This is a useful
test that the code is extracting a good region of the power spectrum.
    
The do_fit function will perform a four pass fit on the cantilever peak and store the
resulting parameters of interest as self.k, self.Q, self.f0, self.force_noise, and
self.detector_noise. The full fit report can be found from self.result[\'brownian'].
The cantilever peak is plotted with the fit function, and the residuals are plotted below.

The residuals_cdf function will plot the CDF of the fit residuals along with the cdf for
a normal distribution with unit width.

The max_logL function will...

The MCMC will...

### Example plots:
![An example plot from the fitting of the data in examples](examples\example_outputs\20250701_baysian_8avgs.png)
![An example residual CDF plot from the fitting of the data in examples](examples\example_outputs\20250701_baysian_8avgs_residual_cdf.png)


## LabView_int

To use the brownian_fit class, import the class from the module as shown below

'''
from MRFM_BrownianFit.LabView_int import LVprocessing
'''

The LVprocessing class takes the following arguments:

n_avg:int, temp:float, x:list, y:list, name:str, path:str
If x and y are given as type ndarray, the _compile_for_fitting function will convert them to
lists to pass to brownian_fit

and the following optional arguments:
fit_range_L, fit_range_U

### User Functions

The do_fit function will format and pass the data to the brownian_fit class and then perform
a 4-pass fit on the cantilever peak and store the resulting parameters of interest as self.k,
self.Q, self.f0, self.force_noise, and self.detector_noise. The full fit report can be found
from self.result[\'brownian'].


The plot_fit function will fit on the cantilever peak and store the resulting parameters of
interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
fit report can be found from self.result[\'brownian'].
A .tex file is generated summarizing the result of the fit. Two figures are included: the 
cantilever peak data and fit with normalized residuals below, and the cumulative distribution 
function of the normalized residuals.

