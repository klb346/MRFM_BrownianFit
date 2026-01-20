# MRFM_BrownianFit 1.1.1
last update date: 01/20/2026

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

```
pip install .
```

To install the package with the optional dependencies for the MCMC module, run the following code in the root directory of the package folder:

```
pip install .[MCMC]
```

## Testing
After installation, run the example NB in the example folder. Outputs can be compared to the html run 01/16/2026 by KLB in the example_outputs folder.

Tests for use with pytest are still under development.

## brownian_fit

To use the brownian_fit class, import the class from the module as shown below

```
from MRFM_BrownianFit.brownian_fit import brownian_fit
```

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

The max_likelihood function will find the parameter values for Gamma tau0, f0, and the baseline that
correspond to the maximum of the likelihood by solving for the minimum of -log(likelihood) using the
scipy minimize function. The minimizer starts at the best fit values from lmfit and therefore one of
the fitting functions within brownian_fit must be run before calling this function. The results are
stored in a dictionary self.bayesian_result.


### Example plots:
![An example plot from the fitting of the data in examples](examples\example_outputs\20250701_baysian_8avgs.png)
![An example residual CDF plot from the fitting of the data in examples](examples\example_outputs\20250701_baysian_8avgs_residual_cdf.png)


## LabView_int

To use the LVprocessing class, import the class from the module as shown below

```
from MRFM_BrownianFit.LabView_int import LVprocessing
```

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

## MCMC

The MCMC module builds upon the bayesian inference in the brownian_fit module. 

To use the MCMC class, import the class from the module as shown below

```
from MRFM_BrownianFit.MCMC import MCMC
```

The brownian_fit class must be initialized and the least squares fitting must be performed prior to initalizing the MCMC class. The MCMC class takes in the brownian fit class as the argument.

### User Functions

The run function preforms a MCMC sampling of the log_likelihood function defined in brownian_fit using emcee. The default parameters chosen work well for Brownian's taken on the fourth generation MRFM probe built by the Marohn group at Cornell University.

```
run(self, param_bounds, walkers = 64, nsteps = 2000, progress = True, moves = emcee.moves.KDEMove(bw_method="silverman"), figpath = None, n=None, ErrorHandling = True)
```

The run function will run a MCMC sampling of the log probability defined in brownian fit, plot the walkers from the sampling, and generate the corner plots for the four fit parameters (Gamma, tau0, f0, and baseline). This function will catch errors related to the auto correlation times, attempt to rerun with longer chain lengths, and abort if the auto correlation times are not defined or too long compared to the chain length. 

Future updates will include changing the method used by emcee to attempt to fix over rejection issues with the sampler.

### Credible intervals
After running the sampling, different confidence intervals can be calculated without rerunning the sampling by using the following functions:

```
_credible_interval_68(self)
```

Calculates the 68% credible interval.

```
_credible_interval_95(self)
```

Calculates the 95% credible interval.

```
_credible_interval_n(self, n:int)
```

Calculates the n% credible interval. n must be any number between 0 and 100 (noninclusive). Floats in this range will be truncated to an integer. If n is outside this range or not of type int or float, the function will raise a ValueError and abort.

### Functions to rerun MCMC Sampling
For instances where you may only wish to rerun parts of the processing, the following functions may be used. The run function will call all of these plus the credible interval specified in the argument n.

```
_run_walkers(self, pos, ndim, param_bounds, walkers, nsteps, progress, moves)
```
May be used to rerun the sampler. To avoid memory leaks, the previous instance of the sampler should be deleted.
Use to change the initial guess, number of walkers (walkers), prior (bounds), chain length (nsteps), or move algorithm (move).

```
_plot_walkers(self, ndim, figpath)
```
Will plot the walker chains, and save if figpath is defined. To skip saving, set figpath = None.


```
_flatten samples(self)
```
Find autocorrelation times of fit parameters. Needs to be run before _gen_corner_plot()


```
_gen_corner_plot(self, figpath)
```
Generates a corner plot showing the distribution of the fit parameters.
