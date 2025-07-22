# MRFM_BrownianFit
The brownian fit package has been built to enable the extraction of cantilever parameters
from the position fluctuation power spectrum for MRFM cantilevers. The current version
houses two modules: brownian_fit and LabView_int. The brownian_fit module has functions
defined for fitting and plotting brownian motion data. The LabView_int module has
functions defined to convert input data into the format that the brownian_fit module
expects.

This package was built to run with python 3.9.12. The author recommends installing in a
new virtual environment to protect any existing dependency installations.

## Installation
To install the package, run the following code in the root directory of the package folder:

'''
pip install .
'''

Tests can be run by running pytest in the root directory of the package
The tests simply verify that the code is running as expected by comparing the output of the
example data to previously performed fits.

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

### Example plots:
![An example plot from the fitting of the data in examples](examples\fit_spec_ex.png)
![An example residual CDF plot from the fitting of the data in examples](examples\example_20250314_residual_cdf.png)


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

