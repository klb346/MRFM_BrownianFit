# brownian_fit_class
The brownian fit package has been built to enable the extraction of cantilever parameters from the position fluctuation power spectrum for MRFM cantilevers. The current version (0.0.1) is a proof of concept for the fitting procedure. Future versions will be integratable into current data collection workflows. 

To install package run the following code in the root directory of the package folder:

'''
pip install .
'''

Tests can be run by running pytest in the root directory of the package
The tests simply verify that the code is running as expected by comparing the output of the example data to previously performed fits.

The 0.0.2 version of brownian_fit houses the definition of a data class that can be used to plot and fit the power spectrum from room temperature position fluctuations for an MRFM cantiliver, and calculate the quality facotr and spring constant of the cantilever from this fit.

To use the class, import the class from the module

'''
from brownian_fit_class.brownian_fit_class import brownian_fit
'''

Call the class with a 5-element tuple:
(N_avgs:int, temp:float, x:list, y:list, file)

The plot_peak function will just plot the cantilever peak in a semilog plot. This is a
useful test that the code is extracting a good region of the power spectrum.
    
The do_fit function will fit on the cantilever peak and store the resulting parameters of
interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
fit report can be found from self.result['brownian'].
The cantilever peak is plotted with the fit function, and the resisuals are plotted below.

Example plot:
![An example plot from the fitting of the data in examples](examples\fit_spec_ex.png)