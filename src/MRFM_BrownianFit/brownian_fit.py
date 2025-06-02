###################################################
# brownian_fit.py
# Defines the data class brownian_fit
# Author: Katrina L. Brown
# 2025/05/20
###################################################

class brownian_fit():
    """
    The brownian_fit class is a data class that will fit a power spectrum from 
    room temperature cantilever position fluctuations to the theoretical power 
    spectrum using lmfit. From the fit parameters, the spring constant and quality
    factor of the cantilever can be calculated.

    The power spectrum and experimental parameters need to be passed as v
    The class defines functions to run the fit with all permutations of plot? and
    report? output.

    There are two main functions to call that will call the necessary interal funtions.

    The plot_peak function will just plot the cantilever peak in a semilog plot. This is a
    useful test that the code is extracting a good region of the power spectrum.
    
    The do_fit function will fit on the cantilever peak and store the resulting parameters of
    interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
    fit report can be found from self.result['brownian'].
    The cantilever peak is plotted with the fit function, and the resisuals are plotted below.

    """
    #import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from lmfit import Model
    import math as m


    def __init__(self, data):
        """
        The __init__ function handles data input and will raise a TypeError if the
        input parameter is not a 5-element tuple.

        If the input is accepted, the x, y, temp, and N_avgs are saved to the class.
        The frequrency data is scaled by a factor of 10^-3 to convert from Hz to kHz.
        The y-data is scaled by a factor of 10^6 to convert from nm^2/Hz to pm^2/Hz.
        These unit conversions will reduce numerical errors in the fitting.

        The colors used in the plotting are also stored in a dictionary called
        'colors' which contains 'pink', 'magenta', and 'darkmagenta'.
        """        
        
        # check input 
        if type(data) == tuple and len(data) == 5:
            tupletypes = list(map(type, data))
            if tupletypes == [int, float, list, list, str]:
                self.N_avgs, self.temp, x, y, self.file= data
                self.fig_file = self.file + ".png"
                self.res_fig_file = self.file + "_residual_cdf.png"
                self.x = self.np.array(x)*10**(-3)     # [kHz]
                self.y = self.np.array(y)*10**(6)      # [pm^2/Hz]
            else:
                raise ValueError("Expecting a 5 element tuple with dtypes (int, float, list, list, str) but got "+str(tupletypes))
        else:
            raise ValueError("Expecting a 5 element tuple with (N_avgs:int, temp:float, x:list, y:list, file:str)")

        self.colors = {
            "magenta" : "#FF24D2",
            "pink" : "#FF80E5",
            "darkmagenta" : "#7D0064"
        }

    def _extract_peak(self):
        """
        The _extract_peak function extracts a 100 Hz wide section of the input data
        centered at the cantiliever peak which is found by finding the second maxium in
        the y-data (the first maximum is at zero).
        """


        #find estimate of canitlever peak (highest non-zero peak)
        f0_estimate_idx, = self.np.where(self.np.isclose(self.np.max(self.y[10:]), self.y, rtol=0.001))
        self.f0_estimate = float(self.x[int(f0_estimate_idx[0])])
        
        #truncate to 1000 Hz centered at f0_estimate **frequencies are in kHz
        #given that sampling is every 0.5 Hz, can do this by index
        
        # upper = self.f0_estimate + 0.500
        # lower = self.f0_estimate - 0.500

        # #raise TypeError(len(self.x),len(self.y), self.f0_estimate)

        # #find indices of upper and lower
        # idx_u, = self.np.where(self.np.isclose(self.x,upper, atol = 0.0005))
        # idx_l, = self.np.where(self.np.isclose(self.x,lower, atol = 0.0005))

        #raise TypeError(len(idx_u),len(idx_l), idx_u, idx_l)
    
        idx_u = int(f0_estimate_idx[0]) + 1000
        idx_l = int(f0_estimate_idx[0]) - 1000

        #truncate data set to indices
        self.x_trunc = self.x[idx_l:idx_u]
        self.y_trunc = self.y[idx_l:idx_u]

    def plot_peak(self):
        """
        The plot_peak function will just plot the cantilever peak in a semilog plot. 
        """
        self._extract_peak()
        
        self.plt.semilogy(self.x_trunc,self.y_trunc,color=self.colors["pink"], label = "data")
        #plt.title("Brownian Motion Power Spectrum")
        self.plt.ylabel('PSD [pm$^2$/Hz]')
        self.plt.xlabel('Frequency [kHz]')
        self.plt.show()

    def _baseline_avg(self):
        """
        The _baseline_avg function will estimate the baseline from the mean of region of
        the spectum at well away from cantilever peak (the last 100 data points of the data).
        """
        return self.np.mean(self.y[-100:-1])

    def _brownian(self, f:np.array, Gamma: float, tau0: float, f0: float, baseline: float):
        r"""
        The _brownian function defines the function of the sum of the power spectrum of position fluctuations of an ideal harmonic oscillator
        and frequency-independent detector noise floor to which the data is fit. 

        The power spectrum of position fluctuations [pm^2/Hz] as a function of frequency f [kHz] for a harmonic oscillator at 
        room temperature is given by

        .. math:: P_{\delta z}(f) = (\frac{kT \tau_0^2}{\Gamma}) frac{1}{((\pi * \tau_0)^4 * (f_0^2-f^2)^2 + (\pi*\tau_0)^2 * f^2))}
        with units of [m^2/Hz]

        kT: 1.38*10^-23 *295 at room temperature [J = N m]
        tau0: Ringdown time [ms]
        Gamma: Friction coefficent [N s/m]
        f0: Resonance frequency [kHz]

        It is helpful to write in terms of unitless frequencies:
        F = pi * tau0 * f
        F0 = pi * tau0 * f0

        P(f) = (kT * tau0^2 / Gamma)*(1 / ((F0^2-F^2)^2 + F^2)) with units of [pm^2/Hz] 
        **tau0 is converted to s for the unitless equation, introducing a factor of 10^-6.
        """

        kT = 1.38*10**(-23) * self.temp    # N m
        F = self.np.pi * tau0 * f         # unitless
        F0 = self.np.pi * tau0 * f0       # unitless

        prefactor = 10**18 * (kT * tau0**2 / Gamma)
        P = prefactor / ((F0**2-F**2)**2 + F**2)  #pm^2/Hz

        return P + baseline
    
    def _fit_power_spec(self):
        """
        _fit_power_spec runs the fit to the power spectrum and saves the result to self.result['brownian'].
        The residuals of the fit are calculated and saved to self.residuals
        """
        #run fit
        gmodel = self.Model(self._brownian)
        y_err = self.y_trunc/self.np.sqrt(self.N_avgs)
        w = 1/y_err

        #define inital guess of resnoance frequency and noise floor
        noise_floor = self._baseline_avg()
        f0_idx, = self.np.where(self.y_trunc == self.np.max(self.y_trunc))
        f0_init = self.x_trunc[f0_idx]

        #initalize the parameters for lmfit
        params = gmodel.make_params(Gamma=1E-11, baseline=noise_floor, f0=f0_init, tau0=200)

        #run the fit
        self.result = {}
        self.result['brownian'] = gmodel.fit(self.y_trunc, params, f=self.x_trunc, weights=w)

        #calculate residuals
        self.residuals = (self.y_trunc - self.result['brownian'].best_fit)*w
        self.resid_mean = self.np.mean(self.residuals)

    def residuals_CDF(self, figpath=None):
        r1 = self.np.sort(self.residuals)
        r2 = self.np.arange(1, len(r1)+1)/len(r1)
        fig, ax1 = self.plt.subplots(1,1,figsize=(8, 6))
        ax1.plot(r1,r2,'.', color = self.colors["darkmagenta"])
        ax1.set_ylabel('CDF')
        ax1.set_ylabel('Normalized Residuals\n[pm$^2$/Hz]')
        self.plt.tight_layout()
        if figpath != None:
            self.plt.savefig(figpath+"\\"+self.res_fig_file)
        return fig

    def plot_fit(self, figpath=None):
        """
        The _plot_fit function will plot the cantilever peak with the fit function on a semilog plot,
        and the normalized residuals underneath.
        """
        fig, (ax1, ax2) = self.plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=(3, 1))

        ax1.semilogy(self.x_trunc, self.y_trunc, color =self.colors["pink"], label='data')
        ax1.plot(self.x_trunc, self.result['brownian'].best_fit, color = self.colors["darkmagenta"], label='fit')

        ax1.set_ylabel('PSD [pm$^2$/Hz]')

        ax2.plot(self.x_trunc, self.residuals, color= self.colors["magenta"])
        ax2.set_xlabel('Frequency [kHz]')
        ax2.set_ylabel('Normalized Residuals\n[pm$^2$/Hz]')

        ax1.legend(loc='best')
        
        self.plt.tight_layout()
        if figpath != None:
            self.plt.savefig(figpath+"\"+self.fig_file)
        return fig
    
    def _find_params(self):
        """
        **need to update!!!!
        The _find_params function will calculate the cantilever spring conatant and quality factor and
        store them as self.k and self.Q respectively.
        The resonance frequency in stored in self.f0
        The force noise is found from :math: $frac{kT * \tau_0^2}{(\pi * \tau_0)^4 f_0 ^4}$ and stored
        as self.force_noise.
        The detector noise is the baseline and is stored as self.Sx
        """
        #calculate canilever spring constant in mN/m
        self.f0 = float(1.0E3*self.result['brownian'].best_values['f0'])      # Hz
        self.f0_stderr = float(1.0E3*self.result['brownian'].params['f0'].stderr) # Hz

        self.tau0 = float(1.0E-3*self.result['brownian'].best_values['tau0']) # s
        self.tau0_stderr = float(1.0E-3*self.result['brownian'].params['tau0'].stderr) # s

        self.k = float(2 * self.np.pi**2 * self.f0**2 * self.tau0 * self.result['brownian'].best_values['Gamma'] * 1e3) # mN/m
        self.k_stderr = self.k * self.m.sqrt(2*(self.tau0_stderr/self.tau0)**2 + (self.f0_stderr/self.f0)**2 + (self.result['brownian'].params['Gamma'].stderr / self.result['brownian'].best_values['Gamma'])**2)

        #calculate cantilever quality factor and propagate error
        self.Q = float(self.np.pi*self.result['brownian'].best_values['tau0']*self.result['brownian'].best_values['f0'])
        self.Q_stderr = self.Q * self.m.sqrt((self.tau0_stderr/self.tau0)**2 +(self.f0_stderr/self.f0)**2)

        #calc force noise
        kT = 1.38*10**(-23) * self.temp    # N m
        self.Sx = float(1E-6 *self.result['brownian'].best_values['baseline'])   #nm^2/Hz
        self.Sx_stderr = float(1E-6 * self.result['brownian'].params['baseline'].stderr)  #nm^2/Hz
        
        self.Sth = float((kT*self.f0)/((0.05E-12)**2 * self.k * 1E03 *self.Q))   #Hz^2/Hz  
        self.Sth_stderr = self.Sth * self.m.sqrt(2*(self.tau0_stderr/self.tau0)**2 + 4*(self.f0_stderr/self.f0)**2) #Hz^2/Hz  

        self.Sdet = (self.Sx)/(0.05**2) #Hz^-1
        self.Sdet_stderr = (self.Sx_stderr)/(0.05**2)  #Hz^-1
        # print("Cantilever quality factor = ", self.Q, " [unitless]")
        # print("Cantilever spring constant = ", self.k, " [mN/m]")

        #calc P(0)
        self.P0 = float(1E18 * kT /(self.result['brownian'].best_values['Gamma'] * self.np.pi**4 *self.tau0**2 *self.f0**4)) #nm^2/Hz
        self.P0_stderr = self.P0 * self.m.sqrt(2*(self.tau0_stderr/self.tau0)**2 + 4*(self.f0_stderr/self.f0)**2 + (self.result['brownian'].params['Gamma'].stderr / self.result['brownian'].best_values['Gamma'])**2) #nm^2/Hz
    
        #calc intrinsic dissipation
        self.intrinsic_diss = 1E9 * self.k/(self.Q*2*self.np.pi*self.f0) #pN s/m
        self.intrinsic_diss_stderr = self.intrinsic_diss * self.m.sqrt((self.k_stderr/self.k)**2 + (self.Q_stderr/self.Q)**2 + (self.f0_stderr/self.f0)**2) #pN s/m

        #calc intrinsic force noise
        self.intrinsic_force_noise = float(1E28 *4 * kT *self.intrinsic_diss) #aN^2/Hz
        self.intrinsic_force_noise_stderr = float(1E28 *4 * kT *self.intrinsic_diss_stderr)  #aN^2/Hz


    def do_fit(self):
        """
        The do_fit function will fit on the cantilever peak and store the resulting parameters of
        interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
        fit report can be found from self.result['brownian'].
        The cantilever peak is plotted with the fit function, and the resisuals are plotted below.
        """
        self._extract_peak()
        self._fit_power_spec()
        self._find_params()
    

