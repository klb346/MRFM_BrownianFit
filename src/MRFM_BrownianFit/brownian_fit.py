###################################################
# brownian_fit.py
# Defines the data class brownian_fit
# Author: Katrina L. Brown
# 2025/05/20

# Updated 2025/07/22
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

    There are two main functions to call that will call the necessary internal functions.

    The plot_peak function will just plot the cantilever peak on a semilogy plot. This is a
    useful test that the code is extracting a good region of the power spectrum.
    
    The do_fit function will fit on the cantilever peak and store the resulting parameters of
    interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
    fit report can be found from self.result['leastsq'].

    The plot_fit function will generate a figure with the cantilever peak plotted with the fit 
    function, and the residuals plotted below.

    The residuals_CDF function will plot the cumulative distribution function of the normalized 4
    residuals **on top of the expected CDF for a normal distribution

    """
    #import libraries
    import matplotlib.pyplot as plt
    import numpy as np
    from lmfit import Model
    import math as m
    import os
    import scipy.stats as ss
    import scipy.special as ssp
    import emcee
    import tqdm
    import os
    import corner



    def __init__(self, data, scale = False):
        """
        The __init__ function handles data input and will raise a TypeError if the
        input parameter is not a 5-element tuple.

        If the input is accepted, the x, y, temp, and N_avgs are saved to the class.
        The frequency data is scaled by a factor of 10^-3 to convert from Hz to kHz.
        The y-data is scaled by a factor of 10^6 to convert from nm^2/Hz to pm^2/Hz.
        These unit conversions will reduce numerical errors in the fitting.

        The colors used in the plotting are also stored in a dictionary called
        'colors' which contains 'pink', 'magenta', and 'darkmagenta'.
        """        
        
        # check input 
        if type(data) == tuple and len(data) == 5:
            tupletypes = list(map(type, data))
            if tupletypes == [int, float, list, list , str]:
                self.N_avgs, self.temp, x, y, self.file= data

                # raise ValueError(len(x), len(y))

                self.fig_file = self.file + ".png"
                self.res_fig_file = self.file + "_residual_cdf.png"
                self.x = self.np.array(x)*10**(-3)     # [kHz]
                self.y = self.np.array(y, dtype = float)*1E6      # [pm^2/Hz]
            else:
                raise ValueError("Expecting a 5 element tuple with dtypes (int, float, list, list, str) but got "+str(tupletypes))
        else:
            raise ValueError("Expecting a 5 element tuple with (N_avgs:int, temp:float, x:list, y:list, file:str)")

        self.colors = {
            "magenta" : "#FF24D2",
            "pink" : "#FF80E5",
            "darkmagenta" : "#7D0064"
        }

        self.scale = scale
        

    def _extract_peak(self, rangeL = None, rangeU = None):
        """
        The _extract_peak function extracts a 1000 Hz wide section of the input data
        centered at the cantilever peak which is found by finding the maximum in
        the y-data (the first maximum is at zero).
        """

        if (rangeL != None) and (rangeU != None):
            # manual-fit fit whole data range
            idx_u=-1
            idx_l=0

            #find indices of fit range
            indices_above, = self.np.where(self.x >= rangeL)
            indices_below, = self.np.where(self.x <= rangeU)
            # raise ValueError(indices_above.shape, indices.shape)
            idx_u = indices_below[-1]
            idx_l = indices_above[0]
        
        else:
            #auto-fit
            #find estimate of canitlever peak (highest non-zero peak)
            # raise ValueError(self.np.shape(self.y), len(self.y), type(self.y))
            f0_estimate_idx = self.np.argmax(self.y)
            # raise ValueError(str(self.np.shape(f0_estimate_idx)))
            self.f0_estimate = float(self.x[int(f0_estimate_idx)])
            
            #truncate to 1000 Hz centered at f0_estimate **frequencies are in kHz
            #given that sampling is every 0.5 Hz, can do this by index
            
            # upper = self.f0_estimate + 0.500
            # lower = self.f0_estimate - 0.500
        
            idx_u = int(f0_estimate_idx) + 1000
            idx_l = int(f0_estimate_idx) - 1000
            if idx_l < 0:
                idx_l = 0
            if (idx_u+1) > len(self.y):
                idx_u = -1

        #truncate data set to indices
        self.x_trunc = self.x[idx_l:idx_u]
        self.y_trunc = self.y[idx_l:idx_u]
        # raise ValueError(len(self.y_trunc), f0_estimate_idx, idx_l, idx_u)

    def plot_peak(self, rangeL = None, rangeU = None):
        """
        The plot_peak function will just plot the cantilever peak in a semilog plot. 
        """
        self._extract_peak(rangeL = rangeL, rangeU = rangeU)

        
        self.plt.semilogy(self.x_trunc,self.y_trunc,color=self.colors["pink"], label = "data")
        #plt.title("Brownian Motion Power Spectrum")
        self.plt.ylabel('PSD [pm$^2$/Hz]')
        self.plt.xlabel('Frequency [kHz]')
        self.plt.show()

    def _baseline_avg(self):
        """
        The _baseline_avg function will estimate the baseline from the mean of region of
        the spectrum at well away from cantilever peak (the last 100 data points of the data).
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
        Gamma: Friction coefficient [N s/m]
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
    
    def _fit_power_spec_constant_noise_floor(self, w=[]):
        #run fit
        gmodeli = self.Model(self._brownian)
        if len(w)==0:
            y_err = self.y_trunc/self.np.sqrt(self.N_avgs)
            w = 1/y_err

        #define initial guess of resonance frequency and noise floor
        self.noise_floor_est = self._baseline_avg()
    
        f0_idx, = self.np.where(self.np.isclose(self.y_trunc,self.np.max(self.y_trunc)))
        f0_init = self.x_trunc[f0_idx]

        #initialize the parameters for lmfit
        params = gmodeli.make_params(Gamma=1E-11, f0=f0_init, tau0=200)

        #run the fit
        self.result = {}
        self.result['pass1'] = gmodeli.fit(self.y_trunc, params, f=self.x_trunc, weights=w, baseline=self.noise_floor_est)

    def _fit_power_spec(self, w=[]):
        """
        _fit_power_spec runs the fit to the power spectrum and saves the result to self.result['leastsq'].
        The residuals of the fit are calculated and saved to self.residuals
        """
        #run fit
        gmodel = self.Model(self._brownian)
        if len(w)==0:
            y_err = self.y_trunc/self.np.sqrt(self.N_avgs)
            w = 1/y_err

        #define initial guess of resnoance frequency and noise floor
        noise_floor = self._baseline_avg()
        # raise ValueError(self.np.shape(self.y_trunc), len(self.y_trunc))
    
        f0_idx, = self.np.where(self.np.isclose(self.y_trunc,self.np.max(self.y_trunc)))
        # raise ValueError(self.np.shape(f0_idx))
        f0_init = self.x_trunc[f0_idx]

        #initalize the parameters for lmfit
        params = gmodel.make_params(Gamma=1E-11, baseline=noise_floor, f0=f0_init, tau0=200)

        #run the fit
        self.result = {}
        self.result['leastsq'] = gmodel.fit(self.y_trunc, params, f=self.x_trunc, weights=w)
        self.fit_y = self.np.array(self.result['leastsq'].best_fit)

        #calculate residuals
        self.residuals = (self.y_trunc - self.result['leastsq'].best_fit)*w
        self.resid_mean = self.np.mean(self.residuals)


        # calculate residual CDF
        self.r1 = self.np.sort(self.residuals)
        self.r2 = self.np.arange(1, len(self.r1)+1)/len(self.r1)

        # calculate normal CDF
        self.x_range = self.np.linspace(-10,10,10000)
        self.norm_cdf = (1 + self.ssp.erf(self.x_range/self.np.sqrt(2)))/2

    # def _two_pass_fit(self):
    #     """
    #     _two_pass_fit runs the fit to the power spectrum and uses the initial fit to get a better estimate of y_std.
    #     The second fitting pass uses the updated standard deviation to fit the data.
    #     """
    #     #run initial fit
    #     self._fit_power_spec()

    #     #recalc weights - will override initial fit
    #     y_err = self.result['leastsq'].best_fit/self.np.sqrt(self.N_avgs)
    #     w_second = 1/y_err
    #     self._fit_power_spec(w=w_second)

    # def _three_pass_fit(self):

    #     # run initial fit with fixed baseline
    #     self._fit_power_spec_constant_noise_floor() 

    #     # update weights
    #     y_err = self.result['pass1'].best_fit/self.np.sqrt(self.N_avgs)
    #     w_second = 1/y_err
    #     # run second fit with updated weights
    #     self._fit_power_spec(w = w_second)

    #     #update weights
    #     y_err = self.result['leastsq'].best_fit/self.np.sqrt(self.N_avgs)
    #     w_third = 1/y_err

    #     #run third fit
    #     self._fit_power_spec(w = w_third)

    def _four_pass_fit(self):
        if self.scale == True:
            self.scale_factor = 10**(round(self.np.log(1/self._baseline_avg())))
            self.y_trunc = self.y_trunc * self.scale_factor

       # run initial fit with fixed baseline
        self._fit_power_spec_constant_noise_floor() 

        # update weights
        y_err = self.result['pass1'].best_fit/self.np.sqrt(self.N_avgs)
        w_second = 1/y_err
        # run second fit with updated weights and fixed baseline
        self._fit_power_spec_constant_noise_floor(w = w_second)

        #update weights
        y_err = self.result['pass1'].best_fit/self.np.sqrt(self.N_avgs)
        w_third = 1/y_err

        #run third fit with floating baseline
        self._fit_power_spec(w = w_third)

        #last update of weights
        self.y_err_f = self.result['leastsq'].best_fit/self.np.sqrt(self.N_avgs)
        w_fourth = 1/self.y_err_f

        #final fit
        self._fit_power_spec(w = w_fourth)

        if self.scale == True:
            self.y_trunc = self.y_trunc / self.scale_factor
            self.result['leastsq'].best_fit = self.result['leastsq'].best_fit / self.scale_factor

            self.result['leastsq'].best_values['baseline'] = self.result['leastsq'].best_values['baseline'] / self.scale_factor
            self.result['leastsq'].params['baseline'].value = self.result['leastsq'].params['baseline'].value / self.scale_factor
            self.result['leastsq'].params['baseline'].stderr = self.result['leastsq'].params['baseline'].stderr / self.scale_factor
            self.result['leastsq'].best_values['Gamma'] = self.result['leastsq'].best_values['Gamma'] * self.scale_factor
            self.result['leastsq'].params['Gamma'].stderr = self.result['leastsq'].params['Gamma'].stderr * self.scale_factor

            self.residuals = self.residuals/self.scale_factor
            self.resid_mean = self.resid_mean/self.scale_factor



    def residuals_CDF(self, figpath=None):
        fig, ax1 = self.plt.subplots(1,1,figsize=(8, 6))
        ax1.plot(self.x_range, self.norm_cdf, "-", color = self.colors["pink"])
        ax1.plot(self.r1,self.r2,'.', color = self.colors["magenta"])
        ax1.axes.text(0.8,0.8, "Normalized Residual Mean = "+str(round(self.resid_mean,4)))
        ax1.set_ylabel('CDF')
        ax1.set_xlabel('Normalized Residual [pm$^2$/Hz]')
        ax1.set_xlim((self.r1[0], self.r1[-1]))
        self.plt.tight_layout()
        if figpath != None:
            self.plt.savefig(self.os.path.join(figpath,self.res_fig_file))
        return fig

    def plot_fit(self, figpath=None):
        """
        The _plot_fit function will plot the cantilever peak with the fit function on a semilogy plot,
        and the normalized residuals underneath.
        """
        fig, (ax1, ax2) = self.plt.subplots(2, 1, figsize=(8, 6), sharex=True, height_ratios=(3, 1))

        ax1.semilogy(self.x_trunc, self.y_trunc, color =self.colors["pink"], label='data')
        ax1.plot(self.x_trunc, self.result['leastsq'].best_fit, color = self.colors["darkmagenta"], label='fit')

        ax1.set_ylabel('PSD [pm$^2$/Hz]')

        ax2.plot(self.x_trunc, self.residuals, color= self.colors["magenta"])
        ax2.set_xlabel('Frequency [kHz]')
        ax2.set_ylabel('Normalized Residuals\n[pm$^2$/Hz]')

        ax1.legend(loc='best')
        
        self.plt.tight_layout()
        if figpath != None:
            self.plt.savefig(self.os.path.join(figpath, self.fig_file))
        # self.plt.show()
        return fig
    
    def _find_params(self):
        """
        **need to update!!!!
        The _find_params function will calculate the cantilever spring constant and quality factor and
        store them as self.k and self.Q respectively.
        The resonance frequency in stored in self.f0
        The force noise is found from :math: $frac{kT * \tau_0^2}{(\pi * \tau_0)^4 f_0 ^4}$ and stored
        as self.force_noise.
        The detector noise is the baseline and is stored as self.Sx
        """

        #calculate canilever spring constant in mN/m
        self.f0 = float(1.0E3*self.result['leastsq'].best_values['f0'])      # Hz
        self.f0_stderr = float(1.0E3*self.result['leastsq'].params['f0'].stderr) # Hz

        self.tau0 = float(1.0E-3*self.result['leastsq'].best_values['tau0']) # s
        self.tau0_stderr = float(1.0E-3*self.result['leastsq'].params['tau0'].stderr) # s

        self.k = float(2 * self.np.pi**2 * self.f0**2 * self.tau0 * self.result['leastsq'].best_values['Gamma'] * 1e3) # mN/m
        self.k_stderr = self.k * self.m.sqrt(2*(self.tau0_stderr/self.tau0)**2 + (self.f0_stderr/self.f0)**2 + (self.result['leastsq'].params['Gamma'].stderr / self.result['leastsq'].best_values['Gamma'])**2)

        #calculate cantilever quality factor and propagate error
        self.Q = float(self.np.pi*self.result['leastsq'].best_values['tau0']*self.result['leastsq'].best_values['f0'])
        self.Q_stderr = self.Q * self.m.sqrt((self.tau0_stderr/self.tau0)**2 +(self.f0_stderr/self.f0)**2)

        #calc force noise
        kT = 1.38*10**(-23) * self.temp    # N m
        self.Sx = float(1E-6 *self.result['leastsq'].best_values['baseline'])   #nm^2/Hz
        self.Sx_stderr = float(1E-6 * self.result['leastsq'].params['baseline'].stderr)  #nm^2/Hz
        
        self.Sth = float((kT*self.f0)/((0.05E-12)**2 * self.k * 1E03 *self.Q))   #Hz^2/Hz  
        self.Sth_stderr = self.Sth * self.m.sqrt(2*(self.tau0_stderr/self.tau0)**2 + 4*(self.f0_stderr/self.f0)**2) #Hz^2/Hz  

        self.Sdet = (self.Sx)/(0.05**2) #Hz^-1
        self.Sdet_stderr = (self.Sx_stderr)/(0.05**2)  #Hz^-1
        # print("Cantilever quality factor = ", self.Q, " [unitless]")
        # print("Cantilever spring constant = ", self.k, " [mN/m]")

        #calc P(0)
        self.P0 = float(1E18 * kT /(self.result['leastsq'].best_values['Gamma'] * self.np.pi**4 *self.tau0**2 *self.f0**4)) #nm^2/Hz
        self.P0_stderr = self.P0 * self.m.sqrt(2*(self.tau0_stderr/self.tau0)**2 + 4*(self.f0_stderr/self.f0)**2 + (self.result['leastsq'].params['Gamma'].stderr / self.result['leastsq'].best_values['Gamma'])**2) #nm^2/Hz
    
        #calc intrinsic dissipation
        self.intrinsic_diss = 1E9 * self.k/(self.Q*2*self.np.pi*self.f0) #pN s/m
        self.intrinsic_diss_stderr = self.intrinsic_diss * self.m.sqrt((self.k_stderr/self.k)**2 + (self.Q_stderr/self.Q)**2 + (self.f0_stderr/self.f0)**2) #pN s/m

        #calc intrinsic force noise
        self.intrinsic_force_noise = float(1E28 *4 * kT *self.intrinsic_diss) #aN^2/Hz
        self.intrinsic_force_noise_stderr = float(1E28 *4 * kT *self.intrinsic_diss_stderr)  #aN^2/Hz


    def do_fit(self):
        """
        The do_fit function will perform a 4-pass fit on the cantilever peak and store the resulting parameters of
        interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
        fit report can be found from self.result['leastsq'].

        """
        self._extract_peak()

        self._four_pass_fit()

        self._find_params()
    
    # define functions for bayesian fit - do NOT change previous functions

    def _log_likelihood(self, theta, x, y):
        #unpack theta
        Gamma, tau0, f0, baseline = theta

        model = self._brownian(x, Gamma, tau0, f0, baseline)
        
        if self.N_avgs == 1:
            return self.np.sum(
            -self.N_avgs*self.np.log(model)
            -(self.N_avgs*y/(model)))

        return self.np.sum(
        -((self.N_avgs-1)*self.np.log((self.N_avgs-1))-(self.N_avgs-1))
        +self.N_avgs*self.np.log(self.N_avgs)
        -self.N_avgs*self.np.log(model)
        +(self.N_avgs-1)*self.np.log(y)
        -(self.N_avgs*y/(model)))
    
    def _log_prior(self, theta, param_bounds):
        Gamma, tau0, f0, baseline = theta

        Gamma_min, Gamma_max = param_bounds[0]
        tau0_min, tau0_max = param_bounds[1]
        f0_min, f0_max = param_bounds[2]
        baseline_min, baseline_max = param_bounds[3]

        if (Gamma_min < Gamma < Gamma_max and 
            tau0_min < tau0 < tau0_max and 
            f0_min < f0 < f0_max and 
            baseline_min < baseline < baseline_max):
            return 0.0
        
        return -self.np.inf
    
    def _log_probability(self, theta, param_bounds, x, y):
        lp = self._log_prior(theta, param_bounds)
        
        if not self.np.isfinite(lp):
            return -self.np.inf
        return lp + self._log_likelihood(theta, x, y)
    
    def max_likelihood(self, param_bounds):
        """
        The max_likelihood function will find the parameter values for Gamma tau0, f0, and the baseline that
        corrpond to the maximum of the likelihood by solving for the minimum of -log(likelihood) using the
        scipy minimize function. The minimizer starts at the best fit values from lmfit and therefore one of
        the fitting functions within brownian_fit must be run before calling this function. The results are
        stored in a dictionary self.bayesian_result.
        """
        from scipy.optimize import minimize
        soln = minimize(
            lambda *args: -self._log_likelihood(*args),
            self.np.array(
                [
                    self.result['leastsq'].best_values['Gamma'], 
                    self.result['leastsq'].best_values['tau0'], 
                    self.result['leastsq'].best_values['f0'], 
                    self.result['leastsq'].best_values['baseline']
                    ])*(1+0.02*self.np.random.randn(1,4)),
            args=(self.x_trunc, self.y_trunc),
            method = 'Nelder-Mead')
        
        self.bayesian_result = {'Gamma': soln.x[0],
                                'tau0': soln.x[1],
                                'f0': soln.x[2],
                                'baseline': soln.x[3],
                                'message': soln.message}
                        

    def MCMC(self, param_bounds, walkers = 64, nsteps = 2000, progress = True, moves = emcee.moves.KDEMove(bw_method="silverman"), figpath = None):
        """
        The MCMC function will run a Markhov Chain Monte Carlo simulation
        
        """
        
        # define inital state within 5% of leastsq fit best values
        pos = self.np.array([
            self.result['leastsq'].best_values['Gamma'], 
            self.result['leastsq'].best_values['tau0'], 
            self.result['leastsq'].best_values['f0'], 
            self.result['leastsq'].best_values['baseline']
             ])*(1+0.05*self.np.random.randn(walkers,4))

        # store ndim
        nwalkers, ndim = pos.shape

        #initalize and run sampler
        print("Running emcee sampler...")
        sampler = self.emcee.EnsembleSampler(walkers, ndim, self._log_probability, moves = moves, args=(param_bounds, self.x_trunc, self.y_trunc))
        sampler.run_mcmc(initial_state=pos, nsteps=nsteps, progress=progress)
        
        print("Plotting walkers...")

        # plot walker path
        fig1, axes = self.plt.subplots(4, figsize=(6.50, 8.0), sharex=True)
        samples = sampler.get_chain()
        labels = ["Gamma", "tau0", "f0", "baseline"]
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("step number");

        fig1.align_ylabels(axes)
        fig1.subplots_adjust(hspace=0.1)
        self.plt.tight_layout()

        if figpath != None:
            self.plt.savefig(self.os.path.join(figpath, (self.file + '_mcmc_walkers.png')), dpi=300)
            self.plt.savefig(self.os.path.join(figpath, (self.file + '_mcmc_walkers.pdf')))

        self.plt.show()

        print("Calculating autocorrelation times...")
        # get autocorrelation time
        self.bayesian_result["tau"]= sampler.get_autocorr_time()
        print(self.bayesian_result["tau"])

        flat_samples = sampler.get_chain(discard=2* int(self.np.max(self.bayesian_result["tau"])), thin=15, flat=True)

        print("Generating corner plots...")
        # plot parameter distributions
        fig2 = self.corner.corner(
            flat_samples, labels=labels, truths=[
                    self.result['leastsq'].best_values['Gamma'], 
                    self.result['leastsq'].best_values['tau0'], 
                    self.result['leastsq'].best_values['f0'], 
                    self.result['leastsq'].best_values['baseline']
                    ])
        if figpath != None:
            fig2.savefig(self.os.path.join(figpath, (self.file + 'mcmc_corner.png')), dpi=300)
            fig2.savefig(self.os.path.join(figpath, (self.file + 'mcmc_corner.pdf')))
        
        print("Done.")