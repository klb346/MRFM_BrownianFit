###################################################
# MCMC.py
# Defines the class for running MCMC sampling for the log likelihood surface of brownian cantilever spectra using the emcee module.
# Requires that the class brownian_fit be initalized, and that the least square fitting be preformed prior to use 
# Author: Katrina L. Brown
# 2026/01/15
###################################################

class MCMC():
    '''
    Docstring for MCMC
    '''

    # import libraries
    from MRFM_BrownianFit.brownian_fit import brownian_fit as bf
    
    try:
        import emcee
    except ModuleNotFoundError: 
        print("emcee is not installed.")
        raise
    try:
        import tqdm
    except ModuleNotFoundError: 
        print("tqdm is not installed.")
        raise
    
    try:
        import corner
    except ModuleNotFoundError: 
        print("corner is not installed.")
        raise
        

    def __init__(self, fit_result, ErrorHandling = True):
        print ("Initializing class MCMC...")
        self.fit_result = fit_result
        self.ErrorHandling = ErrorHandling
        self.proc_prog = 0

        if ErrorHandling == True:
            print("Auto error handling is turned on, error correcting steps will be run. This can be turned off by setting self.ErrorHandling to False")
        else:
            print("Auto error handling is turned off. This can be turned on by setting self.ErrorHandling to True")


    # def _run_walkers(self, pos, ndim, param_bounds, walkers, nsteps, progress, moves):
    def _run_walkers(self):

        # define inital state within 5% of leastsq fit best values
        self.pos = self.bf.np.array([
                self.fit_result.result['leastsq'].best_values['Gamma'], 
                self.fit_result.result['leastsq'].best_values['tau0'], 
                self.fit_result.result['leastsq'].best_values['f0'], 
                self.fit_result.result['leastsq'].best_values['baseline']
                ])*(1+0.05*self.bf.np.random.randn(self.walkers,4))            
        # store ndim
        nwalkers, ndim = self.pos.shape
        self.ndim = ndim 
        del nwalkers

        print("Running emcee sampler...")
        # print(param_bounds)
        # print(self.fit_result.x_trunc)
        # print(self.fit_result.y_trunc)
        self.sampler = self.emcee.EnsembleSampler(self.walkers, self.ndim, self.fit_result._log_probability, moves = self.moves, args=(self.param_bounds, self.fit_result.x_trunc, self.fit_result.y_trunc))
        self.sampler.run_mcmc(initial_state=self.pos, nsteps=self.nsteps, progress=self.progress)

        # update progress var
        self.proc_prog = 1

    def _plot_walkers(self, figpath):

        if self.proc_prog == 0:
            print("Walkers must be run before plotting.")
            if self.ErrorHandling == True:
                print("Auto handling error...")

                self._run_walkers()
            else:
                raise RuntimeError("Failed to plot walkers as sampler was not run.")
                
        print("Plotting walkers...")

        # plot walker path
        fig1, axes = self.bf.plt.subplots(4, figsize=(6.50, 8.0), sharex=True)
        samples = self.sampler.get_chain()
        labels = ["Gamma", "tau0", "f0", "baseline"]
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)

        axes[-1].set_xlabel("step number");

        fig1.align_ylabels(axes)
        fig1.subplots_adjust(hspace=0.1)
        self.bf.plt.tight_layout()

        if figpath != None:
            self.bf.plt.savefig(self.bf.os.path.join(figpath, (self.fit_result.file + '_mcmc_walkers.png')), dpi=300)
            self.bf.plt.savefig(self.bf.os.path.join(figpath, (self.fit_result.file + '_mcmc_walkers.pdf')))

        self.bf.plt.show()

    def _find_auto_corr_time(self):
        
        if self.proc_prog == 0:
            print("Walkers must be run before autocorrelation times can be calculated.")
            if self.ErrorHandling == True:
                print("Auto handling error...")

                self._run_walkers()
            else:
                raise RuntimeError("Failed to calculate autocorrelation times as sampler was not run.")

        print("Calculating autocorrelation times...")
        # get autocorrelation time - catch error

        self.fit_result.bayesian_result["tau"]= self.sampler.get_autocorr_time()

        print(self.fit_result.bayesian_result["tau"])

        # catch error for long auto correlation time
        #TBD

    def _flatten_samples(self):
        
        if self.proc_prog == 0:
            print("Walkers must be run before samples can be flattened.")
            if self.ErrorHandling == True:
                print("Auto handling error...")

                self._run_walkers()
            else:
                raise RuntimeError("Failed to flatten samples as sampler was not run.")

        # get autocorrelation time - catch autocorr error
        try:
            self._find_auto_corr_time()
        except self.emcee.autocorr.AutocorrError:
            if self.ErrorHandling == True:
                hold = self.nsteps
                self.nsteps = hold*4
                
                self._run_walkers()
                try:
                    self._find_auto_corr_time()
                except self.emcee.autocorr.AutocorrError:
                    self.nsteps = hold
                    return ValueError("Auto correlation times are too long. Try changing the move method.")
                finally:
                    del hold


        # catch nan    
        if any(self.bf.np.isnan(self.fit_result.bayesian_result["tau"])):
            return ValueError("Autocorrelation times could not be calculated. Make sure that the bounds are set properly and try again.")

        # flatten samples
        print("Flattening samples...")
        self.flat_samples = self.sampler.get_chain(discard=2* int(self.bf.np.max(self.fit_result.bayesian_result["tau"])), thin=15, flat=True)

        # update progress var
        self.proc_prog = 2

    def _gen_corner_plot(self, figpath):

        # catch error if _flatten_samples has not been run
        if self.proc_prog != 2:
            print("Cannot generate corner plots without flattened sample array.")
            if self.ErrorHandling == True:
                print("Auto handling error...")
                self._flatten_samples()
            else:
                raise RuntimeError("Failed to generate corner plot. Check that the sampler has been run and that the samples have been flattened")

        print("Generating corner plots...")
        # plot parameter distributions
        labels = ["Gamma", "tau0", "f0", "baseline"]
        fig2 = self.corner.corner(
            self.flat_samples, labels=labels, truths=[
                    self.fit_result.result['leastsq'].best_values['Gamma'], 
                    self.fit_result.result['leastsq'].best_values['tau0'], 
                    self.fit_result.result['leastsq'].best_values['f0'], 
                    self.fit_result.result['leastsq'].best_values['baseline']
                    ])
        if figpath != None:
            fig2.savefig(self.bf.os.path.join(figpath, (self.fit_result.file + '_mcmc_corner.png')), dpi=300)
            fig2.savefig(self.bf.os.path.join(figpath, (self.fit_result.file + '_mcmc_corner.pdf')))

    def _credible_interval_68(self):
        # catch error if _flatten_samples has not been run
        if self.proc_prog != 2:
            print("Cannot calculate credible interval without flattened sample array.")
            if self.ErrorHandling == True:
                print("Auto handling error...")
                self._flatten_samples()
            else:
                raise RuntimeError("Failed to calculate credible interval. Check that the sampler has been run and that the samples have been flattened")
            
        print(f"Calculating 68% credible intervals...")
        # Gamma
        CI_68_Gamma = [self.bf.np.percentile(self.flat_samples[:,0], 16), self.bf.np.percentile(self.flat_samples[:,0], 84)]
        print(f'Gamma [N s/m]: [{CI_68_Gamma[0]:.3e}, {CI_68_Gamma[1]:.3e}]')
        
        # tau0
        CI_68_tau0 = [self.bf.np.percentile(self.flat_samples[:,1], 16), self.bf.np.percentile(self.flat_samples[:,1], 84)]
        print(f'tau0 [ms]: [{CI_68_tau0[0]:.3e}, {CI_68_tau0[1]:.3e}]')

        # f0
        CI_68_f0 = [self.bf.np.percentile(self.flat_samples[:,2], 16), self.bf.np.percentile(self.flat_samples[:,2], 84)]
        print(f'f0 [kHz]: [{CI_68_f0[0]:.4e}, {CI_68_f0[1]:.4e}]')
        
        # baseline
        CI_68_baseline = [self.bf.np.percentile(self.flat_samples[:,3], 16), self.bf.np.percentile(self.flat_samples[:,3], 84)]
        print(f'baseline [nm^2/Hz]: [{CI_68_baseline[0]:.3e}, {CI_68_baseline[1]:.3e}]')
        
        #store in dictionary
        self.fit_result.bayesian_result["CI_68per_Gamma"] = CI_68_Gamma
        self.fit_result.bayesian_result["CI_68per_tau0"] = CI_68_tau0
        self.fit_result.bayesian_result["CI_68per_f0"] = CI_68_f0
        self.fit_result.bayesian_result["CI_68per_baseline"] = CI_68_baseline

    def _credible_interval_95(self):
        # catch error if _flatten_samples has not been run
        if self.proc_prog != 2:
            print("Cannot calculate credible interval without flattened sample array.")
            if self.ErrorHandling == True:
                print("Auto handling error...")
                self._flatten_samples()
            else:
                raise RuntimeError("Failed to calculate credible interval. Check that the sampler has been run and that the samples have been flattened")
            
        print(f"Calculating 95% credible intervals...")
        # Gamma
        CI_95_Gamma = [self.bf.np.percentile(self.flat_samples[:,0], 2.5), self.bf.np.percentile(self.flat_samples[:,0], 97.5)]
        print(f'Gamma [N s/m]: [{CI_95_Gamma[0]:.3e}, {CI_95_Gamma[1]:.3e}]')
        
        # tau0
        CI_95_tau0 = [self.bf.np.percentile(self.flat_samples[:,1], 2.5), self.bf.np.percentile(self.flat_samples[:,1], 97.5)]
        print(f'tau0 [ms]: [{CI_95_tau0[0]:.3e}, {CI_95_tau0[1]:.3e}]')

        # f0
        CI_95_f0 = [self.bf.np.percentile(self.flat_samples[:,2], 2.5), self.bf.np.percentile(self.flat_samples[:,2], 97.5)]
        print(f'f0 [kHz]: [{CI_95_f0[0]:.4e}, {CI_95_f0[1]:.4e}]')
        
        # baseline
        CI_95_baseline = [self.bf.np.percentile(self.flat_samples[:,3], 2.5), self.bf.np.percentile(self.flat_samples[:,3], 97.5)]
        print(f'baseline [nm^2/Hz]: [{CI_95_baseline[0]:.3e}, {CI_95_baseline[1]:.3e}]')
        
        #store in dictionary
        self.fit_result.bayesian_result["CI_95per_Gamma"] = CI_95_Gamma
        self.fit_result.bayesian_result["CI_95per_tau0"] = CI_95_tau0
        self.fit_result.bayesian_result["CI_95per_f0"] = CI_95_f0
        self.fit_result.bayesian_result["CI_95per_baseline"] = CI_95_baseline

    def _credible_interval_n(self, n:int):
        # catch error if _flatten_samples has not been run
        if self.proc_prog != 2:
            print("Cannot calculate credible interval without flattened sample array.")
            if self.ErrorHandling == True:
                print("Auto handling error...")
                self._flatten_samples()
            else:
                raise RuntimeError("Failed to calculate credible interval. Check that the sampler has been run and that the samples have been flattened")

        if n <=0 or n >=100:
            raise ValueError("n must be between 0 and 100")

        if type(n) != int:
            if type(n) == float:
                hold = n
                n = int(hold)
                print("n must be an integer. Truncating "+str(hold)+" to "+str(n)+"%.")
                del hold
            else:
                raise ValueError("n must be an integer")

        x=50-(n/2)

        print(f"Calculating {n}% credible intervals...")
        # Gamma
        self.CI_n_Gamma = [self.bf.np.percentile(self.flat_samples[:,0], x), self.bf.np.percentile(self.flat_samples[:,0], 100-x)]
        print(f'Gamma [N s/m]: [{self.CI_n_Gamma[0]:.3e}, {self.CI_n_Gamma[1]:.3e}]')
        
        # tau0
        self.CI_n_tau0 = [self.bf.np.percentile(self.flat_samples[:,1], x), self.bf.np.percentile(self.flat_samples[:,1], 100-x)]
        print(f'tau0 [ms]: [{self.CI_n_tau0[0]:.3e}, {self.CI_n_tau0[1]:.3e}]')

        # f0
        self.CI_n_f0 = [self.bf.np.percentile(self.flat_samples[:,2], x), self.bf.np.percentile(self.flat_samples[:,2], 100-x)]
        print(f'f0 [kHz]: [{self.CI_n_f0[0]:.4e}, {self.CI_n_f0[1]:.4e}]')
        
        # baseline
        self.CI_n_baseline = [self.bf.np.percentile(self.flat_samples[:,3], x), self.bf.np.percentile(self.flat_samples[:,3], 100-x)]
        print(f'baseline [nm^2/Hz]: [{self.CI_n_baseline[0]:.3e}, {self.CI_n_baseline[1]:.3e}]')
        
        #store in dictionary
        self.fit_result.bayesian_result["CI_"+str(n)+"per_Gamma"] = self.CI_n_Gamma
        self.fit_result.bayesian_result["CI_"+str(n)+"per_tau0"] = self.CI_n_tau0
        self.fit_result.bayesian_result["CI_"+str(n)+"per_f0"] = self.CI_n_f0
        self.fit_result.bayesian_result["CI_"+str(n)+"per_baseline"] = self.CI_n_baseline

    def run(self, param_bounds, walkers = 64, nsteps = 2000, progress = True, moves = emcee.moves.KDEMove(bw_method="silverman"), figpath = None, n=None):
        """
        Docstring for run
        
        :param param_bounds: Description
        :param walkers:  Description
        :param nsteps: Description
        :param progress: Description
        :param moves: Description
        :param figpath: Description
        :param n: Description
        """
        
        # save variables
        self.param_bounds = param_bounds
        self.walkers = walkers
        self.nsteps = nsteps
        self.progress = progress
        self.moves = moves
        
        # run sampler
        self._run_walkers()

        # plot walkers
        self._plot_walkers(figpath)

        self._flatten_samples()
        
        # plot parameter distributions
        self._gen_corner_plot(figpath)

        # calculate credible intervals
        if n == None:
            return
        
        if n == 68:
            self._credible_interval_68()
        elif n == 95:
            self._credible_interval_95()
        else:
            self._credible_interval_n(n)

        print("Done.")

        


