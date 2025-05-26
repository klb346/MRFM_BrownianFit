###################################################
# LabVIEW_int.py
# Defines the data class LVprocessing
# Author: Katrina L. Brown
# 2025/05/20
###################################################

class LVprocessing():

    #import libraries
    from MRFM_BrownianFit.brownian_fit import brownian_fit
    from pylatex import (
        Alignat,
        Axis,
        Document,
        Figure,
        Math,
        Matrix,
        NoEscape,
        Plot,
        Section,
        Subsection,
        Tabular,
        TikZ,
        )
    from datetime import datetime

    def __init__(self, N_avg:int, temp:float, x:list, y:list, name:str, path:str):
        
        self.N_avg = N_avg
        self.temp = temp
        self.name = name
        self.x = x
        self.y = y
        self.path = path
        self.save = str(path +'\\'+ name)

    def _compile_for_fitting(self):
        self.datatuple = tuple((self.N_avg, self.temp, self.x, self.y, self.name))

    def _call_fitting_class(self):
        self._compile_for_fitting()
        self.fit = self.brownian_fit(self.datatuple)

    def do_fit(self):
        """
        The do_fit function will fit on the cantilever peak and store the resulting parameters of
        interest as self.k, self.Q, self.f0, self.force_noise, and self.detector_noise. The full
        fit report can be found from self.result['brownian'].
        The cantilever peak is plotted with the fit function, and the resisuals are plotted below.
        """
        self._call_fitting_class()

        self.fit._extract_peak()
        self.fit._fit_power_spec()
        self.fit._find_params()

    def plot_report(self):
        self._call_fitting_class()

        self.fit._extract_peak()
        self.fit._fit_power_spec()
        self.fit._find_params()

        ##build the report using pylatex
        doc = self.Document(self.save, geometry_options = {"right": "2cm", "left": "2cm"})
        
        #add time stamp to top of doc
        today = self.datetime.now().replace(microsecond=0)
        fd = today.strftime("%Y-%m-%d %H:%M:%S")
        doc.append(str(fd))
        
        #section header - Fit Report
        with doc.create(self.Section("Brownian Motion Fit Report")):
            with doc.create(self.Subsection("Summary")):
                #the math
                doc.append(self.NoEscape("The position-fluctuation power spectrum was fit to the equation:"))
                doc.append("\n")
                doc.append(self.NoEscape("$P(f) = \\frac{kT \\tau_0^2}{\Gamma}\\frac{1}{(\pi\\tau_0)^4(f_0^2-f^2)^2+(\pi\\tau_0)^2f^2} + S_{x}$"))
                doc.append("\n\n")
                doc.append(self.NoEscape("From the optimized fit paramters $\\tau_{0}$, $f_{0}$, $/Gamma$, and $S_{x}$:"))
                doc.append("\n")
                doc.append(self.NoEscape("$Q = \pi  \\tau_0  f_0$"))
                doc.append("\n")
                doc.append(self.NoEscape("$k = 2  \pi^2  f_0^2  \\tau_0  \Gamma$"))
                doc.append("\n")
                #doc.append(self.NoEscape("$P(0) = frac{kt}{\Gamma Q^4}"))
                doc.append(self.NoEscape("$S_{th} = \\frac{kT \\tau_0^2}{(\pi  \\tau_0)^4  f_0^4}$"))
                doc.append("\n")
                doc.append(self.NoEscape("$S_{det} = \\frac{S_x}{z_{rms}^2}$ where $z_{rms} = 0.05  pm$"))
                doc.append("\n")

                #number averages
                doc.append("N_avgs = " + str(self.fit.N_avgs))
                doc.append("\n")
                #temperature
                doc.append("Temperature = " + str(self.fit.temp) + " K")
                doc.append("\n")

                #freq range for fit
                doc.append("Power Spectrum fit to frequency range " + str(self.fit.x_trunc[0]) + " to " + str(self.fit.x_trunc[-1]) + " Hz")
                doc.append("\n\n")
                
                #the fit parameters
                doc.append("From curve fitting the position-fluctuation power spectrum:")
                doc.append("\n")
                #tau_0 with error - ringdown time
                doc.append(self.NoEscape("$\\tau_0 =$ " + str(round(self.fit.tau0,4)) + "$\pm$" + str(round(self.fit.tau0_stderr,4)) + " s"))
                doc.append("\n")
                #f0 with error - resonance freq
                doc.append(self.NoEscape("$f_0 =$ " + str(round(self.fit.f0,2)) + "$\pm$" + str(round(self.fit.f0_stderr, 2)) + " Hz"))
                doc.append("\n")
                #gamma with error - friction coefficent
                doc.append(self.NoEscape("$\Gamma =$ " + str(round(self.fit.result['brownian'].best_values['Gamma'] * 1E12,5)) + "$\pm$" + str(round(self.fit.result['brownian'].params['Gamma'].stderr * 1E12,5)) + " pN s/m"))
                doc.append("\n")
                #Sx with error -(detector noise floor)
                doc.append(self.NoEscape("$S_x =$ " + str(round(self.fit.Sx * 1E6,5)) + "$\pm$" + str(round(self.fit.Sx_stderr * 1E6, 5)) + " ${\mu m}^2$/Hz"))
                doc.append("\n")

                #Q with error
                doc.append(self.NoEscape("$Q =$ " + str(round(self.fit.Q,2)) + "$\pm$" + str(round(self.fit.Q_stderr,2))))
                doc.append("\n")
                #Px(0) with error - zero frequency
                doc.append(self.NoEscape("$P_0 =$ " + str(self.fit.P0) + "$\pm$" + str(self.fit.P0_stderr) + " $nm^2$/Hz"))
                doc.append("\n\n")

                doc.append("From analysis of the area under the position-fluctuation power spectrum:")
                doc.append("\n")
                #k with error
                doc.append(self.NoEscape("$k =$ " + str(round(self.fit.k,5)) + "$\pm$" + str(round(self.fit.k_stderr,5)) + " mN/m"))
                doc.append("\n")
                # intrinsic dissipation
                doc.append((self.NoEscape("Intrinsic dissapation = " + str(self.fit.intrinsic_diss) + "$\pm$" + str(self.fit.intrinsic_diss_stderr) + " pN s/m")))
                doc.append("\n")
                #intrinsic force noise
                doc.append((self.NoEscape("Intrinsic force noise = " + str(self.fit.intrinsic_force_noise) + "$\pm$" + str(self.fit.intrinsic_force_noise_stderr)  + " $aN^2$/Hz")))
                doc.append("\n\n")

                doc.append("Predicted frequency noise:")
                doc.append("\n")
                #zrms
                doc.append((self.NoEscape("$z_{rms} =$ " + str(50) + "nm")))
                doc.append("\n")
                #Sth
                doc.append(self.NoEscape("$S_{th} =$ " + str(self.fit.Sth) + "$\pm$" + str(self.fit.Sth_stderr) + " $Hz^2$/Hz"))
                doc.append("\n")
                #Sdet
                doc.append(self.NoEscape("$S_{det} =$ " + str(self.fit.Sdet) + "$\pm$" + str(self.fit.Sdet_stderr) + " Hz"))
                doc.append("\n")

            #the fit report
            with doc.create(self.Subsection("LMFit Full Report")):
                doc.append(self.NoEscape(self.fit.result['brownian'].fit_report()))

        #section header - Plot fit w/residuals
        with doc.create(self.Section("Brownian Motion Fit Plot")):
            with doc.create(self.Figure(position="htbp")) as plot:
                self.fit.plot_fit(figpath = self.path)
                plot.add_image(filename = self.path+"\\"+self.name+"_residual_cdf.png",width = self.NoEscape(r"1\textwidth"))

        #section header - Plot the CDF
        with doc.create(self.Section("Brownian Residuals Cumulative Distribution Function Plot")):
            doc.append("Residual average is" + str(self.fit.resid_mean))
            with doc.create(self.Figure(position="htbp")) as plot:
                self.fit.residuals_CDF(figpath = self.path)
                plot.add_image(filename = self.path+"\\"+self.name+".png",width = self.NoEscape(r"1\textwidth"))
        
        doc.generate_tex(self.save)
