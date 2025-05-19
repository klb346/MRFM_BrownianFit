###################################################
# LabVIEW_integration.py
# Defines the data class LVprocessing
# Author: Katrina L. Brown
# 2025/05/04
###################################################

class LVprocessing():

    #import libraries
    from MRFM_BrownainFit.brownian_fit import brownian_fit


    def __intit__(self, N_avg:int, temp:float, plot:bool, report:bool, x:list, y:list):
        
        self.N_avg = N_avg
        self.temp = temp
        self.plot = plot
        self.report = plot
        self.x = x
        self.y = y

    def _compile_for_fitting(self):
        self.datatuple = tuple(self.N_avg, self.temp, self.plot, self.report, self.x, self.y)

    def _call_fitting_class(self):
        self._compile_for_fitting
        self.fit = self.brownian_fit(self.datatuple)

    def plot_report(self):
        if self.plot == True:
            if self.report == True:
                self.fit.do_fit_and_plot_with_report()
            else:
                self.fit.do_fit_and_plot()
        else:
            if self.report == True:
                self.fit.do_fit_with_report()
            else:
                self.fit.do_fit()
