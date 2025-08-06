import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from MRFM_BrownianFit.LabVIEW_int import LVprocessing
import pandas as pd
import numpy as np

#read in test data
df=pd.read_csv("tests/example_brownian_031425.csv")

x = list(df['Frequency [Hz] - Plot 0']*1E-3)
y = list(df['nm^2/Hz - Plot 0']*1e6)

temp = float(df['temp [K]'][0])
N_avgs = int(df['averages'][0])


#test LVprocessing
data = LVprocessing(N_avgs, temp, x, y, name="test", path = "tests")
data.do_fit()
data.fit._extract_peak()
data.fit._fit_power_spec()
data.fit._find_params()

def test_Q_example(data):
    """
    checks that fit approximates Q close to the expected value 6937.27 for the example brownian fit
    """
    assert np.allclose(data.Q, float(6937.27), atol = 5e2), "Q is not close to expected value for example data"

def test_f0_example(data):
    """
    checks that fit approximates f0 close to the expected value 8981.79 Hz for the example brownian fit
    """
    assert np.allclose(data.f0, float(8981.79), atol = 5e2), "f0 is not close to expected value for example data"

def test_k_example(data):
    """
    checks that fit approximates k close to the expected value 4.868 mN/m for the example brownian fit
    """
    assert np.allclose(data.k, float(4.868), atol = 1), "k is not close to expected value for example data"

def test_tau0_example(data):
    """
    checks that fit approximates f0 close to the expected value 0.246 s for the example brownian fit
    """
    assert np.allclose(data.tau0*10**-3, float(0.246), atol = 5e-2), "tau0 is not close to expected value for example data"

def test_residuals_reasonable(data):
    """
    Checks that the normalized residuals average to close to zero
    """
    assert np.allclose(np.mean(data.r1), 0, atol=5e-1)