from MRFM_BrownianFit.LabVIEW_int import LVprocessing
import h5py
import os
from MRFM_BrownianFit.MCMC import MCMC

# # read in test data
# path = r".\\"
# os.chdir(path)

# # read h5 file for average power spectrum

# file = h5py.File(r'brownian_k_20250924_145412_.h5', 'r')

# temp = float(file['temp'][()])
# n_avgs = int(file['n_avgs'][()])
# freq = file['x'][:]
# power = file['y'][:]

# # check that data is read as expected
# print( len(freq), freq.shape)
# print(temp, n_avgs)


# #test LVprocessing
# data = LVprocessing(N_avgs, temp, x, y, name="test", path = "tests")
# data.do_fit()
# data.fit._extract_peak()
# data.fit._fit_power_spec()
# data.fit._find_params()


def test_data_import():
    try:
        file = h5py.File(r'brownian_k_20250924_145412_.h5', 'r')
    except:
        assert False, ".h5 file does not exist"

def test_data_read_1():
    try:
        temp = float(file['temp'][()])
    except:
        assert False, "file has no attribute 'temp'"

def test_data_read_2():
    try:
        n_avgs = int(file['n_avgs'][()])
    except:
        assert False, "file has no attribute 'n_avgs'"
    
def test_data_read_3():
    try:
        freq = file['x'][:]
        assert freq.shape()==(2000,), "frequency list is not a 1D array with shape (2000,)"
    except:
        assert False, "file has no attribute 'x'"
    
def test_data_read_4():
    try:
        power = file['y'][:]
    except:
        assert False, "file has no attribute 'x'"



def test_initialize_LV_class():
    try:
        # call LVprocessing class
        data = LVprocessing(n_avgs, temp, list(freq), list(power), name="test", path=r".\\example_outputs")
    except:
        assert False, "Failed to initialize LV processing class"

def test_initialize_brownian_class():
    try:
        data







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