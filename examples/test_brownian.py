import brownian_fit_class.brownian_fit_class
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Model
import brownian_fit_class.brownian_fit_class as bf


data = bf.brownian_fit("example_brownian_031425.csv")

data.plot_peak()
data.do_fit()