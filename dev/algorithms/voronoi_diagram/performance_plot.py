import matplotlib.pyplot as plt
import numpy as np


nvertices = [
    100,
    250,
    500,
    1000,
    2500,
    5000,
    10000,
    25000,
    50000,
    100000,
]

times = [
    0.3819,
    1.1854,
    1.9937,
    4.0562,
    11.3035,
    23.2941,
    50.043,
    128.465,
    276.036,
    587.705
]

plt.plot(nvertices[:6], times[:6], color="r")
plt.title("Compute time in ms")
plt.show()