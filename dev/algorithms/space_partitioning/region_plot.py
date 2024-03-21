import json
import random
import matplotlib.pyplot as plt

f = open("build/dev/algorithms/space_partitioning/region_partitions.json", "r")
data = json.load(f)["matrix"]
xsize = len(data)
ysize = len(data[0])

plt.matshow(data)
plt.show()