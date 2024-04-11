import json
import matplotlib.pyplot as plt


def plot_border(ax, xsize, ysize):
    xs = [0, xsize, xsize, 0, 0]
    ys = [0, 0, ysize, ysize, 0]
    ax.plot(xs, ys, c="black")


def plot_region_map(ax, region_data, title):
    ax.matshow(region_data)
    ax.set_title(title)


fig, ax = plt.subplots(2, 2)
f1 = open("build/dev/algorithms/voronoi_diagram/initial_diagram.json", "r")
f2 = open("build/dev/algorithms/voronoi_diagram/iterated_diagram.json", "r")
f3 = open("build/dev/algorithms/voronoi_diagram/iterated_x3_diagram.json", "r")
f4 = open("build/dev/algorithms/voronoi_diagram/iterated_x20_diagram.json", "r")
plot_region_map(ax[0][0], json.load(f1)["matrix"], "Initial Diagram")
plot_region_map(ax[0][1], json.load(f2)["matrix"], "Iterated x1 Diagram")
plot_region_map(ax[1][0], json.load(f3)["matrix"], "Iterated x3 Diagram")
plot_region_map(ax[1][1], json.load(f4)["matrix"], "Iterated x20 Diagram")
plt.show()
