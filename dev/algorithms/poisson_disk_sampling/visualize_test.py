import json
import matplotlib.pyplot as plt


def plot_points_with_radius(ax, title, points, radius):

    for point in points:
        circle = plt.Circle(
            (point[0], point[1]), 
            radius/2, 
            color="k", 
            fill=False
        )
        ax.add_patch(circle)
        ax.scatter(point[0], point[1], color="b")

    ax.set_title(title)


fig, ax = plt.subplots(1, 2)

f1 = open("build/dev/algorithms/poisson_disk_sampling/unit_fill.json", "r")
f2 = open("build/dev/algorithms/poisson_disk_sampling/grid_fill.json", "r")
data1 = json.load(f1)["vector"]
data2 = json.load(f2)["vector"]

plot_points_with_radius(ax[0], "Unit Fill", data1, 0.15)
plot_points_with_radius(ax[1], "Grid Fill", data2, 128.0)

plt.show()

