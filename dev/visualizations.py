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


def plot_border(ax, xsize, ysize):
    xs = [0, xsize, xsize, 0, 0]
    ys = [0, 0, ysize, ysize, 0]
    ax.plot(xs, ys, c="black")


def plot_region_map(ax, region_data, title):
    ax.matshow(region_data)
    ax.set_title(title)


fig, ax = plt.subplots(1, 2)

f1 = open("build/dev/poisson_voronoi_diagram.json", "r")
f2 = open("build/dev/nested_poisson_voronoi_diagram.json", "r")
data1 = json.load(f1)["matrix"]
data2 = json.load(f2)["matrix"]

plot_region_map(ax[0], data1, "Voronoi from Poisson disks")
plot_region_map(ax[1], data2, "Nested Voronoi")

plt.show()

