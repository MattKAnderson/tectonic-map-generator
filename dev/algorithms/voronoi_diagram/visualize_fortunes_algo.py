import json
import matplotlib.pyplot as plt


xmin = 0.0
ymin = 0.0
xmax = 4096.0
ymax = 4096.0


def plot_line_segments(ax, lines, color=None):
    for line in lines:
        xs = [p[0] for p in line]
        ys = [p[1] for p in line]
        if color:
            ax.plot(xs, ys, color=color)
        else:
            ax.plot(xs, ys)


def plot_points(ax, points, color): 
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.scatter(xs, ys, color=color)


def plot_bounds(ax, xmin, xmax, ymin, ymax, color="r"):
    xs = [xmin, xmin, xmax, xmax, xmin]
    ys = [ymin, ymax, ymax, ymin, ymin]
    ax.plot(xs, ys, color=color)


def plot_region(ax, region, bounds_color="r", centroid_color="b"):
    xs = [p[0] for p in region["vertices"]] + [region["vertices"][0][0]]
    ys = [p[1] for p in region["vertices"]] + [region["vertices"][0][1]]
    ax.plot(xs, ys, color=bounds_color)
    ax.scatter(region["centroid"][0], region["centroid"][1], color=centroid_color)


fig, ax = plt.subplots(1)
base_dir = "build/dev/algorithms/voronoi_diagram/"
#f1 = open("vertex_line_segments.json")
#f2 = open("diagram_seeds.json")
f1 = open(base_dir + "vertex_line_segments.json")
f2 = open(base_dir + "diagram_seeds.json")
f3 = open(base_dir + "region_adjacency_line_segments.json")
f4 = open(base_dir + "regions.json")
line_segments = json.load(f1)["vector"]
region_connections = json.load(f3)["vector"]
diagram_seeds = json.load(f2)["vector"]
regions = json.load(f4)["vector"]

plot_line_segments(ax, line_segments, "r")
plot_line_segments(ax, region_connections, "b")
#plot_points(ax, diagram_seeds, "b")
#for i in range(int(len(regions) / 8)):
#    plot_region(ax, regions[i])
#plot_bounds(ax, xmin, xmax, ymin, ymax)
#plot_points(ax, manual_points, "r")
#for p in sorted(diagram_seeds, key=lambda x: x[0]):
#    print(p)
# print(len(diagram_seeds))

ax.set_xlim(-100.0, 4196.0)
ax.set_ylim(-100.0, 4196.0)
plt.show()
