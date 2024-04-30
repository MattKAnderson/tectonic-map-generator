import json
import matplotlib.pyplot as plt


xmin = 0.0
ymin = 0.0
xmax = 4096.0
ymax = 4096.0


def plot_line_segments(ax, lines, color):
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


fig, ax = plt.subplots(1)
f1 = open("build/dev/algorithms/voronoi_diagram/vertex_line_segments.json")
f2 = open("build/dev/algorithms/voronoi_diagram/diagram_seeds.json")
f3 = open("build/dev/algorithms/voronoi_diagram/region_adjacency_line_segments.json")
line_segments = json.load(f1)["vector"]
region_connections = json.load(f3)["vector"]
diagram_seeds = json.load(f2)["vector"]



plot_line_segments(ax, line_segments, "r")
#plot_line_segments(ax, region_connections, "b")
plot_points(ax, diagram_seeds, "b")
plot_bounds(ax, xmin, xmax, ymin, ymax)
#plot_points(ax, manual_points, "r")
#for p in sorted(diagram_seeds, key=lambda x: x[0]):
#    print(p)
# print(len(diagram_seeds))

ax.set_xlim(-200.0, 4296.0)
ax.set_ylim(-200.0, 4296.0)
plt.show()
