import json
import matplotlib.pyplot as plt


def plot_line_segments(ax, lines):
    for line in lines:
        xs = [p[0] for p in line]
        ys = [p[1] for p in line]
        ax.plot(xs, ys)


def plot_points(ax, points, color): 
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.scatter(xs, ys, color=color)


fig, ax = plt.subplots(1)
f1 = open("build/dev/algorithms/voronoi_diagram/line_segments.json")
f2 = open("build/dev/algorithms/voronoi_diagram/diagram_seeds.json")
line_segments = json.load(f1)["vector"]
diagram_seeds = json.load(f2)["vector"]
plot_line_segments(ax, line_segments)
plot_points(ax, diagram_seeds, "b")
#plot_points(ax, manual_points, "r")
#for p in sorted(diagram_seeds, key=lambda x: x[0]):
#    print(p)
print(len(diagram_seeds))

ax.set_xlim(-100.0, 4196.0)
ax.set_ylim(-100.0, 4196.0)
plt.show()
