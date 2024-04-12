import math
import json
import matplotlib.pyplot as plt


def plot_centroid(ax, centroid):
    circle = plt.Circle(
        centroid,
        1.25,
        color="b",
        fill=False
    )
    ax.scatter(centroid[0], centroid[1], color="b")
    ax.add_patch(circle)


def plot_sample_points(ax, points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.scatter(xs, ys, color="k")


def plot_cost_history(ax, cost_history):
    xs = [i for i in range(len(cost_history))]
    ys = [p for p in cost_history]
    ax.plot(xs, ys)


def euclidean_centroid(sample_points):
    xsum = sum([p[0] for p in sample_points])
    ysum = sum([p[1] for p in sample_points])
    n = len(sample_points)
    return [int(xsum / n), int(ysum / n)]


def plot_euclidean_centroid(ax, centroid):
    circle = plt.Circle(
        centroid,
        1,
        color="r",
        fill=False
    )
    ax.add_patch(circle)


def compute_cost(point, points):
    cost = 0
    for p in points:
        dx = point[0] - p[0]
        dy = point[1] - p[1]
        cost += dx * dx + dy * dy
    return cost


f1 = open("build/dev/algorithms/polygon_center/visit_history.json")
f2 = open("build/dev/algorithms/polygon_center/cost_history.json")
f3 = open("build/dev/algorithms/polygon_center/sample_points.json")
f4 = open("build/dev/algorithms/polygon_center/centroid.json")

visit_history = json.load(f1)["vector"]
cost_history = json.load(f2)["vector"]
sample_points = json.load(f3)["vector"]
centroid = json.load(f4)["vector"][0]
ecentroid = euclidean_centroid(sample_points)

fig, ax = plt.subplots(1, 2)
plot_centroid(ax[0], centroid)
plot_euclidean_centroid(ax[0], ecentroid)
plot_sample_points(ax[0], sample_points)
plot_cost_history(ax[1], cost_history)
ax[0].set_title("Centroid and sample points")
ax[1].set_title("Cost history")
plt.show()

print(
    "Cost of approx centroid: ", 
    compute_cost(centroid, sample_points)
)
print(
    "Cost of euclidean centroid: ",
    compute_cost(ecentroid, sample_points)
)
