import json
import matplotlib.pyplot as plt


def plot_polygon(ax, polygon_vertices):
    xs = [p[0] for p in polygon_vertices]
    ys = [p[1] for p in polygon_vertices]

    # add logic to break this up at the boundaries

    ax.plot(xs, ys)


def plot_points(ax, points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    ax.scatter(xs, ys)


def plot_centroid(ax, centroid):
    circle = plt.Circle(
        centroid,
        1.25,
        color="r",
        fill=False
    )
    ax.scatter(centroid[0], centroid[1], color="r")
    ax.add_patch(circle)


def plot_cost_history(ax, cost_history):
    xs = [i for i in range(len(cost_history))]
    ys = [p for p in cost_history]
    ax.plot(xs, ys)


class VectorLoader:
    def __init__(self, base_path):
        self.base_path = base_path

    def load(self, filename):
        return json.load(open(self.base_path + filename))["vector"]



loader = VectorLoader("build/dev/algorithms/polygon_center/")

polygon_vertices = loader.load("example_vertices.json")
polygon_samples = loader.load("example_polygon_samples.json")
centroid = loader.load("example_centroid.json")[0]
cost_history = loader.load("example_cost_history.json")
visit_history = loader.load("example_visit_history.json")

translated_samples = loader.load("example_translated_samples.json")
translated_centroid = loader.load("example_translated_centroid.json")[0]
translated_visit_history = loader.load("example_translated_visit_history.json")
translated_cost_history = loader.load("example_translated_cost_history.json")

fig, ax = plt.subplots(2, 2)

plot_polygon(ax[0][0], polygon_vertices)
plot_points(ax[0][0], polygon_samples)
plot_centroid(ax[0][0], centroid)
plot_cost_history(ax[0][1], cost_history)

plot_points(ax[1][0], translated_samples)
plot_centroid(ax[1][0], translated_centroid)
plot_cost_history(ax[1][1], translated_cost_history)
plt.show()