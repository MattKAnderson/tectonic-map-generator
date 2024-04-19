import math
import matplotlib.pyplot as plt
import numpy as np


def triangle_circumcenter(p1, p2, p3):
    a = p1[0] - p2[0]
    b = p1[1] - p2[1]
    c = p2[0] - p3[0]
    d = p2[1] - p3[1]
    u = 0.5 * (p1[0] * p1[0] + p1[1] * p1[1] - p2[0] * p2[0] - p2[1] * p2[1])
    v = 0.5 * (p2[0] * p2[0] + p2[1] * p2[1] - p3[0] * p3[0] - p3[1] * p3[1])

    inv_denominator = 1.0 / (a * d - b * c);
    x = (d * u - b * v)  * inv_denominator
    y = (a * v - c * u) * inv_denominator

    return [x, y]


def euclidean_distance(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    return math.sqrt(dx * dx + dy * dy)


def quadratic_eqn(a, b, c):
    inner = b * b - 4 * a * c
    if inner < 0:
        return []
    elif inner == 0:
        return [-b / (2 * a)]
    else:
        root_inner = math.sqrt(inner)
        return [
            (-b + math.sqrt(inner)) / (2 * a),
            (-b - math.sqrt(inner)) / (2 * a)
        ]


def plot_circle(ax, origin, radius, xstep=0.001):
    xmin = origin[0] - radius * 1.33
    xmax = origin[0] + radius * 1.33

    xs = list(np.arange(xmin, xmax, xstep))
    pts = []
    for x in xs:
        b = -2 * origin[1]
        c = origin[1]**2 + (x - origin[0])**2 - radius**2
        ys = quadratic_eqn(1.0, b, c)
        for y in ys:
            pts.append([x, y])

    ax.scatter(
        [pt[0] for pt in pts],
        [pt[1] for pt in pts],
        s=5
    )        
    

def plot_point(ax, pt, color):
    ax.scatter(pt[0], pt[1], color=color, marker="X", s=100)


p1 = [12, 18]
p2 = [8, 10]
p3 = [14, 13]
circumcenter = triangle_circumcenter(p1, p2, p3)
radius = euclidean_distance(circumcenter, p1)

ax = plt.gca()
plot_point(ax, p1, "r")
plot_point(ax, p2, "r")
plot_point(ax, p3, "r")
plot_point(ax, circumcenter, "g")
plot_circle(ax, circumcenter, radius)

plt.show()