import math
import matplotlib.pyplot as plt
import numpy as np


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
    

def is_between(a, x, y):
    return min(x, y) < a < max(x, y)


def parabola_x_from_y(directrix, focus, y):
    t1 = ((y - focus[1]) ** 2) / (2 * (focus[0] - directrix))
    t2 = (focus[0] + directrix) / 2
    return t1 + t2


def parabola_intersection(directrix, focus_1, focus_2):
    xf1, yf1 = focus_1
    xf2, yf2 = focus_2
    xd = directrix

    a = xf2 - xf1
    b = 2 * (yf2 * (xf1 - xd) - yf1 * (xf2 - xd))
    c = (yf1 ** 2) * (xf2 - xd) - (yf2 ** 2) * (xf1 - xd) \
      - (xf2 - xf1) * (xf1 - xd) * (xf2 - xd)
    
    intersects = quadratic_eqn(a, b, c)
    print(intersects)
    if len(intersects) == 2 and is_between(intersects[1], yf1, yf2):
        return intersects[1]
    else:
        return intersects[1]


def generate_parabola(focus, directrix, xs):
    lower_points = []
    upper_points = []
    middle = []
    a = 1
    b = - 2 * focus[1]
    partial_c = focus[0] ** 2 + focus[1] ** 2 - directrix ** 2
    for x in xs:
        c = partial_c + 2 * x * (directrix - focus[0])
        ys = quadratic_eqn(a, b, c)
        if len(ys) == 1:
            middle = [[x, ys[0]]]
        elif len(ys) == 2:
            lower_points.append([x, min(ys)])
            upper_points.append([x, max(ys)])
    
    upper_points.reverse()
    return lower_points + middle + upper_points


def plot_parabola(ax, focus, directrix, xmin, xmax, step=0.001):
    xs = list(np.arange(xmin, xmax, step))
    parabola_pts = generate_parabola(focus, directrix, xs)
    pxs = [p[0] for p in parabola_pts]
    pys = [p[1] for p in parabola_pts]
    ax.plot(pxs, pys)


directrix = 20.0
focus_1 = [14.0, 7.0]
focus_2 = [11.0, 10.0]
xmin = -10.0
xmax = 20.0

ax = plt.gca()
plot_parabola(ax, focus_1, directrix, xmin, xmax)
plot_parabola(ax, focus_2, directrix, xmin, xmax)

intersection_y = parabola_intersection(directrix, focus_1, focus_2)
ax.scatter(
    parabola_x_from_y(directrix, focus_1, intersection_y), 
    intersection_y,
    color = "r"
)

plt.show()