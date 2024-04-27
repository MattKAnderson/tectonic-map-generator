import math
import matplotlib.pyplot as plt
import numpy as np

first=3
directrix = 2270.001
xmin = -100
xmax = directrix + 50.0
focii = sorted([
[268, 54],
[1879, 787],
[2220, 263]
])


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


def plot_parabola(ax, focus, directrix, xmin, xmax, color, step=0.1):
    xs = list(np.arange(xmin, xmax, step))
    parabola_pts = generate_parabola(focus, directrix, xs)
    pxs = [p[0] for p in parabola_pts]
    pys = [p[1] for p in parabola_pts]
    ax.plot(pxs, pys)


def plot_point(ax, point, color):
    ax.scatter(point[0], point[1], color=color)


def plot_horizontal(ax, y, xmin, xmax, color):
    ax.plot([xmin, xmax], [y, y], color=color)


def main():
    fig, ax = plt.subplots(1)
    #plot_horizontal(ax, 245, xmin, xmax, "r")
    for focus in focii[:first]:
        plot_parabola(ax, focus, directrix, xmin, xmax, "b")
        plot_point(ax, focus, "g")
    plt.show()

main()