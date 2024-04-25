import math
import matplotlib.pyplot as plt
import numpy as np

first=128
directrix = 1013.001
xmin = -100
xmax = directrix + 50.0
focii = sorted([
[216, 900],
[257, 288],
[406, 611],
[704, 1005],
[452, 529],
[513, 716],
[351, 120],
[23, 428],
[27, 697],
[811, 591],
[890, 691],
[227, 514],
[190, 734],
[584, 353],
[787, 924],
[118, 396],
[941, 129],
[254, 964],
[175, 407],
[925, 749],
[223, 900],
[470, 499],
[614, 286],
[880, 291],
[801, 22],
[814, 408],
[77, 442],
[102, 379],
[367, 366],
[957, 550],
[763, 357],
[377, 214],
[516, 377],
[460, 101],
[219, 434],
[770, 460],
[146, 694],
[919, 202],
[248, 138],
[886, 887],
[948, 355],
[891, 923],
[923, 171],
[394, 1006],
[837, 963],
[815, 853],
[994, 824],
[882, 848],
[651, 853],
[287, 603],
[334, 226],
[596, 848],
[29, 22],
[782, 584],
[823, 870],
[478, 42],
[783, 979],
[643, 636],
[940, 262],
[80, 406],
[973, 652],
[484, 232],
[575, 198],
[822, 733],
[941, 789],
[26, 647],
[687, 171],
[859, 695],
[707, 402],
[754, 965],
[785, 63],
[92, 152],
[218, 402],
[558, 755],
[129, 73],
[784, 845],
[91, 726],
[164, 516],
[498, 541],
[98, 206],
[589, 264],
[981, 509],
[521, 543],
[727, 709],
[507, 595],
[746, 928],
[846, 451],
[724, 1017],
[229, 977],
[372, 187],
[395, 482],
[677, 692],
[626, 114],
[859, 752],
[1013, 390],
[382, 959],
[0, 882],
[374, 272],
[981, 851],
[719, 176],
[50, 537],
[689, 338],
[232, 161],
[565, 129],
[411, 854],
[397, 371],
[365, 499],
[1, 34],
[190, 279],
[719, 817],
[152, 1002],
[999, 893],
[1019, 994],
[374, 140],
[91, 905],
[573, 854],
[892, 89],
[1, 47],
[599, 79],
[367, 82],
[319, 482],
[380, 472],
[899, 934],
[522, 926],
[80, 990],
[427, 578],
[389, 422],
[23, 43],
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