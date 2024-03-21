import json
import random
import matplotlib.pyplot as plt


def plot_border(ax, xsize, ysize):
    xs = [0, xsize, xsize, 0, 0]
    ys = [0, 0, ysize, ysize, 0]
    ax.plot(xs, ys, c="black")

fig, ax = plt.subplots(2, 2)

f = open("build/dev/algorithms/space_partitioning/region_partitions.json", "r")
data = json.load(f)["matrix"]
xsize = len(data)
ysize = len(data[0])

ax[0][0].matshow(data)

f = open("build/dev/algorithms/space_partitioning/border_map.json", "r")
data = json.load(f)["matrix"]

ax[0][1].matshow(data)

f = open("build/dev/algorithms/space_partitioning/traversal_history.json", "r")
traversal_data = json.load(f)["vector"]

print(f"Number of nodes touched: {len(traversal_data)}")

traversal_split_lines = [[traversal_data[0]]]
for i in range(1, len(traversal_data)):
    x1 = traversal_data[i-1][0]
    y1 = traversal_data[i-1][1]
    x2 = traversal_data[i][0]
    y2 = traversal_data[i][1]
    if (abs(x1 - x2) > 1 or abs(y1 - y2) > 1):
        traversal_split_lines.append([traversal_data[i]])
    else:
        traversal_split_lines[-1].append(traversal_data[i])

for line in traversal_split_lines:
    xs = [e[1] for e in line]
    ys = [e[0] for e in line]
    ax[1][0].plot(xs, ys, c="red")

plot_border(ax[1][0], xsize, ysize)

f = open("build/dev/algorithms/space_partitioning/borderlines.json", "r")
data = json.load(f)["matrix"]

print(f"num lines: {len(data)}")

cmap = plt.colormaps["hsv"]
for i, line in enumerate(data):
    print(f"Vertices in line: {len(line)}")
    cc = cmap(random.random())
    split_lines = [[line[0]]]
    for p in line[1:]:
        x1 = p[0]
        y1 = p[1]
        x2 = split_lines[-1][-1][0]
        y2 = split_lines[-1][-1][1]
        if (abs(x1 - x2) > 1 or abs(y1 - y2) > 1):
            split_lines.append([p])
        else:
            split_lines[-1].append(p)

    for sl in split_lines:
        xs = [e[1] for e in sl]
        ys = [e[0] for e in sl]
        ax[1][1].plot(xs, ys, c=cc)
    # plt.text(line[-1][1], line[-1][0], f"line {i}")

ax[1][1].invert_yaxis()
ax[1][0].invert_yaxis()
#plot_border(ax[1][1], xsize, ysize)
# plt.legend()
plt.show()

# TODO: implement logic so that touching lines all have different colors