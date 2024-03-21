import json
import random
import matplotlib.pyplot as plt


f = open("build/dev/algorithms/boundary_identification/plate_partitioning.json", "r")
data = json.load(f)["matrix"]
xsize = len(data)
ysize = len(data[0])


fig, ax = plt.subplots(1, 2)
ax[0].matshow(data)

# f2 = open("build/dev/algorithms/plate_boundaries.json", "r")
# data = json.load(f2)["matrix"]

# groups = dict()

# print(f"data len: {len(data)}")

# for i in range(len(data)):
#     for k in range(len(data[i])):
#         if data[i][k] == 0:
#             continue
#         if i % 2:
#             p1 = ((i-1)/2+1, k)
#             p2 = ((i-1)/2+1, k+1)
#         else:
#             p1 = (i/2, k+1)
#             p2 = (i/2+1, k+1)
#         print(f"p1: {p1}, p2: {p2}")
#         if data[i][k] in groups:
#             groups[data[i][k]].append([p1, p2])
#         else:
#             groups[data[i][k]] = [[p1, p2]]

# print(len(groups))

# cmap = plt.cm.get_cmap("hsv", len(groups))
# for i, (k, v) in enumerate(groups.items()):
#     for l in v:
#         ax[1].plot([l[0][1], l[1][1]], [l[0][0], l[1][0]], c=cmap(i), label=f"l({k})")

# f3 = open("build/dev/algorithms/boundary_edges.json", "r")
# data = json.load(f3)["edges"]
# print(f"Number of edges: {len(data)}")

# for edge in data:
#     ax[1].plot([edge[0][1], edge[1][1]], [edge[0][0], edge[1][0]], color="red")

f4 = open("build/dev/algorithms/boundary_identification/labelled_boundary_edges.json", "r")
data = json.load(f4)["labelled_edges"]

groups = {}
for d in data:
    if d["label"] in groups:
        groups[d["label"]].append(d["vertices"])
    else:
        groups[d["label"]] = [d["vertices"]]

print(f"Number of distinct labels: {len(groups)}")

# cmap = plt.cm.get_cmap("hsv", len(groups))
cmap = plt.colormaps["hsv"]
for i, (label, group) in enumerate(groups.items()):
    cc = cmap(random.random())
    for segment in group:
        xs = [vertex[1] for vertex in segment]
        ys = [vertex[0] for vertex in segment]
        ax[1].plot(xs, ys, c=cc)


ax[1].set_xlim([0, xsize])
ax[1].set_ylim([ysize, 0])
plt.show()