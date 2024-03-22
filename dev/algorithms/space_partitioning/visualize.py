import json
import random
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor


def plot_border(ax, xsize, ysize):
    xs = [0, xsize, xsize, 0, 0]
    ys = [0, 0, ysize, ysize, 0]
    ax.plot(xs, ys, c="black")


def not_adjacent(p1, p2):
    return abs(p1[0] - p2[0]) > 1 or abs(p1[1] - p2[1]) > 1


def split_on_disconnect(split_lines, p1, p2):
    if not_adjacent(p1, p2):
        split_lines.append([p2])
    else:
        split_lines[-1].append(p2)


def separate_travel_history_lines(traversal_data):
    traversal_split_lines = [[traversal_data[0][0]]]
    split_on_disconnect(traversal_split_lines, traversal_split_lines[0][0], traversal_data[0][1])
    for i in range(1, len(traversal_data)):
        split_on_disconnect(traversal_split_lines, traversal_split_lines[-1][-1], traversal_data[i][0])
        split_on_disconnect(traversal_split_lines, traversal_data[i][0], traversal_data[i][1])
    return traversal_split_lines


def plot_edges(ax, traversal_data):
    for edge in traversal_data:
        if not_adjacent(edge[0], edge[1]):
            continue
        xs = [p[1] for p in edge]
        ys = [p[0] for p in edge]
        ax.plot(xs, ys, c="red")


def plot_region_map(ax, region_data):
    ax.matshow(region_data)


# function returns a color for each line, such that at each intersection
# no two lines have the same color. 
def select_colors(split_lines):
    # these could be tuned to amke the plot nicer, but specifically 
    # choosing colors is much better than using a colormap
    css_named_colors = [
        "gray", "lightcoral", "darkred", "chocolate",
        "orange", "tan", "darkkhaki", "yellow", "chartreuse",
        "forestgreen", "aquamarine", "lightseagreen", "darkcyan",
        "royalblue", "navy", "slateblue", "blueviolet", "indigo",
        "violet",  "magenta"
    ]
    ncolor = 20
    # cmap = plt.get_cmap("tab20", 20)
    counter = 0
    interesting_vertices = set()
    vertex_to_lines = {}
    for i in range(len(split_lines)):
        for j in range(len(split_lines[i])):
            p = (split_lines[i][j][0], split_lines[i][j][1])
            if p in vertex_to_lines:
                interesting_vertices.add(p)
                vertex_to_lines[p].append(i)
            else:
                vertex_to_lines[p] = [i]

    line_intersections = [[] for i in range(len(split_lines))]
    for v in interesting_vertices:
        for i in set(vertex_to_lines[v]):
            line_intersections[i].append(v)
    
    intersections = {}
    chosen_color_ids = []
    for i in range(len(split_lines)):
        intersecting_colors = []
        for v in line_intersections[i]:
            if v in intersections:
                intersecting_colors.extend(intersections[v])
            else:
                intersections[v] = []
        color_id = random.randint(0, ncolor-1)
        while (color_id in intersecting_colors):
            counter += 1
            color_id = random.randint(0, ncolor-1)
        chosen_color_ids.append(color_id)
        for v in line_intersections[i]:
            intersections[v].append(color_id)

    print(f"Stopped {counter} instances of color clashing")

    # return [
    #    cmap(cid/ncolor) for cid in chosen_color_ids
    #]
    return [
        mcolor.CSS4_COLORS[css_named_colors[cid]] for cid in chosen_color_ids
    ]


def plot_colored_borderlines(ax, borderlines_data):
    all_split_lines = []
    for i, line in enumerate(borderlines_data):
        split_lines = [[line[0]]]
        for p in line[1:]:
            split_on_disconnect(split_lines, split_lines[-1][-1], p)
        all_split_lines.extend(split_lines)
    colors = select_colors(all_split_lines)
    for sl, c in zip(all_split_lines, colors):
        xs = [e[1] for e in sl]
        ys = [e[0] for e in sl]
        ax.plot(xs, ys, c=c)


f = open("build/dev/algorithms/space_partitioning/region_partitions.json", "r")
region_data = json.load(f)["matrix"]
fig, ax = plt.subplots(1, 3)
xsize = len(region_data)
ysize = len(region_data[0])
plot_region_map(ax[0], region_data)


f = open("build/dev/algorithms/space_partitioning/border_map.json", "r")
# border_map_data = json.load(f)["matrix"]
# ax[0][1].matshow(border_map_data)


f = open("build/dev/algorithms/space_partitioning/traversal_history.json", "r")
traversal_data = json.load(f)["vector"]
print(f"Number of nodes touched: {len(traversal_data)}")
plot_edges(ax[1], traversal_data)
plot_border(ax[1], xsize, ysize)


f = open("build/dev/algorithms/space_partitioning/borderlines.json", "r")
borderlines_data = json.load(f)["matrix"]
print(f"num lines: {len(borderlines_data)}")
plot_colored_borderlines(ax[2], borderlines_data)
plot_border(ax[2], xsize, ysize)

ax[2].invert_yaxis()
ax[1].invert_yaxis()
plt.show()