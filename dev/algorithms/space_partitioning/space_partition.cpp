#include "space_partition.hpp"

using namespace space_partition;

std::vector<Impl::GridCell> Impl::adjacent_ids(
    GridCell& cell, size_t size_x, size_t size_y
) {   
    int above = cell.y - 1 < 0 ? size_y - 1 : cell.y - 1;
    int below = cell.y + 1 < size_y ? cell.y + 1 : 0;
    int left = cell.x - 1 < 0 ? size_x - 1 : cell.x - 1;
    int right = cell.x + 1 < size_x ? cell.x + 1 : 0;

    return {
        {left, above}, {cell.x, above}, {right, above}, {left, cell.y}, 
        {right, cell.y}, {left, below}, {cell.x, below}, {left, below}
    };    
}

std::vector<Impl::GridCellAdjacencyCounter> Impl::unassigned_adjacent_grid_cells(
    GridCell& cell, std::vector<std::vector<int>>& region_affiliation
) {
    int affiliated_plate = region_affiliation[cell.x][cell.y];
    int size_x = region_affiliation.size();
    int size_y = region_affiliation[0].size();

    std::vector<GridCellAdjacencyCounter> unassigned_adjacent;
    auto adjacent_cells = adjacent_ids(cell, size_x, size_y); 
    for (GridCell& c : adjacent_cells) {
        if (region_affiliation[c.x][c.y] == -1) {
            auto adj_adj_cells = adjacent_ids(c, size_x, size_y);
            int adj_count = 0;
            for (GridCell& aac : adj_adj_cells) {
                adj_count += (
                    region_affiliation[aac.x][aac.y] == affiliated_plate
                );
            }
            unassigned_adjacent.emplace_back(c, adj_count);
        }
    }

    return unassigned_adjacent;
}

int Impl::find_parent(std::vector<int>& parent_labels, int label) {
    if (parent_labels[label] == label) {
        return label;
    }
    else {
        int parent = find_parent(parent_labels, parent_labels[label]);
        parent_labels[label] = parent;
        return parent;
    }
}

void Impl::join_labels(
    std::vector<int>& parent_labels, int label_a, int label_b
) {
    int parent_a = find_parent(parent_labels, label_a);
    int parent_b = find_parent(parent_labels, label_b);
    parent_labels[label_a] = parent_b;
}

std::vector<int> Impl::relabel(std::vector<int>& parent_labels) {
    std::unordered_map<int, int> parent_map;
    std::vector<int> relabelled(parent_labels.size());
    for (int i = 0;  i < parent_labels.size(); i++) {   
        int parent = Impl::find_parent(parent_labels, i);
        auto search = parent_map.find(parent);
        if (search == parent_map.end()) {
            int new_label = parent_map.size() + 1;
            relabelled[i] = new_label;
            parent_map[parent] = new_label;
        }
        else {
            relabelled[i] = search->second;
        }
    }
    return relabelled;
}