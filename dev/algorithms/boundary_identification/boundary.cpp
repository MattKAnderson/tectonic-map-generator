#include "boundary.hpp"
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <utility>

struct LabelNode {
    LabelNode(int value) : value(value), parent(nullptr) {};
    LabelNode(int value, LabelNode* parent) : value(value), parent(parent) {};
    int value;
    LabelNode* parent;
};

LabelNode* find_root_parent(LabelNode* label) {
    if (label->parent == nullptr) {
        return label;
    }
    else {
        LabelNode* root_parent = find_root_parent(label->parent);
        label->parent = root_parent;
        return root_parent;
    }
}

void label_join(std::vector<LabelNode*>& labels, int label_a_val, int label_b_val) {
    std::cout << "joining labels: " << label_a_val << "   and   " << label_b_val << "\n";
    LabelNode* label_a = labels[label_a_val - 1];
    LabelNode* label_b = labels[label_b_val - 1];
    find_root_parent(label_a)->parent = find_root_parent(label_b);
}

void add_boundary_edges(std::vector<std::vector<std::vector<edge*>>>& vertices, std::vector<std::vector<int>>& group_map, int x1, int y1, int x2, int y2) {
    if (group_map[x1][y1] != group_map[x2][y2]) {
        edge* new_edge = new edge({x1, y1}, {x2, y2}, group_map[x1][y1], group_map[x2][y2]);
        vertices[x1][y1].push_back(new_edge);
        vertices[x2][y2].push_back(new_edge);
    }    
}

/*
 *   This one works great for identifying the boundaries, but it does not label the boundaries
 *   based on the regions on either side
 */
std::vector<edge> boundary_outlines(std::vector<std::vector<int>>& group_map) {
    
    int xsize = group_map.size();
    int ysize = group_map[0].size();

    // construct a graph from the vertices with labelled edges. Edges are labelled with the 
    // groups on either side of the boundary. Construct boundary segments by exploring the 
    // produced graph and 'exhausting' edges as they are crossed in a DFS.

    std::vector<edge> edges;

    for (int i = 0; i < xsize - 1; i++) {       
        for (int k = 0; k < ysize - 1; k++) {
            if (group_map[i][k] != group_map[i][k + 1]) {
                edges.push_back(edge({i, k + 1}, {i + 1, k + 1}, group_map[i][k], group_map[i][k + 1]));
            }
            if (group_map[i][k] != group_map[i + 1][k]) {
                edges.push_back(edge({i + 1, k}, {i + 1, k + 1}, group_map[i][k], group_map[i + 1][k]));
            }
        }
        if (group_map[i][0] != group_map[i][ysize - 1]) {
            edges.push_back(edge({i, 0}, {i + 1, 0}, group_map[i][0], group_map[i][ysize - 1]));
        }
    }
    for (int k = 0; k < ysize - 1; k++) {
        if (group_map[xsize - 1][k] != group_map[xsize - 1][k + 1]) {
            edges.push_back(edge({xsize - 1, k + 1}, {xsize, k + 1}, group_map[xsize - 1][k], group_map[xsize - 1][k + 1]));
        }
        if (group_map[0][k] != group_map[xsize - 1][k]) {
            edges.push_back(edge({xsize, k}, {xsize, k + 1}, group_map[0][k], group_map[xsize - 1][k]));
        }
    }
    if (group_map[xsize - 1][0] != group_map[xsize - 1][ysize - 1]) {
        //edges.push_back(edge({xsize - 1, 0}, {xsize, 0}, group_map[xsize - 1][0], group_map[xsize - 1][ysize - 1]));
    }

    return edges; 
}

int find_labelled_edge_parent(std::vector<std::pair<int, edge>>& labelled_edges, int label) {
    if (label == labelled_edges[label].first) {
        return label;
    }
    else {
        int parent = find_labelled_edge_parent(labelled_edges, labelled_edges[label].first);
        labelled_edges[label].first = parent;
        return parent;
    }
}

void join_labelled_edges(std::vector<std::pair<int, edge>>& labelled_edges, int label_a, int label_b) {
    int parent_a = find_labelled_edge_parent(labelled_edges, label_a);
    int parent_b = find_labelled_edge_parent(labelled_edges, label_b);
    labelled_edges[parent_a].first = parent_b;
}

//find a better name in a moment
void label_edges_visit_vertex(
    std::unordered_map<Vector2D<int>, std::vector<std::pair<int, edge*>>, hash_int_Vector2D>& vertices, 
    std::vector<std::pair<int, edge>>& labelled_edges,
    Vector2D<int>& vertex,
    int i
) {
    auto& e = labelled_edges[i].second;
    if (auto vertex_itr = vertices.find(vertex); vertex_itr != vertices.end()) {
        // check if you need to union with any edges and then add this edge the vertex
        auto our_boundary = std::make_pair(std::min(e.group_a, e.group_b), std::max(e.group_a, e.group_b));
        for (auto& other_edge : vertex_itr->second) {
            auto their_boundary = std::make_pair(
                std::min(other_edge.second->group_a, other_edge.second->group_b), 
                std::max(other_edge.second->group_a, other_edge.second->group_b)
            );
            if (our_boundary == their_boundary) {
                join_labelled_edges(labelled_edges, i, other_edge.first);
                break; // if 2 other edges are part of the union then they are already joined on a previous iter
            }
        }
        vertices[vertex].emplace_back(i, &e);
    }
    else {
        vertices[vertex] = {{i, &e}};
    }
}

std::vector<std::pair<int, edge>> label_edges(std::vector<edge>& edges) {
    std::vector<std::pair<int, edge>> labelled_edges;
    labelled_edges.reserve(edges.size());
    for (int i = 0; i < edges.size(); i++) {
        labelled_edges.emplace_back(i, edges[i]);
    }
    
    std::unordered_map<Vector2D<int>, std::vector<std::pair<int, edge*>>, hash_int_Vector2D> vertices;
    for (int i = 0; i < labelled_edges.size(); i++) {
        label_edges_visit_vertex(vertices, labelled_edges, labelled_edges[i].second.p1, i);
        label_edges_visit_vertex(vertices, labelled_edges, labelled_edges[i].second.p2, i);
    }

    // re-label at the end so that the labels increase monotonically
    std::unordered_set<int> parents;
    for (int i = 0; i < labelled_edges.size(); i++) {
        parents.insert(find_labelled_edge_parent(labelled_edges, i));
    }
    std::unordered_map<int, int> parent_to_relabel;
    int next_label = 1;
    for (int p : parents) {
        parent_to_relabel[p] = next_label++;
    }
    for (int i = 0; i < labelled_edges.size(); i++) {
        labelled_edges[i].first = parent_to_relabel[labelled_edges[i].first];
    }

    return labelled_edges;
}


std::vector<std::vector<int>> return_boundaries_2(std::vector<std::vector<int>>& group_map) {

    int xsize = group_map.size();
    int ysize = group_map[0].size();

    // construct a graph from the vertices with labelled edges. Edges are labelled with the 
    // groups on either side of the boundary. Construct boundary segments by exploring the 
    // produced graph and 'exhausting' edges as they are crossed in a DFS.

    std::vector<std::vector<std::vector<edge*>>> vertices;
    std::vector<std::vector<std::vector<Vector2D<int>>>> boundary_lines; 

    for (int i = 0; i < xsize - 1; i++) {       
        for (int k = 0; k < ysize - 1; k++) {
            add_boundary_edges(vertices, group_map, i, k, i, k + 1);
            add_boundary_edges(vertices, group_map, i, k, i + 1, k);
        }
        add_boundary_edges(vertices, group_map, i, 0, i, ysize - 1);
    }
    for (int k = 0; k < ysize - 1; k++) {
        add_boundary_edges(vertices, group_map, xsize - 1, k, xsize - 1, k + 1);
        add_boundary_edges(vertices, group_map, 0, k, xsize - 1, k);
    }
    add_boundary_edges(vertices, group_map, xsize - 1, 0, xsize - 1, ysize - 1);

    /*
     *   Iterate over the vertices, each time a vertex is found that has at least 1 edge do a 
     *   DFS and build a linked list representation of the boundary. This linked list will 
     *   provide a sequence of points defining the line and it can be used to recover the 
     *   boundary grid as well. 
     */
    int next_label = 1;
    std::vector<std::vector<std::pair<int, int>>> lines;

    for (int i = 0; i < xsize; i++) {
        for (int k = 0; k < ysize; k++) {
            if (vertices[i][k].size() == 0) {
                continue;
            }

            for (int j = 0; j < vertices[i][k].size(); j++) {

                // check if same boundary as any of the other edges at this vertex
                int index = -1;
                edge* e1 = vertices[i][k][j];
                for (int l = j - 1; l >= 0; l--) {
                    auto other_pair = std::make_pair(vertices[i][k][l]->group_a, vertices[i][k][l]->group_b);
                    if (std::make_pair(e1->group_a, e1->group_b) == other_pair || std::make_pair(e1->group_b, e1->group_a) == other_pair) {
                        index = lines.size() + j - l;
                    }
                }
                if (index == -1) {
                    index = lines.size();
                    lines.push_back({});
                }

                // do a dfs

            }

            

        }
    }
    return std::vector<std::vector<int>>();
}

std::vector<std::vector<int>> return_boundaries(std::vector<std::vector<int>>& group_map) {

    int xsize = group_map.size();
    int ysize = group_map[0].size();
    
    /*
     *   Data structure used to join labels (a boundary could be given more than 1 label while)
     *   iterating over the group_map i.e. upside down u for example would be read as 2 boundaries
     *   initially, and only at the bottom would it be apparent).  Post processing will re-label 
     *   all boundaries so each has 1 label and label number will increase monotonically
     */
    std::vector<LabelNode*> label_tree;

    /*
     *   boundary map elements are labels of a unique boundary between two groups in group_map. 
     *   Labelling starts at 1, 0 is used to represent no boundary. The boundary between the i 
     *   and i+1 row in group_map is mapped to the 2i row in boundary map (k=k). The boundary between 
     *   the k and k+1 column is mapped to the 2i+1 row (k=k). 
     */
    
    std::vector<std::vector<int>> boundary_map(2*xsize, std::vector<int>(ysize));

    for (int i = 0; i < xsize; i++) {
        for (int k = 0; k < ysize; k++) {
            int below = i + 1 < xsize ? i + 1 : 0;
            int right = k + 1 < ysize ? k + 1 : 0;

            if (group_map[i][k] != group_map[i][right]) {
                int above = i - 1 >= 0 ? i - 1 : xsize - 1;
                if (group_map[i][k] != group_map[above][k] && boundary_map[2 * above + 1][k] != 0) {
                    boundary_map[2 * i][k] = boundary_map[2 * above + 1][k];
                }
                else if (group_map[i][k] != group_map[above][right] && boundary_map[2 * above][k] != 0) {
                    boundary_map[2 * i][k] = boundary_map[2 * above][k];
                }
                else if (boundary_map[2 * above + 1][right] != 0) {
                    boundary_map[2 * i][k] = boundary_map[2 * above + 1][right];
                }
                else {
                    LabelNode* next_label = new LabelNode(label_tree.size() + 1, nullptr);
                    label_tree.push_back(next_label); 
                    boundary_map[2 * i][k] = label_tree.size();
                }
            }

            if (group_map[i][k] != group_map[below][k]) {
                int left = k - 1 >= 0 ? k - 1 : ysize - 1;

                if (group_map[i][k] != group_map[i][left]) {
                    if (group_map[below][k] == group_map[i][left] && boundary_map[2 * i][left] != 0) {
                        boundary_map[2 * i + 1][k] = boundary_map[2 * i][left];
                        if (group_map[i][right] == group_map[i][left]) {
                            label_join(label_tree, boundary_map[2* i][k], boundary_map[2 * i + 1][k]);
                        }
                    }
                    else if (group_map[below][k] == group_map[i][right]) {
                        boundary_map[2 * i + 1][k] = boundary_map[2 * i][k];
                    }
                    else {
                        LabelNode* next_label = new LabelNode(label_tree.size() + 1, nullptr);
                        label_tree.push_back(next_label); 
                        boundary_map[2 * i + 1][k] = label_tree.size();  
                    }
                }
                else if (group_map[i][k] != group_map[below][left] && boundary_map[2 * i + 1][left] != 0) {
                    boundary_map[2 * i + 1][k] = boundary_map[2 * i + 1][left];
                    if (boundary_map[2 * i][k] != 0 && boundary_map[2 * i][k] != boundary_map[2 * i + 1][k] && group_map[below][k] == group_map[i][right]) {
                        label_join(label_tree, boundary_map[2* i][k], boundary_map[2 * i + 1][k]);
                    }                    
                }
                else if (boundary_map[2 * i][k] != 0) {
                    boundary_map[2 * i + 1][k] = boundary_map[2 * i][k];
                }
                else {
                    LabelNode* next_label = new LabelNode(label_tree.size() + 1, nullptr);
                    label_tree.push_back(next_label); 
                    boundary_map[2 * i + 1][k] = label_tree.size();                   
                }
            }
        }
        
    }

    std::vector<int> label_map(label_tree.size());
    for (int i = 0; i < label_map.size(); i++) {
        label_map[i] = find_root_parent(label_tree[i])->value;
    }
    for (int i = 0; i < 2 * xsize; i++) {
        for (int k = 0; k < ysize; k++) {
            boundary_map[i][k] = boundary_map[i][k] == 0 ? 0 : label_map[boundary_map[i][k]];
        }
    }
    

    return boundary_map;
    
}