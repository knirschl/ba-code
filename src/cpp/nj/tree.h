//
// Created by knirschl on 02.05.23.
//

#ifndef BA_TREE_H
#define BA_TREE_H

#include <vector>
#include <string>
#include <format>
#include "../misc/DynamicBitset.h"
#include "../misc/meta.h"
#include "../misc/common_types.h"

/**
 * Representation of a tree as a vector.
 */
struct Tree {
    /*
     * Public and private definitions
     */
private:
    struct Node;
public:
    std::vector<Node> tree;
    std::vector<int> leaf_indices{}; // leaf array indices
    std::vector<std::string> idx2leafname{}; // map node array index to locus name ("12_0_1")
    int root;

private:
    /**
     * Internal representation of a node.
     */
    struct Node {
        /**
         * Indices in `tree`
         */
        int idx{-1};
        int left_child_idx{-1};
        int right_child_idx{-1};
        int parent_idx{-1};
        int score{-1};
        bool is_dup{};
        bool is_leaf{}; // can be incorrect
        DynamicBitset covered_groups{};
        dist_t branch_length{-1};
    };

    /**
     * Adds an "empty" node to the tree.
     *
     * @return the index of the new node in the tree vector
     */
    int make_node() {
        int node_idx = tree.size();
        tree.emplace_back();
        return node_idx;
    }

    /**
     * Adds a node to the tree with `left` and `right` as its children.
     *
     * @param left  the left child
     * @param right the right child
     * @return the index of the new node in the tree vector
     */
    int make_node(int left, int right) {
        return make_node(-1, left, -1, right, -1);
    }

    /**
     * Resets the node at index `node_idx` or adds a new one at the next free.
     *
     * @param node_idx  the preferred index of the new node
     * @param left      the left child
     * @param bl_left   the branch length from the left child to the new node
     * @param right     the right child
     * @param bl_right  the branch length from the right child to the new node
     * @return the index of the new node in the tree vector
     */
    int make_node(int node_idx, int left, dist_t bl_left, int right, dist_t bl_right) {
        if (node_idx >= tree.size() || node_idx < 0) {
            node_idx = make_node();
        }
        tree[node_idx].idx = node_idx;
        tree[node_idx].left_child_idx = left;
        tree[node_idx].right_child_idx = right;
        // TODO needed?
        //tree[node_idx].covered_groups = tree[left].covered_groups | tree[right].covered_groups;

        tree[left].parent_idx = node_idx;
        tree[left].branch_length = bl_left;
        tree[right].parent_idx = node_idx;
        tree[right].branch_length = bl_right;

        return node_idx;
    }

    /**
     * Builds a tree from a set of leafs. Overwrites `Tree::tree`.
     *
     * @param leafnames the names of the leafs
     */
    void make_tree(const std::vector<std::string> &leafnames) {
        make_leafs(leafnames);
        int i{};
        while (i + 1 < tree.size()) {
            Tree::root = make_node(i, i + 1);
            i += 2;
        }
    }

public:
    /**
     * Constructs an empty tree.
     */
    Tree() : tree{} {}

    /**
     * Constructs a tree from a set of leafs.
     *
     * @param leafnames the names of the leafs
     */
    Tree(const std::vector<std::string> &leafnames) {
        make_tree(leafnames);
    }

    /**
     * Creates n = leafnames.size() many leaf nodes and set mappings accordingly.
     * Overwrites `Tree::tree`.
     *
     * @param leafnames names of the leafs
     */
    void make_leafs(const std::vector<std::string> &leafnames) {
        // reset structs
        idx2leafname.clear();
        this->leaf_indices.clear();
        idx2leafname.resize(leafnames.size());
        tree = std::vector<Node>{leafnames.size()};
        root = -1;

        for (int leaf_idx{}; leaf_idx < leafnames.size(); leaf_idx++) {
            tree[leaf_idx].idx = leaf_idx;
            tree[leaf_idx].is_leaf = true;
            tree[leaf_idx].score = 0;
            //if (idx2leafname.size() <= leaf_idx) {
            //    idx2leafname.resize(leaf_idx + 1);
            //}
            idx2leafname[leaf_idx] = leafnames.at(leaf_idx);
            this->leaf_indices.emplace_back(leaf_idx);
            // give associated groupname an id or use existing one to set bitset
            std::string groupname = leafname2groupname[leafnames[leaf_idx]];
            if (!(groupname2id.contains(groupname))) {
                groupname2id.emplace(groupname, groupname2id.size());
            }
            tree[leaf_idx].covered_groups.set(groupname2id.at(groupname));
        }
    }

    /**
     * Adds a node to the tree with `left` and `right` as its children.
     *
     * @param left      the left child
     * @param bl_left   the branch length from the left child to the new node
     * @param right     the right child
     * @param bl_right  the branch length from the right child to the new node
     * @return the index of the new node in the tree vector
     */
    int make_node(int left, dist_t bl_left, int right, dist_t bl_right) {
        return make_node(-1, left, bl_left, right, bl_right);
    }

    /**
     * Reverses all branches from `new_parent` upwards until the root. The parent of `new_parent`
     * will be set as the left or right child of `new_child` depending on the `overwrite_left` flag.
     * The parent of the parent of `new_parent` will be set to `new_parent` and the root is skipped.
     *
     * @param new_parent        the new parent (old child)
     * @param overwrite_left    set the old parent as left or right child of the new parent
     */
    void reverse_branch(int new_parent, bool overwrite_left) {
        Node &new_parent_node{tree[new_parent]};
        int old_parent{new_parent_node.parent_idx};
        Node &old_parent_node{tree[old_parent]};
        if (old_parent_node.parent_idx == -1) {
            int old_sibling{
                    old_parent_node.left_child_idx == new_parent ? old_parent_node.right_child_idx
                                                                 : old_parent_node.left_child_idx};
            tree[old_sibling].parent_idx = new_parent;
            tree[old_sibling].branch_length += new_parent_node.branch_length;
            if (overwrite_left) {
                new_parent_node.left_child_idx = old_sibling;
            } else {
                new_parent_node.right_child_idx = old_sibling;
            }
            return;
        }
        if (overwrite_left) {
            new_parent_node.left_child_idx = old_parent;
        } else {
            new_parent_node.right_child_idx = old_parent;
        }
        reverse_branch(old_parent, old_parent_node.left_child_idx == new_parent);
        old_parent_node.parent_idx = new_parent;
        old_parent_node.branch_length = new_parent_node.branch_length;
    }

    /**
     * Inserts a root along the branch (`root_child`, parent(`root_child`)).
     * If `root_child` or its parent is already the root, nothing changes. Else, a new root node is
     * inserted.
     *
     * @param root_child one child of the new root
     * @return the index of the new root in the tree vector
     */
    int reroot(int root_child) {
        int old_parent{tree[root_child].parent_idx};
        if (old_parent == -1) {
            // already root
            return root_child;
        }
        int old_parent_parent{tree[old_parent].parent_idx};
        if (old_parent_parent == -1) {
            // old_parent is root
            return old_parent;
        }
        int new_root{make_node()};
        bool overwrite_left{tree[old_parent].left_child_idx == root_child};
        reverse_branch(old_parent, overwrite_left);
        dist_t new_bl{tree[root_child].branch_length / 2};
        if (overwrite_left) {
            make_node(new_root, root_child, new_bl, old_parent, new_bl);
        } else {
            make_node(new_root, old_parent, new_bl, root_child, new_bl);
        }

        Tree::root = new_root;
        return root;
    }

    /**
     *
     * https://github.com/SchadePaul/MasterThesis/blob/master/njst/src/tree.c
     * https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp
     *
     * @param cur
     * @return
     */
    int tag_APro(int cur) {
        /* TODO
         * if (tree[cur].score != -1) {
            return tree[cur].score;
        }
        */
        if (tree[cur].is_leaf) {
            return tree[cur].score = 0;
        }
        tree[cur].score = tag_APro(tree[cur].left_child_idx) + tag_APro(tree[cur].right_child_idx);
        DynamicBitset &covered_left = tree[tree[cur].left_child_idx].covered_groups;
        DynamicBitset &covered_right = tree[tree[cur].right_child_idx].covered_groups;
        tree[cur].covered_groups = covered_left | covered_right;

        // from astral-pro project
        if (!covered_left.isDisjointTo(covered_right)) {
            tree[cur].is_dup = true;
            tree[cur].score++;
            if (tree[cur].covered_groups != covered_left) {
                tree[cur].score++;
            }
            if (tree[cur].covered_groups != covered_right) {
                tree[cur].score++;
            }
        } else {
            tree[cur].is_dup = false;
        }
        return tree[cur].score;

        // my implementation of pauls thesis
        /*
        if (covered_left.isDisjointTo(covered_right)) {
            // speciation
            tree[cur].is_dup = false;
            return tree[cur].score;
        }
        tree[cur].is_dup = true;
        if (covered_left == tree[cur].covered_groups || covered_right == tree[cur].covered_groups) {
            if (covered_left == covered_right) {
                return tree[cur].score++;
            }
            return tree[cur].score += 2;
        }
        return tree[cur].score += 3;
        */
    }

    /**
     *
     * https://github.com/chaoszhang/ASTER/blob/master/src/astral-pro.cpp
     *
     * @return
     */
    int reroot_APro() {
        int best_score{std::numeric_limits<int>::max()};
        int best_root{-1};
        int best_counter{};

        std::vector<Node> backup_tree{Tree::tree};

        for (int test_root{}; test_root < Tree::tree.size(); test_root++) {
            // TODO Make it work without rerooting the tree all the time
            reroot(test_root);
            int score{tag_APro(test_root)};
            if (score == best_score) {
                best_counter++;
                if (rand() % best_counter == 0) {
                    best_root = test_root;
                }
            } else if (score < best_score) {
                best_score = score;
                best_root = test_root;
                best_counter = 1;
            }
            Tree::tree = backup_tree;
        }

        best_root = reroot(best_root);
        tag_APro(best_root);
        return best_root;
    }

    /**
     * Computes the depth of the node in the tree.
     * @link https://www.geeksforgeeks.org/lowest-common-ancestor-in-a-binary-tree-using-parent-pointer/
     * @param node
     * @return
     */
    int depth(int node) {
        int d{-1};
        while (node >= 0) {
            d++;
            node = tree[node].parent_idx;
        }

        return d;
    }

    /**
     * Easy but not efficient for many queries Lowest Common Ancestor algorithm.
     * Time: O(h)
     * Space: O(1)
     * @link https://www.geeksforgeeks.org/lowest-common-ancestor-in-a-binary-tree-using-parent-pointer/
     *
     * @param low_node
     * @param high_node
     * @return
     */
    int lca(int low_node, int high_node) {
        // find difference of depths of the two nodes
        int depth_difference{depth(low_node) - depth(high_node)};
        // if high_node is deeper, swap low_node and high_node
        if (depth_difference < 0) {
            int tmp{low_node};
            low_node = high_node;
            high_node = tmp;
            depth_difference *= -1;
        }

        // move low_node up until it reaches the same level as high_node
        while (depth_difference--) {
            low_node = tree[low_node].parent_idx;
        }

        while (low_node >= 0 && high_node >= 0) {
            if (low_node == high_node) {
                return low_node;
            }
            low_node = tree[low_node].parent_idx;
            high_node = tree[high_node].parent_idx;
        }

        // no common ancestor
        return -1;
    }

    /**
    * Computes all leaf pairs whose lowest common ancestor was not tagged as a duplication event.
    *
    * @return
    */
    std::vector<std::pair<int, int>> get_speciation_pairs() {
        std::vector<std::pair<int, int>> pairs{};
        for (int i{}; i < this->leaf_indices.size(); i++) {
            int leaf1_id{tree[this->leaf_indices[i]].idx};
            for (int j{i + 1}; j < this->leaf_indices.size(); j++) {
                int leaf2_id{tree[this->leaf_indices[j]].idx};
                if (!tree[lca(leaf1_id, leaf2_id)].is_dup) {
                    pairs.emplace_back(leaf1_id, leaf2_id);
                }
            }
        }

        return pairs;
    }

    std::string get_name_or_idx(int node) const {
        if (node < 0 || node >= idx2leafname.size()) {
            return std::to_string(node);
        }
        return idx2leafname[node];

    }

    std::string to_string(int cur) const {
        if (tree[cur].is_leaf) {
            return idx2leafname.at(cur);
        }
        return std::format("({}:{}, {}:{})",
                           to_string(tree[cur].left_child_idx),
                           tree[tree[cur].left_child_idx].branch_length,
                           to_string(tree[cur].right_child_idx),
                           tree[tree[cur].right_child_idx].branch_length);
    }

    std::string to_newick() const {
        return std::format("{};", to_string(root));
    }

    std::string node_info(int cur) const {
        if (tree[cur].is_leaf) {
            return std::format("LEAF {}[par= {}, bl= {}, dup= {}]", idx2leafname.at(cur),
                               get_name_or_idx(tree[cur].parent_idx), tree[cur].branch_length,
                               tree[cur].is_dup);
        }
        return std::format("NODE {}[par= {}, l= {}, r= {}, bl= {}, dup= {}];    {};    {}",
                           get_name_or_idx(cur), get_name_or_idx(tree[cur].parent_idx),
                           get_name_or_idx(tree[cur].left_child_idx),
                           get_name_or_idx(tree[cur].right_child_idx),
                           tree[cur].branch_length, tree[cur].is_dup,
                           node_info(tree[cur].left_child_idx),
                           node_info(tree[cur].right_child_idx));
    }
};

#endif //BA_TREE_H
