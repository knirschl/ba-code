import os
import sys
from threading import Thread
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/msa')
sys.path.insert(0, 'tools/trees')
import fam
import metrics
import analyze_msa
import rf_distance


def ttask(datadir, tree, avg_abs_dico, avg_rel_dico, families_dico):
    for family in fam.get_families_list(datadir):
        if (not analyze_msa.has_distinct_seqs(
                fam.get_alignment_file(fam.get_family_path(datadir, family)))):
            # not enough distinct sequences
            # !! -> there shouldn't be any trees except if left over from old runs
            continue
        true_tree = fam.get_true_tree(datadir, family)
        abs_path_tree = os.path.join(fam.get_gene_tree_dir(datadir, family), tree)
        if (not os.path.isfile(abs_path_tree)):
            picked_abs_path_tree = abs_path_tree.replace(".geneTree.newick",
                                                         ".generax_pick.geneTree.newick")
            if (not os.path.isfile(picked_abs_path_tree)):
                # print("Tree file does not exist:", family, tree)
                continue
            abs_path_tree = picked_abs_path_tree
            # print("Picked", tree, "  (path", abs_path_tree, ")")
        # CALCULATIONS
        try:
            #dist_abs, dist_rel = rf_distance.raxmlng_rf(abs_path_tree, true_tree)
            dist_abs, dist_rel = rf_distance.ete3_rf_files(abs_path_tree, true_tree)
            #print(dist_abs, dist_rel)
        except Exception as exc:
            # during converting or computing
            print("I've catched an excpetion, maybe look into this?\n", exc)
            print(abs_path_tree, true_tree, dist_abs, dist_rel, "\n")
            continue
        # single distance
        # This should be thread-safe
        families_dico[family][tree] = dist_rel
        # method average
        if not tree in avg_abs_dico:
            # if tree (= method) not in dico, add
            avg_abs_dico[tree] = []
            avg_rel_dico[tree] = []
        avg_abs_dico[tree].append(dist_abs)
        avg_rel_dico[tree].append(dist_rel)

def compare_all(datadir):
    avg_abs_dico = {}
    avg_rel_dico = {}
    families_dico = {}
    for family in fam.get_families_list(datadir):
        families_dico[family] = {}
    threads = []
    possible_trees = set().union(*[fam.get_gene_tree_list(datadir, family) for family in fam.get_families_list(datadir)])
    for tree in possible_trees:
        if (tree == fam.get_true_gene_tree_name()):
            continue
        tree = tree.replace(".generax_pick", "")
        fam_thread = Thread(target=ttask, args=(datadir, tree, avg_abs_dico, avg_rel_dico, families_dico))
        threads.append(fam_thread)
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    # save method average distance
    for tree in avg_abs_dico:
        avg_abs_dico[tree] = sum(avg_abs_dico[tree]) / len(avg_abs_dico[tree])
        avg_rel_dico[tree] = sum(avg_rel_dico[tree]) / len(avg_rel_dico[tree])
    metrics.save_dico(datadir, avg_abs_dico, "rf_distance_avg-abs")
    metrics.save_dico(datadir, avg_rel_dico, "rf_distance_avg-rel")
    for family in fam.get_families_list(datadir):
        metrics.save_dico(fam.get_family_path(datadir, family), families_dico[family], "rf_distance-rel")