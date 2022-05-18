import uuid
from copy import deepcopy
from threading import Thread

import pm4py
from pm4py.algo.filtering.dfg import dfg_filtering
from pm4py.evaluation.replay_fitness.variants import alignment_based

import dfg_util_functions as utils
from pm4py.objects.dfg.utils import dfg_utils


def run_iter(bdfg, indfg, insa, inea, inactcount, log, untouchable_edge, no_gain_dfg_edge_sets,
             max_act_fm, edge_max_act_fm, discovered_dfg):
    DFGs_results = {}
    alignment_result = {}
    alignment_result_edge_name = {}
    max_act_fm_iter = 0
    edge_max_act_fm_iter = ('null', 'null')

    # no_gain_dfg_edge_sets = []
    ingraph, instart_node, inend_node = dfg_filtering.generate_nx_graph_from_dfg(indfg, insa, inea, inactcount)

    threads_fitness = [None] * len(indfg.keys())
    threads_precision = [None] * len(indfg.keys())
    fitness_results = [None] * len(indfg.keys())
    precision_results = [None] * len(indfg.keys())
    all_dfgs = [None] * len(indfg.keys())

    i = -1
    for edge in indfg.keys():
        i = i + 1
        if edge in untouchable_edge:
            print('Debug: Untouchable Edge')
            continue

        if is_a_no_gain_edge(bdfg, indfg, edge, no_gain_dfg_edge_sets):
            print('Debug: No gain Edge')
            continue

        # Remove edge and calculate alignment
        dfg, graph, untouchable_edge, activities_count = utils.remove_edge(indfg, ingraph, insa, inea, inactcount, edge,
                                                                           untouchable_edge)
        # If there is no changes!
        if not dfg:
            continue

        # discovered_dfg contains removed edges of all discovered dfg. Here we prevent repetitions.
        removed_edges = set(bdfg.keys()).difference(set(dfg.keys()))
        if is_already_discovered(discovered_dfg, removed_edges):
            continue

        all_dfgs[i] = edge, dfg, activities_count

        threads_fitness[i] = Thread(target=calc_fitness, args=(dfg, insa, inea, log, fitness_results, i))
        threads_fitness[i].start()

        threads_precision[i] = Thread(target=cal_precision, args=(dfg, insa, inea, log, precision_results, i))
        threads_precision[i].start()

        # threads[i] = Thread(target=utils.calculate_dfg_alignment, args=(dfg, insa, inea, log, results, i))
        # threads[i].start()
        # fitness, prec, actFM, possFM = utils.calculate_dfg_alignment(dfg, insa, inea, log)

    for j in range(i):
        if threads_fitness[j] is not None:
            threads_fitness[j].join()
            threads_precision[j].join()

    for j in range(i):
        if threads_fitness[j] is None:
            continue

        # fitness, prec, actFM, possFM = results[j]
        fitness = fitness_results[j]
        prec = precision_results[j]
        edge, dfg, activities_count = all_dfgs[j]
        actFM = (2 * fitness * prec) / (fitness + prec)

        # Update Max Actual FM
        if actFM > max_act_fm:
            max_act_fm = actFM
            edge_max_act_fm = edge

        # Do not Add the new dfg to the list if it's PossFM is less than max_act_fm
        poss_FM = (2 * fitness) / (fitness + 1.0)
        if poss_FM < max_act_fm:
            removed_edges = set(bdfg.keys()).difference(set(indfg.keys()))
            removed_edges.add(edge)
            no_gain_dfg_edge_sets.append(removed_edges)
            continue

        if actFM > max_act_fm_iter:
            max_act_fm_iter = actFM
            edge_max_act_fm_iter = edge

        # Otherwise, add it to the lists
        removed_edges = set(bdfg.keys()).difference(set(dfg.keys()))
        DFGs_results.update({edge: [dfg, insa, inea, activities_count, removed_edges]})
        alignment_result.update({i: [fitness, prec, actFM, poss_FM, len(removed_edges)]})
        alignment_result_edge_name.update({edge: [fitness, prec, actFM, poss_FM, len(removed_edges)]})

    if not alignment_result:
        return DFGs_results, alignment_result, max_act_fm, edge_max_act_fm_iter, untouchable_edge, alignment_result_edge_name
    return DFGs_results, alignment_result, max_act_fm, edge_max_act_fm_iter, untouchable_edge, alignment_result_edge_name


def update_no_gain_dfg_edge_sets(alignment_result_edge_name, DFGs_result, no_gain_dfg_edge_sets, max_act_fm, bdfg):
    alignment_result = deepcopy(alignment_result_edge_name)
    # print('DEBUG: Update no_gain_dfg_edge_sets, max Act FM: ', max_act_fm)
    for item in alignment_result.items():
        edge = item[0]
        possFm = item[1][3]
        if possFm < max_act_fm:
            print('DEBUG: poss_FM < max_ActFM --> edge, possFM:', edge, possFm)
            dfg = DFGs_result[edge][0]
            # removed_edges = set(bdfg.keys()).difference(set(dfg.keys()))
            removed_edges = DFGs_result[edge][4]
            # if not check_edge_sets(no_gain_dfg_edge_sets, removed_edges):
            #     no_gain_dfg_edge_sets.append(removed_edges)
            no_gain_dfg_edge_sets.append(removed_edges)
            del DFGs_result[edge]
            del alignment_result_edge_name[edge]


def check_edge_sets(no_gain_dfg_edge_sets, removed_edges):
    for e in no_gain_dfg_edge_sets:
        if e.issubset(removed_edges):
            return True
    return False


def update_untouchables_edge_node(edge_aligns, DFGs_result, bdfg, bacts, maxActFM, iter_num):
    result = {}
    untouchble_edge = []
    alignment_result = deepcopy(edge_aligns)
    dfgs = deepcopy(DFGs_result)
    untouchble_nodes = set()
    for item in edge_aligns.items():
        edge = item[0]
        possFm = item[1][3]
        if possFm < maxActFM:
            print('DEBUG: poss_FM < max_ActFM ', edge, possFm)
            del alignment_result[edge]
            del dfgs[edge]
            untouchble_edge.append(edge)
            dfg = DFGs_result[edge][0]
            removed_edges = set(bdfg.keys()).difference(set(dfg.keys()))
            act_count = set(DFGs_result[edge][3].keys())
            removed_nodes = set(bacts.keys()).difference(act_count)
            if removed_nodes:
                # IF these nodes are already completely were not in untouchble_nodes. Include the result.
                if not removed_nodes.issubset(untouchble_nodes):
                    result.update({str(uuid.uuid4().fields[-1])[:5]: [removed_nodes, removed_edges, iter_num]})
                #                 any(x in A for x in D)
                untouchble_nodes.update(removed_nodes)
            #                 print('***********')
    return result, dfgs, alignment_result, untouchble_edge, untouchble_nodes


def is_a_no_gain_edge(bdfg, indfg, edge, no_gain_dfg_edge_sets):
    will_removed_edges = set(bdfg.keys()).difference(set(indfg.keys()))
    will_removed_edges.add(edge)
    # .add(edge)
    for item in no_gain_dfg_edge_sets:
        if item.issubset(will_removed_edges):
            return True
        # elif will_removed_edges.issubset(item):
        #     return True
    return False


def is_already_discovered(discovered_dfg, removed_edges):
    # removed_edges = set(bdfg.keys()).difference(set(dfg.keys()))
    if bool(discovered_dfg) and len(removed_edges) in discovered_dfg.keys():
        for item in discovered_dfg[len(removed_edges)]:
            if item == removed_edges:
                return True
        (discovered_dfg[len(removed_edges)]).append(removed_edges)
    else:
        discovered_dfg[len(removed_edges)] = [removed_edges]
    return False


def calc_fitness(indfg, insa, inea, log, result, index):
    aligned_traces = pm4py.conformance_diagnostics_alignments(log, indfg, insa, inea)
    fitness_dfg = alignment_based.evaluate(aligned_traces)
    result[index] = fitness_dfg["averageFitness"]
    return fitness_dfg["averageFitness"]


def cal_precision(indfg, insa, inea, log, result, index):
    prec = dfg_utils.etconformance_precision(indfg, insa, log)
    result[index] = prec
    return prec
