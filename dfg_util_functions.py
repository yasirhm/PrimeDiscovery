from copy import deepcopy

import networkx as nx
import pm4py
from pm4py.algo.filtering.dfg import dfg_filtering
from pm4py.evaluation.replay_fitness.variants import alignment_based
from pm4py.objects.dfg.utils import dfg_utils


def calculate_dfg_alignment(indfg, insa, inea, log, result, index):
    aligned_traces = pm4py.conformance_diagnostics_alignments(log, indfg, insa, inea)
    fitness_dfg = alignment_based.evaluate(aligned_traces)
    prec = dfg_utils.etconformance_precision(indfg, insa, log)
    actFM = (2 * fitness_dfg["averageFitness"] * prec) / (fitness_dfg["averageFitness"] + prec)
    possFM = (2 * fitness_dfg["averageFitness"]) / (fitness_dfg["averageFitness"] + 1.0)
    result[index] = [fitness_dfg["averageFitness"], prec, actFM, possFM]
    return fitness_dfg["averageFitness"], prec, actFM, possFM


def remove_edge(dfg, graph, insa, inea, activities_count, edge, untouchable_edge):
    if edge[0] == edge[1]:  # if it is a loop
        if edge in dfg:
            new_Dfg = deepcopy(dfg)
            del new_Dfg[edge]
            new_graph, instart_node, inend_node = dfg_filtering.generate_nx_graph_from_dfg(new_Dfg, insa, inea,
                                                                                           activities_count)
            return new_Dfg, new_graph, untouchable_edge, activities_count

    new_graph = nx.DiGraph(graph)
    new_graph.remove_edge(edge[0], edge[1])
    if not nx.is_weakly_connected(new_graph):
        # print('Not Connected')
        untouchable_edge.append(edge)
        return {}, {}, untouchable_edge, activities_count

    # reachable_from_start = set(nx.descendants(new_graph, "START"))
    start = list(insa.items())[0][0]
    end = list(inea.items())[0][0]
    reachable_from_start = set(nx.descendants(new_graph, start))
    # reachable_to_end = set(nx.ancestors(new_graph, "END"))
    reachable_to_end = set(nx.ancestors(new_graph, end))
    reachable_start_end = reachable_from_start.intersection(reachable_to_end)
    reachable_start_end.add(start)   #('START')
    reachable_start_end.add(end) #('END')

    activities_set = set(activities_count.keys())
    non_reachable_activities = activities_set.difference(reachable_start_end)
    new_dfg = deepcopy(dfg)
    if edge in new_dfg:
        del new_dfg[edge]
    # remove these non-reachable activities
    new_act_count = deepcopy(activities_count)
    for act in non_reachable_activities:
        new_dfg = {x: y for x, y in new_dfg.items() if x[0] != act and x[1] != act}
        del new_act_count[act]
    #         if act in start_activities:
    #             del start_activities[act]
    #         if act in end_activities:
    #             del end_activities[act]

    if non_reachable_activities:
        set(new_act_count.keys())

        # make sure that the DFG contains only edges between these activities
    new_dfg = {x: y for x, y in new_dfg.items() if x[0] in new_act_count and x[1] in new_act_count}
    return new_dfg, new_graph, untouchable_edge, new_act_count
