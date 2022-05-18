import pm4py
import numpy as np
from copy import deepcopy
import networkx as nx
import xlsxwriter
import sys
import os
import time

from pm4py.algo.filtering.dfg import dfg_filtering
from pm4py.algo.evaluation.replay_fitness import algorithm as replay_fitness_evaluator
from pm4py.algo.evaluation.precision import algorithm as precision_evaluator
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery

from pm4py.evaluation.replay_fitness.variants import alignment_based

from pm4py.objects.conversion.dfg import converter as df_converter
from pm4py.objects.dfg.utils import dfg_utils
from pm4py.objects.dfg.exporter import exporter as dfg_exporter
from pm4py.objects.dfg.importer import importer as dfg_importer

from pm4py.visualization.dfg import visualizer as dfg_visualization
import dfg_util_functions as utils


def calculate_dfg_alignment(indfg, insa, inea, log):
    aligned_traces = pm4py.conformance_diagnostics_alignments(log, indfg, insa, inea)
    fitness_dfg = alignment_based.evaluate(aligned_traces)
    prec = dfg_utils.etconformance_precision(indfg, insa, log)
    actFM = (2 * fitness_dfg["averageFitness"] * prec) / (fitness_dfg["averageFitness"] + prec)
    possFM = (2 * fitness_dfg["averageFitness"]) / (fitness_dfg["averageFitness"] + 1.0)
    return fitness_dfg["averageFitness"], prec, actFM, possFM


def export_visual_dfg(indfg, log, name, output_path):
    parameters = {dfg_visualization.Variants.PERFORMANCE.value.Parameters.FORMAT: "jpg"}
    gviz = dfg_visualization.apply(indfg, log=log, variant=dfg_visualization.Variants.FREQUENCY, parameters=parameters)
    dfg_visualization.save(gviz, output_path + '/' + name)


def export_dfg_file(dfg, sac, eac, output_path, name):
    print('Debug: Begin Export DFG!')
    export_parameters = {dfg_exporter.Variants.CLASSIC.value.Parameters.START_ACTIVITIES: sac,
                         dfg_exporter.Variants.CLASSIC.value.Parameters.END_ACTIVITIES: eac}
    dfg_exporter.apply(dfg, output_path + '/' + name, variant=dfg_exporter.Variants.CLASSIC,
                       parameters=export_parameters)
    print('Debug: END Export DFG!')


def export_final_result(alignment_result, outname, output_path):
    workbook = xlsxwriter.Workbook(output_path + '/frequency_' + outname + '_final.xlsx')
    worksheet = workbook.add_worksheet()
    col_num = 0
    worksheet.write(0, col_num, 'Iteration')
    worksheet.write(1, col_num, 'Removed Edge')
    worksheet.write_column(2, col_num, ['Fitness', 'Precision', 'Actual FM ', 'Possible FM ', 'len of removed edges'])
    col_num = 1
    for key, value in alignment_result.items():
        worksheet.write(0, col_num, key)
        col_key = value[0][0] + ',' + value[0][1]
        worksheet.write(1, col_num, col_key)
        worksheet.write_column(2, col_num, value[1])
        col_num += 1

    workbook.close()


def run_iter(indfg, insa, inea, inactcount, log, logname, output_path):
    print('Debug: begin frequency-based discovery, size of dfg:', len(indfg))
    DFGs_results = {}
    alignment_result = {}
    untouchable_edge = []
    max_actual_fm = 0
    edge_max_actual_fm = ('null', 'null')
    iter_max_act_fm = 0
    highest_dfg = deepcopy(indfg)

    graph, instart_node, inend_node = dfg_filtering.generate_nx_graph_from_dfg(indfg, insa, inea, inactcount)
    dfg = deepcopy(indfg)
    start_activities = deepcopy(insa)
    end_activities = deepcopy(inea)
    activities_count = deepcopy(inactcount)

    i = -1
    # Sort edges of DFG based on frequency.
    all_edges = sorted(indfg.items(), key=lambda x: x[1], reverse=False)
    start_time = time.time()
    for item in all_edges:
        edge = item[0]
        i = i + 1
        if edge in untouchable_edge:
            # print('Debug: Untouchable Edge:', edge)
            continue

        # If source/end of edge is already removed.
        if edge[0] not in activities_count:
            continue
        if edge[1] not in activities_count:
            continue
        else:
            try:
                new_dfg, new_graph, untouchable_edge, activities_count = utils.remove_edge(dfg, graph, insa, inea,
                                                                                           inactcount,
                                                                                           edge,
                                                                                           untouchable_edge )
                # If there is no changes! It was not a connected graph.
                if not new_dfg:
                    continue
                else:
                    dfg = new_dfg
                    graph = new_graph
            except Exception as e:
                print("Error: While removing edge: ", edge, ", iteration: ", i)
                print("Highest DFG Alignment: ", alignment_result.get(iter_max_act_fm))
                export_dfg_file(dfg, insa, inea, output_path, 'dfg_error')
                raise

        DFGs_results.update({edge: [dfg, start_activities, end_activities, activities_count]})
        fitness, prec, actFM, possFM = calculate_dfg_alignment(dfg, insa, inea, log)
        if max_actual_fm < actFM:
            max_actual_fm = actFM
            highest_dfg = dfg
            edge_max_actual_fm = edge
            iter_max_act_fm = i
        removed_edges = set(indfg.keys()).difference(set(dfg.keys()))
        alignment_result.update({i: [edge, [fitness, prec, actFM, possFM, len(removed_edges)]]})
        if prec == 1.0:
            break

    print("--- %s seconds ---" % (time.time() - start_time))
    print("Size of Prime DFG: ", len(highest_dfg))
    print("Prime DFG Alignment: ", alignment_result.get(iter_max_act_fm))
    # export_visual_dfg(dfg, log, logname + '_final_Dfg', output_path)
    # export_visual_dfg(DFGs_results.get(edge_max_actual_fm)[0], log,
    #                   logname + '_best_actFM_Dfg_' + str(edge_max_actual_fm[0]) + '-' + str(edge_max_actual_fm[1]),
    #                   output_path)
    export_dfg_file(highest_dfg, insa, inea, output_path, 'prime_dfg')
    return DFGs_results, alignment_result


def main(logname, dfgName=None):
    # Import Log
    parent_path = '/home/yasi/MI/thesis/a_first_heuristics/'
    file_path = parent_path + "event_logs/" + logname
    # file_path = "event_logs/" + logname
    logname = logname.replace('.xes', '')
    print("file_path: ", file_path)
    log = pm4py.read_xes(file_path)
    base_activities_count = pm4py.get_event_attribute_values(log, "concept:name")

    if dfgName:
        print("dfg_Name: ", dfgName)
        bdfg, bsa, bea = dfg_importer.apply("dfgs/" + dfgName, variant=dfg_importer.Variants.CLASSIC, parameters=None)

    output_path = 'complete_frequency_method_output/' + logname + '_freq_output'

    try:
        os.makedirs(output_path, exist_ok=True)
        print("Directory '%s' created successfully" % output_path)
    except OSError as error:
        print("Directory '%s' can not be created" % output_path)

    # Discover DFG
    bdfg, bsa, bea = pm4py.discover_directly_follows_graph(log)
    # pm4py.view_dfg(bdfg, bsa, bea)
    # export_visual_dfg(bdfg, log, logname + '_base_Dfg', output_path)

    best_result = {}
    # Compute Alignment of base Model
    fitness, prec, actFM, possFM = calculate_dfg_alignment(bdfg, bsa, bea, log)
    print("Base DFG alignment: ", fitness, prec, actFM, possFM)
    best_result.update({-1: [['', ''], [fitness, prec, actFM, possFM]]})

    DFGs_results, alignment_result = run_iter(bdfg, bsa, bea, base_activities_count, log, logname, output_path)

    best_result.update(alignment_result)
    export_final_result(best_result, logname, output_path)


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        main(sys.argv[1])
