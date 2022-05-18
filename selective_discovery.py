import pm4py
import numpy as np
from copy import deepcopy
from threading import Thread
import networkx as nx
import xlsxwriter
import sys
import os
import csv
import time

from pm4py.algo.evaluation.replay_fitness import algorithm as replay_fitness_evaluator
from pm4py.algo.evaluation.precision import algorithm as precision_evaluator
from pm4py.algo.conformance.alignments.dfg import algorithm as dfg_alignment
from pm4py.algo.filtering.dfg import dfg_filtering
from pm4py.algo.discovery.dfg import algorithm as dfg_discovery

from pm4py.objects.conversion.dfg import converter as df_converter
from pm4py.objects.dfg.utils import dfg_utils
from pm4py.objects.dfg.exporter import exporter as dfg_exporter
from pm4py.objects.dfg.importer import importer as dfg_importer

from pm4py.evaluation.replay_fitness.variants import alignment_based

from pm4py.visualization.dfg import visualizer as dfg_visualization

import dfg_util_functions as utils

parameters = {dfg_visualization.Variants.FREQUENCY.value.Parameters.FORMAT: "jpg"}

report_file_name = ''


def create_report(log_name, output_path):
    header = ['Fitness ', 'Precision ', 'ActualFM ', 'PossibleFM ', '#Removed edges']
    global report_file_name
    report_file_name = output_path + '/selective_' + log_name[:-10] + '.csv'
    with open(report_file_name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
    return report_file_name


def update_report(fitness, precision, actualFM, possibleFM, len_removed_Edges, removed_edge):
    with open(report_file_name, mode='a+') as report:
        report_writer = csv.writer(report, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        edge = str(removed_edge[0]) + ',' + str(removed_edge[1])
        report_writer.writerow([fitness, precision, actualFM, possibleFM, len_removed_Edges, edge])


def export_visual_dfg(indfg, log, name, output_path):
    gviz = dfg_visualization.apply(indfg, log=log, variant=dfg_visualization.Variants.FREQUENCY, parameters=parameters)
    dfg_visualization.save(gviz, output_path + '/' + name)


def export_dfg_file(dfg, sac, eac, output_path, name):
    export_parameters = {dfg_exporter.Variants.CLASSIC.value.Parameters.START_ACTIVITIES: sac,
                         dfg_exporter.Variants.CLASSIC.value.Parameters.END_ACTIVITIES: eac}
    dfg_exporter.apply(dfg, output_path + '/' + name, variant=dfg_exporter.Variants.CLASSIC,
                       parameters=export_parameters)


def export_excel_result(alignment_result, keys, outname, output_path, itera):
    workbook = xlsxwriter.Workbook(output_path + '/heuristics_' + outname + '.xlsx')
    worksheet = workbook.add_worksheet()
    col_num = 0
    worksheet.write(0, col_num, '')
    worksheet.write_column(1, col_num,
                           ['Fitness ' + str(itera), 'Precision ' + itera, 'ActualFM ' + itera,
                            'PossibleFM ' + itera, ])
    col_num = 1
    for key, value in alignment_result.items():
        col_key = keys[key][0] + ',' + keys[key][1]
        worksheet.write(0, col_num, col_key)
        worksheet.write_column(1, col_num, value)
        col_num += 1

    workbook.close()


def export_final_result(alignment_result, outname, output_path):
    workbook = xlsxwriter.Workbook(output_path + '/heuristics_' + outname + '_final.xlsx')
    worksheet = workbook.add_worksheet()
    col_num = 0
    worksheet.write(0, col_num, 'Iteration')
    worksheet.write(1, col_num, 'Removed Edge')
    worksheet.write_column(2, col_num, ['Fitness', 'Precision', 'Actual FM ', 'Possible FM ', '#removed edges'])
    col_num = 1
    for key, value in alignment_result.items():
        worksheet.write(0, col_num, key)
        col_key = value[0][0] + ',' + value[0][1]
        worksheet.write(1, col_num, col_key)
        worksheet.write_column(2, col_num, value[1])
        col_num += 1

    workbook.close()


def calc_fitness(indfg, insa, inea, log, result, index):
    aligned_traces = pm4py.conformance_diagnostics_alignments(log, indfg, insa, inea)
    fitness_dfg = alignment_based.evaluate(aligned_traces)
    result[index] = fitness_dfg["averageFitness"]
    return fitness_dfg["averageFitness"]


def cal_precision(indfg, insa, inea, log, result, index):
    prec = dfg_utils.etconformance_precision(indfg, insa, log)
    result[index] = prec
    return prec


def calculate_alignment(indfg, log):
    petriNet, iniMarking, finalMarking = df_converter.apply(indfg)
    fitness = replay_fitness_evaluator.apply(log, petriNet, iniMarking, finalMarking,
                                             variant=replay_fitness_evaluator.Variants.ALIGNMENT_BASED)
    prec = precision_evaluator.apply(log, petriNet, iniMarking, finalMarking,
                                     variant=precision_evaluator.Variants.ALIGN_ETCONFORMANCE)
    actFM = (2 * fitness['averageFitness'] * prec) / (fitness['averageFitness'] + prec)
    possFM = (2 * fitness['averageFitness']) / (fitness['averageFitness'] + 1.0)
    return fitness['averageFitness'], prec, actFM, possFM


def calculate_dfg_alignment(indfg, insa, inea, log):
    aligned_traces = pm4py.conformance_diagnostics_alignments(log, indfg, insa, inea)
    fitness_dfg = alignment_based.evaluate(aligned_traces)
    prec = dfg_utils.etconformance_precision(indfg, insa, log)
    actFM = (2 * fitness_dfg["averageFitness"] * prec) / (fitness_dfg["averageFitness"] + prec)
    possFM = (2 * fitness_dfg["averageFitness"]) / (fitness_dfg["averageFitness"] + 1.0)
    return fitness_dfg["averageFitness"], prec, actFM, possFM


def prepare(bdfg, in_sa, in_ea, in_act_count):
    all_edges = sorted(bdfg.items(), key=lambda x: x[1], reverse=False)
    indfg = deepcopy(bdfg)
    # if len(all_edges) > 50:
    low_freq_edges = list(filter(lambda x: x[1] < 5, all_edges))
    ingraph, instart_node, inend_node = dfg_filtering.generate_nx_graph_from_dfg(indfg, in_sa, in_ea,
                                                                                 in_act_count)
    for item in low_freq_edges:
        edge = item[0]
        freq = item[1]
        if freq < 5:
            new_dfg, new_graph, untouchable_edge, in_act_count = utils.remove_edge(indfg, ingraph, in_sa, in_ea,
                                                                                   in_act_count, edge, [])
            if not new_dfg:
                continue
            else:
                indfg = new_dfg
                ingraph = new_graph

    high_freq_edges = list(filter(lambda x: x[1] > 4, all_edges))
    return indfg, in_sa, in_ea, in_act_count


def timber(bdfg, in_dfg, insa, inea, in_act_count, log, untouchable_edge, no_gain_edges, itera=0):
    print('Debug: Iteration :', itera)
    i = -1
    # edges = indfg.keys() - untouchable_edge
    DFGs_results = {}
    alignment_result = {}
    alignment_result_edge_name = {}
    max_act_fm = 0
    max_edge_act_fm = ('null', 'null')

    all_edges = sorted(in_dfg.items(), key=lambda x: x[1], reverse=True)
    graph, gstart_node, gend_node = dfg_filtering.generate_nx_graph_from_dfg(in_dfg, insa, inea, in_act_count)
    indfg = deepcopy(in_dfg)
    inactcount = deepcopy(in_act_count)

    threads_fitness = [None] * len(indfg.keys())
    threads_precision = [None] * len(indfg.keys())
    fitness_results = [None] * len(indfg.keys())
    precision_results = [None] * len(indfg.keys())
    all_dfgs = [None] * len(indfg.keys())

    for item in all_edges:
        i = i + 1
        edge = item[0]
        # edge = edge_tuple[0]
        if edge in untouchable_edge:
            continue
        if edge in no_gain_edges:
            print("Debug: NO Gain Edge")
            continue

        # If source/end of edge is already removed.
        if edge[0] not in in_act_count or edge[1] not in in_act_count:
            continue
        else:
            new_dfg, new_graph, untouchable_edge, activities_count = utils.remove_edge(indfg, graph, insa, inea,
                                                                                       inactcount, edge,
                                                                                       untouchable_edge)
            # If there is no changes! It was not a connected graph.
            if not new_dfg:
                continue
            else:
                dfg = new_dfg

            all_dfgs[i] = edge, dfg, activities_count

            threads_fitness[i] = Thread(target=calc_fitness, args=(dfg, insa, inea, log, fitness_results, i))
            threads_fitness[i].start()

            threads_precision[i] = Thread(target=cal_precision, args=(dfg, insa, inea, log, precision_results, i))
            threads_precision[i].start()
            # fitness, precision, act_fm, poss_fm = calculate_dfg_alignment(dfg, insa, inea, log)

    for j in range(i):
        if threads_fitness[j] is not None:
            threads_fitness[j].join()
            threads_precision[j].join()

    for j in range(i):
        if threads_fitness[j] is None:
            continue

        # fitness, prec, actFM, possFM = results[j]
        fitness = fitness_results[j]
        precision = precision_results[j]
        edge, dfg, activities_count = all_dfgs[j]
        act_fm = (2 * fitness * precision) / (fitness + precision)

        if precision == 1.0:
            no_gain_edges.append(edge)
            continue

        poss_fm = (2 * fitness) / (fitness + 1.0)
        if poss_fm < max_act_fm:
            no_gain_edges.append(edge)
            continue

        DFGs_results.update({edge: [dfg, insa, inea, activities_count]})
        removed_edges = set(bdfg.keys()).difference(set(dfg.keys()))
        alignment_result_edge_name.update({edge: [fitness, precision, act_fm, poss_fm, len(removed_edges)]})
        alignment_result.update({i: [fitness, precision, act_fm, poss_fm, len(removed_edges)]})

        if act_fm > max_act_fm:
            max_act_fm = act_fm
            max_edge_act_fm = edge

    return DFGs_results, alignment_result, max_act_fm, max_edge_act_fm, untouchable_edge, no_gain_edges, alignment_result_edge_name


def apply_heuristics(in_dfg, in_sa, in_ea, in_act_count, log, log_name, output_path):
    print('Debug: Apply Selective Heuristics!')
    max_edges = {}
    best_result = {}
    alignment_results = {}
    dfgs_result = {}
    untouchable_edges = []
    no_gain_edges = []
    maxActFM = 0
    maxEdgeActFM = 0
    DFG_iter_result = {0: 0}
    j = 0
    bdfg = deepcopy(in_dfg)
    try:
        while DFG_iter_result:

            DFG_iter_result, alignment_iter_result, maxActFM_curr, maxEdgeActFM_curr, untouchable_edges, no_gain_edges, \
            alignment_result_edge_name = timber(bdfg, in_dfg, in_sa, in_ea, in_act_count, log, untouchable_edges,
                                                no_gain_edges, j)

            if not DFG_iter_result:
                # name = str(j) + "_dfg_" + str(maxEdgeActFM) + ".jpg"
                # export_visual_dfg(DFG_iter_result.get(maxEdgeActFM)[0], log, name, output_path)
                # print("Prime DFG Alignment: ", alignment_result_edge_name.get(maxEdgeActFM))
                # export_dfg_file(DFG_iter_result.get(maxEdgeActFM)[0], in_sa, in_ea, output_path, 'prime_dfg')
                return best_result

            maxEdgeActFM = maxEdgeActFM_curr
            maxActFM = maxActFM_curr
            alignment_results.update({j: alignment_iter_result})
            max_edges.update({j: [maxActFM, maxEdgeActFM]})
            best_result.update({j: [maxEdgeActFM, alignment_result_edge_name.get(maxEdgeActFM)]})
            dfgs_result.update({j: DFG_iter_result.get(maxEdgeActFM)})
            fitness, prec, actFM, possFM, removed_edges = alignment_result_edge_name.get(maxEdgeActFM)
            update_report(fitness, prec, actFM, possFM, removed_edges, maxEdgeActFM)

            # export_excel_result(alignment_iter_result, list(in_dfg.keys()), log_name + '[' + str(j) + ']',
            # output_path, str(j))

            fitness = alignment_result_edge_name.get(maxEdgeActFM)[0]
            precision = alignment_result_edge_name.get(maxEdgeActFM)[1]
            act_FM = alignment_result_edge_name.get(maxEdgeActFM)[2]
            poss_FM = alignment_result_edge_name.get(maxEdgeActFM)[3]
            if precision == 1 or poss_FM < act_FM:
                print("Prime DFG Alignment: ", alignment_result_edge_name.get(maxEdgeActFM))
                export_dfg_file(DFG_iter_result.get(maxEdgeActFM)[0], in_sa, in_ea, output_path, 'prime_dfg')
                return best_result

            # if j > 1:
            # prev_fitness, prev_prec, prev_act_fm, prev_poss = alignment_results.get(j - 1)

            j = j + 1
            in_dfg = deepcopy(DFG_iter_result.get(maxEdgeActFM)[0])
            in_sa = deepcopy(DFG_iter_result.get(maxEdgeActFM)[1])
            in_ea = deepcopy(DFG_iter_result.get(maxEdgeActFM)[2])
            in_act_count = deepcopy(DFG_iter_result.get(maxEdgeActFM)[3])
    except Exception as e:
        print("Prime DFG Alignment: ", alignment_result_edge_name.get(maxEdgeActFM))
        export_final_result(best_result, log_name, output_path)
        export_dfg_file(DFG_iter_result.get(maxEdgeActFM)[0], DFG_iter_result.get(maxEdgeActFM)[1],
                        DFG_iter_result.get(maxEdgeActFM)[2], output_path, 'last_dfg')
        raise
    else:
        pass
    finally:
        pass

    return best_result


def main(logname, dfgName=None):
    # Import Log
    parent_path = '/'
    file_path = "event_logs/" + logname
    log = pm4py.read_xes(file_path)
    base_activities_count = pm4py.get_event_attribute_values(log, "concept:name")

    # Discover DFG
    bdfg, bsa, bea = pm4py.discover_directly_follows_graph(log)
    # pm4py.view_dfg(bdfg, bsa, bea)

    if dfgName:
        print("dfg_Name: ", dfgName)
        bdfg, bsa, bea = dfg_importer.apply("dfgs/" + dfgName, variant=dfg_importer.Variants.CLASSIC, parameters=None)

    output_path = 'f_measure_output/' + logname + '_selective'

    try:
        os.makedirs(output_path, exist_ok=True)
        print("Directory '%s' created successfully" % output_path)
    except OSError as error:
        print("Directory '%s' can not be created" % output_path)

    # export_visual_dfg(bdfg, log, logname + '_base_Dfg', output_path)

    best_result = {}
    # Compute Alignment of base Model
    fitness, precision, act_fm, poss_fm = calculate_dfg_alignment(bdfg, bsa, bea, log)
    print("Base DFG alignment: ", fitness, precision, act_fm, poss_fm)
    best_result.update({-1: [['', ''], [fitness, precision, act_fm, poss_fm]]})
    start_time = time.time()
    create_report(logname, output_path)
    result = apply_heuristics(bdfg, bsa, bea, base_activities_count, log, logname, output_path)
    print("--- %s seconds ---" % (time.time() - start_time))
    best_result.update(result)
    export_final_result(best_result, logname, output_path)


if __name__ == "__main__":
    print(len(sys.argv))
    print(sys.argv)
    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        main(sys.argv[1])
