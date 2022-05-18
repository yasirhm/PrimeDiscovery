import sys
import os
import csv
import logging
from time import perf_counter
import pm4py
from pm4py.objects.dfg.importer import importer as dfg_importer

import xlsxwriter
import BFS
import in_out
from pm4py.objects.dfg.exporter import exporter as dfg_exporte

report_file_name = ''


def create_report(log_name):
    header = ['Fitness ', 'Precision ', 'ActualFM ', 'PossibleFM ', '#Removed edges']
    global report_file_name
    report_file_name = 'greedy_' + log_name[:-4] + '.csv'
    logging.info('REPORT NAME:',report_file_name )
    with open(report_file_name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        # write the header
        writer.writerow(header)
    return report_file_name


def update_report(fitness, precision, actualFM, possibleFM, len_removed_Edges, removed_edge):
    with open(report_file_name, mode='a+') as report:
        report_writer = csv.writer(report, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        edge = str(removed_edge[0])+','+str(removed_edge[1])
        report_writer.writerow([fitness, precision, actualFM, possibleFM, len_removed_Edges, edge])


def export_excel_final_result(alignment_result, out_name, output_path, max_act_fm):
    print('DEBUG: Export Begin!')
    workbook = xlsxwriter.Workbook(output_path + '/final_heuristic_' + out_name[:-4] + '.xlsx')
    worksheet = workbook.add_worksheet()
    col_num = 0
    worksheet.write(0, col_num, '')
    worksheet.write_column(1, col_num,
                           ['Fitness ', 'Precision ', 'ActualFM ', 'PossibleFM ', '#Removed edges'])

    try:
        out_dfgs = {}
        col_num = 1
        for key, value in alignment_result.items():
            if not all(key):
                continue
            if not all(value):
                continue
            col_key = str(key[0]).join(' ').join(str(key[1]))  # Write edge
            worksheet.write(0, col_num, col_key)
            worksheet.write_column(1, col_num, value[1])
            col_num += 1
            act_fm = value[1][2]
            if max_act_fm <= act_fm:
                dfg = value[0][0]
                name = '_' + col_key+'_' + str(col_num) + '.jpg'
                out_dfgs.update({name: dfg})
        workbook.close()
        return out_dfgs
    except NameError:
        print("Error: Variable x is not defined!")
    except TypeError:
        print("Error: Variable does not have value!")


def apply_heuristic(bdfg, indfg, insa, inea, inactcount, log, logname, untouchable_edge, best_result,
                    max_act_fm, edge_max_act_fm, no_gain_dfg_edge_sets, discovered_dfg, itera=0):
    print('Debug: Apply Heuristics!', itera)
    untouchable_nodes = []
    edge_max_act_fm_iter = ('null', 'null')

    try:
        dfg_iter_result, alignment_iter_result, max_act_fm, edge_max_act_fm_iter, untouchable_edge, \
        alignment_result_edge_name = BFS.run_iter(bdfg, indfg, insa, inea, inactcount, log, untouchable_edge,\
                                            no_gain_dfg_edge_sets, max_act_fm, edge_max_act_fm, discovered_dfg)

        # Sort the iteration result based on ActFM
        # !!!!! dfg_iter_result is not sorted
        sorted_x = dict(sorted(alignment_result_edge_name.items(), key=lambda kv: kv[1][2], reverse=True))
        alignment_result_edge_name = sorted_x
        # final_dict = dict(zip(alignment_iter_result.keys(), list(sorted_x.values())))
        # alignment_iter_result = final_dict
        BFS.update_no_gain_dfg_edge_sets(alignment_result_edge_name, dfg_iter_result, no_gain_dfg_edge_sets, max_act_fm,
                                         bdfg)

        # Should consider that there are maybe more than one max FM
        if all(edge_max_act_fm_iter) and alignment_result_edge_name:
            best_result.update({(itera, edge_max_act_fm_iter): [dfg_iter_result.get(edge_max_act_fm_iter),
                                                                alignment_result_edge_name.get(edge_max_act_fm_iter)]})
            fitness, prec, actFM, possFM, removed_edges = alignment_result_edge_name.get(edge_max_act_fm_iter)
            update_report(fitness, prec, actFM, possFM, removed_edges, edge_max_act_fm_iter)
            name_dfg = str(edge_max_act_fm_iter)
            # logging.DEBUG('DFG:', dfg_iter_result.get(edge_max_act_fm_iter)[0])
            in_out.export_dfg_file(dfg_iter_result.get(edge_max_act_fm_iter)[0], insa, inea,'greedy_method_output',name_dfg)
        # DFGs_results.update({edge: [dfg, insa, inea, activities_count, removed_edges]})
        # prime_dfg = DFGs_results[]
    except Exception as e:
        logging.info("Error: itera '{0}' :".format(itera))
        logging.error(e)
        raise

    return alignment_result_edge_name, dfg_iter_result, untouchable_edge, max_act_fm


def run_service(log, log_name, bdfg, bsa, bea, base_activities_count, output_path):
    # output_path = "/home/yasi/MI/thesis/a_first_heuristics/new_greedy_result/";
    best_result = {}  # A dictionary of best alignment result of each iteration
    visited = []  # List to keep track of visited nodes.
    queue = []  # Initialize a queue
    # queue.append({ ('','') :[bdfg, bsa, bea, srlg_acts]})
    # e = ('', [bdfg, bsa, bea, base_activities_count])
    e = ([bdfg, bsa, bea, base_activities_count])
    queue.append(e)
    untouchable_edge = []
    itera = 0
    max_act_fm = 0
    edge_max_act_fm = ('null', 'null')
    no_gain_dfg_edge_sets = []
    discovered_dfg = {}
    start_time = perf_counter()
    try:
        while queue:
            _curr = queue.pop(0) #Get the first element of queue
            curr_dfg = _curr[0]   #[1][0]  # list(_curr.values())[0][0]
            curr_act_count = _curr[3]  # list(_curr.values())[0][3]

            curr_precision = _curr[5][1] if (len(_curr) > 4) else 0
            if curr_precision == 1:
                visited.append(curr_dfg)
                continue

            iter_best_aligns, iter_dfgs_result, untouchable_edge, max_act_fm = apply_heuristic(bdfg, curr_dfg, bsa, bea,
                                                                                               curr_act_count, log,
                                                                                               log_name,
                                                                                               untouchable_edge,
                                                                                               best_result,
                                                                                               max_act_fm,
                                                                                               edge_max_act_fm,
                                                                                               no_gain_dfg_edge_sets,
                                                                                               discovered_dfg,
                                                                                               itera)

            visited.append(curr_dfg)
            for e in iter_best_aligns.items():
                # add all the alignment value and dfg to the queue
                alaki = iter_dfgs_result.get(e[0])
                alaki.append(e[1])
                queue.append(alaki)
            itera = itera + 1

        # out_dfgs = export_excel_final_result(best_result, log_name, output_path, max_act_fm)
        # in_out.export_visual_dfg(out_dfgs.items()[1], log, out_dfgs.items()[0], output_path)

        end_time = perf_counter()
        logging.info(f'END--It took {end_time - start_time :0.2f} second(s) to complete.')
    except KeyboardInterrupt:
        logging.debug('Hello KeyboardInterrupt')
        # in_out.export_dfg_file(dfg,bsa, bea, output_path, log_name[:-4]+'_greedy_'+max_act_fm)
        export_excel_final_result(best_result, log_name, output_path, max_act_fm)
    except Exception as e:
        print('Error: Exporting current result!')
        out_dfgs = export_excel_final_result(best_result, log_name+'+Error+', output_path, max_act_fm)
        logging.debug(e)
        raise
    finally:
        # print('Finally: Exporting current result!')
        # export_excel_final_result(best_result, log_name, output_path)

        pass


def start(log_name, dfg_Name=None):
    # Import Log
    parent_path = '/'
    # file_path = "event_logs/" + log_name
    file_path = parent_path + log_name
    log = pm4py.read_xes(file_path)
    base_activities_count = pm4py.get_event_attribute_values(log, "concept:name")

    # Discover DFG
    bdfg, bsa, bea = pm4py.discover_directly_follows_graph(log)
    # pm4py.view_dfg(bdfg, bsa, bea)

    # If want to begin from a DFG
    if dfg_Name:
        print("dfg_Name: ", dfg_Name)
        bdfg, bsa, bea = dfg_importer.apply("dfgs/" + dfg_Name, variant=dfg_importer.Variants.CLASSIC, parameters=None)

    output_path = 'greedy_method_output/' + log_name[:-4] + '_greedy_output'

    try:
        os.makedirs(output_path, exist_ok=True)
        print("Directory '%s' created successfully" % output_path)
    except OSError as error:
        print("Directory '%s' can not be created" % output_path)

    # in_out.export_visual_dfg(bdfg, log, log_name + '_base.jpg', output_path)
    create_report(log_name)
    run_service(log, log_name, bdfg, bsa, bea, base_activities_count, output_path)
    # export_visual_dfg(bdfg, log, log_name + '_base_Dfg', output_path)


if __name__ == '__main__':
    print(len(sys.argv))
    print(sys.argv)
    if len(sys.argv) > 2:
        start(sys.argv[1], sys.argv[2])
    else:
        start(sys.argv[1])
