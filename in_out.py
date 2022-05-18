import xlsxwriter
from pm4py.visualization.dfg import visualizer as dfg_visualization
from pm4py.objects.dfg.exporter import exporter as dfg_exporter
from pm4py.objects.dfg.importer import importer as dfg_importer

parameters = {dfg_visualization.Variants.FREQUENCY.value.Parameters.FORMAT: "jpg"}


def export_excel_result(alignment_result, keys, outname, output_path, itera):
    workbook = xlsxwriter.Workbook(output_path + '/final_heuristics_' + outname + '.xlsx')
    worksheet = workbook.add_worksheet()
    col_num = 0
    worksheet.write(0, col_num, '')
    worksheet.write_column(1, col_num,
                           ['Fitness ' + str(itera), 'Precision ' + itera, 'ActualFM ' + itera, 'PossibleFM ' + itera])
    col_num = 1
    for key, value in alignment_result.items():
        col_key = keys[key][0] + ',' + keys[key][1]
        worksheet.write(0, col_num, col_key)
        worksheet.write_column(1, col_num, value)
        col_num += 1

    workbook.close()


def export_excel_final_result(alignment_result, out_name, output_path):
    workbook = xlsxwriter.Workbook(output_path + '/final_heuristic_' + out_name + '.xlsx')
    worksheet = workbook.add_worksheet()
    col_num = 0
    worksheet.write(0, col_num, '')
    worksheet.write_column(1, col_num,
                           ['Fitness ', 'Precision ', 'ActualFM ', 'PossibleFM '])
    col_num = 1
    for key, value in alignment_result.items():
        # col_key = keys[key][0] + ',' + keys[key][1]
        col_key = str(key[0]).join(' ').join(key[1])
        worksheet.write(0, col_num, col_key)
        worksheet.write_column(1, col_num, value[1])
        col_num += 1

    workbook.close()


def export_dfg_file(dfg, sac, eac, output_path, name):
    export_parameters = {dfg_exporter.Variants.CLASSIC.value.Parameters.START_ACTIVITIES: sac,
                         dfg_exporter.Variants.CLASSIC.value.Parameters.END_ACTIVITIES: eac}
    dfg_exporter.apply(dfg, output_path + '/' + name, variant=dfg_exporter.Variants.CLASSIC,
                       parameters=export_parameters)


def export_visual_dfg(indfg, log, name, output_path):
    gviz = dfg_visualization.apply(indfg, log=log, variant=dfg_visualization.Variants.FREQUENCY, parameters=parameters)
    dfg_visualization.save(gviz, output_path + '/' + name)
