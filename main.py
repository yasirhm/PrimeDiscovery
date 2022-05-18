
import optimum_discovery
import frequency_based_heuristic as freq_discovery
import selective_discovery
import selective_heuristic


if __name__ == '__main__':

    method = ''
    log = ''
    if len(sys.argv) > 2:
        method = sys.argv[1]
        log = sys.argv[2]
    else :
        print("The command is: python main.py <method name> <event log>")

    
    logs_artificial = ['BPI_2020_Travel_Permit_arificial.xes', 'Hospital-Billing-Arttificial.xes',
                   "Sepsis Artificial Events.xes", "Road_Traffic_Artificial_Events.xes",
                   "BPI_Challenge_2019_Artificial_evenets.xes", "financial_log_2012_Artificial_Events.xes",
                   "RequestForPaymentArtificial.xes"]

    generated_logs = ["l1.xes", "l-choice.xes", "le3.xes", "le12.xes", "lchoice.xes"]
    path2 = "../generated-log/"

    if method = 'optimum_discovery':
        optimum_discovery.start(log)
    elif method = 'frequency_discovery':
        freq_discovery.main(log)
    elif method = 'selective_discovery':
        selective_discovery.main(log)