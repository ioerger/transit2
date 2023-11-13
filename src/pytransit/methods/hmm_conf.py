import sys, numpy
import scipy.stats

states = ["ES","GD","NE","GA"]
new_column_names = [
    "consis",
    "probES",
    "probGD",
    "probNE",
    "probGA",
    "conf",
    "flag",
]
def calc_probability(sat, non_zero_mean, mean_sat, stdev_sat, mean_non_zero_mean, stdev_non_zero_mean,):
    a = scipy.stats.norm.pdf(sat, loc=mean_sat, scale=stdev_sat)
    b = scipy.stats.norm.pdf(non_zero_mean, loc=mean_non_zero_mean, scale=stdev_non_zero_mean)
    return a * b


def get_rows(genes_file_path):
    rows = []
    for line in open(genes_file_path):
        line = line.strip()
        if line[0] == "#":
            continue
        row_elements = line.split("\t")
        rows.append(row_elements)
    
    return rows

def process_rows(rows):
    consistency_values = []
    orf_ids            = []
    mean_insertions_per_gene    = {}
    nz_means_per_gene           = {}
    calls_per_gene              = {}
    for orf, gene_name, description, total_sites, es_count, gd_count, ne_count, ga_count, mean_insertions, mean_reads, state_call in rows:
        total_sites     = int(total_sites)
        votes           = [ int(es_count), int(gd_count), int(ne_count), int(ga_count) ]
        mean_insertions_per_gene[orf] = float(mean_insertions)
        nz_means_per_gene[orf]        = float(mean_reads)
        calls_per_gene[orf]           = state_call
        if total_sites == 0:
            continue
        orf_ids.append(orf)
        consistency = max(votes) / float(total_sites)
        consistency_values.append(consistency)
    
    return consistency_values, calls_per_gene, orf_ids, nz_means_per_gene, mean_insertions_per_gene

def first_pass(all, calls, data, headers, nz_means, sats):
    print("# avg gene-level consistency of HMM states: %s" % (round(numpy.mean(all), 4)))
    id_s = [x[0] for x in data]

    mean_sats, mean_non_zero_means = {}, {}
    std_sats, stdev_non_zero_means = {}, {}

    print("# state posterior probability distributions:")
    for st in ["ES", "GD", "NE", "GA"]:
        sub                 = list(filter(lambda id: calls[id] == st, id_s))
        mean_sat            = numpy.mean([sats[id] for id in sub])
        mean_non_zero_mean  = numpy.mean([nz_means[id] for id in sub])
        stdev_sat           = numpy.std([sats[id] for id in sub])
        stdev_non_zero_mean = numpy.std([nz_means[id] for id in sub])
        print(
            "#  %s: genes=%s, meanSat=%s, stdSat=%s, meanNZmean=%s, stdNZmean=%s"
            % (
                st,
                len(sub),
                round(mean_sat, 3),
                round(stdev_sat, 3),
                round(mean_non_zero_mean, 1),
                round(stdev_non_zero_mean, 1),
            )
        )
        mean_sats[st] = mean_sat
        mean_non_zero_means[st] = mean_non_zero_mean
        std_sats[st] = mean_sat
        stdev_non_zero_means[st] = mean_non_zero_mean
    
    return mean_non_zero_means, mean_sats, std_sats, stdev_non_zero_means
    
# first pass...
def compute(genes_file_path):
    all, calls, data, headers, nz_means, sats = read_hmm_genes_file(genes_file_path)
    
    mean_non_zero_means, mean_sats, std_sats, stdev_non_zero_means = first_pass(all, calls, data, headers, nz_means, sats)
    
    
    new_rows = []
    # second pass...
    
    new_rows.append(
        headers[-1].split("\t")+new_column_names  # assume last comment line is column headers
    )
    for line in open(genes_file_path):
        line = line.strip()
        if line[0] == "#":
            continue
        row_elements = line.split("\t")
        id                 = row_elements[0]
        number_of_ta_sites = int(row_elements[3])
        votes              = [ int(x) for x in row_elements[4:8] ]
        sat                = float(row_elements[-3])
        non_zero_mean      = float(row_elements[-2])
        call               = row_elements[-1]
        
        if number_of_ta_sites == 0:
            rows.append(row_elements)
            continue
        
        probabilities = [
            calc_probability(
                sat=sat,
                non_zero_mean=non_zero_mean,
                mean_sat=mean_sats[each_state],
                stdev_sat=std_sats[each_state],
                mean_non_zero_mean=mean_non_zero_means[each_state],
                stdev_non_zero_mean=stdev_non_zero_means[each_state],
            )
                for each_state in states
        ]
        to_t_prob = float(sum(probabilities))
        relative_probabilites = [each / to_t_prob for each in probabilities]
        confidence = relative_probabilites[states.index(call)]
        flag = ""
        if confidence < 0.5:
            flag = "low-confidence"
        if max(relative_probabilites) < 0.7:
            if (call == "ES" or call == "GD") and (
                relative_probabilites[0] > 0.25 and relative_probabilites[1] > 0.25
            ):
                flag = "ambiguous"
            if (call == "GD" or call == "NE") and (
                relative_probabilites[1] > 0.25 and relative_probabilites[2] > 0.25
            ):
                flag = "ambiguous"
            if (call == "NE" or call == "GA") and (
                relative_probabilites[2] > 0.25 and relative_probabilites[3] > 0.25
            ):
                flag = "ambiguous"
        
        consistency = max(votes) / float(number_of_ta_sites)
        rounded_relative_probabilites = [ round(each, 6) for each in relative_probabilites ]
        
        new_row = list(row_elements)
        new_row += [ round(consistency, 3) ]
        new_row += [ *rounded_relative_probabilites ]
        new_row += [ round(confidence, 4) ]
        new_row += [ flag ]
        
        new_rows.append(
            new_row
        )


compute(genes_file_path=sys.argv[1])