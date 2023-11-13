import sys, numpy
import scipy.stats

def get_rows(genes_file_path):
    rows = []
    for line in open(genes_file_path):
        line = line.strip()
        if line[0] == "#":
            continue
        row_elements = line.split("\t")
        rows.append(
            (
                orf,
                gene_name,
                description,
                int(total_sites),
                int(es_count),
                int(gd_count),
                int(ne_count),
                int(ga_count),
                float(mean_insertions),
                float(mean_reads),
                state_call
            )
        )
    
    return rows

class HmmConfidenceHelper:
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
    
    @staticmethod
    def extract_data_from_rows(rows):
        consistency_values = []
        orf_ids            = []
        mean_insertions_per_gene    = {}
        nz_means_per_gene           = {}
        calls_per_gene              = {}
        for orf, gene_name, description, total_sites, es_count, gd_count, ne_count, ga_count, mean_insertions, mean_reads, state_call in rows:
            mean_insertions_per_gene[orf] = mean_insertions
            nz_means_per_gene[orf]        = mean_reads
            calls_per_gene[orf]           = state_call
            if total_sites == 0:
                continue
            orf_ids.append(orf)
            consistency = max([ es_count, gd_count, ne_count, ga_count ]) / float(total_sites)
            consistency_values.append(consistency)
        
        return consistency_values, calls_per_gene, orf_ids, nz_means_per_gene, mean_insertions_per_gene
    
    @staticmethod
    def first_pass(rows, debug=False):
        consistency_values, calls_per_gene, data, nz_means_per_gene, mean_insertions_per_gene = HmmConfidenceHelper.extract_data_from_rows(rows)
        
        debug and print("# avg gene-level consistency of HMM states: %s" % (round(numpy.mean(consistency_values), 4)))
        mean_sats, mean_non_zero_means = {}, {}
        stdev_sats, stdev_non_zero_means = {}, {}
        
        debug and print("# state posterior probability distributions:")
        for each_state in ["ES", "GD", "NE", "GA"]:
            relevent_orf_ids = [ each_orf_id for each_orf_id in orf_ids if calls_per_gene[orf_id] == each_state ]
            mean_sat            = numpy.mean([mean_insertions_per_gene[orf_id] for orf_id in relevent_orf_ids])
            mean_non_zero_mean  = numpy.mean([nz_means_per_gene[orf_id]        for orf_id in relevent_orf_ids])
            stdev_sat           = numpy.std([mean_insertions_per_gene[orf_id]  for orf_id in relevent_orf_ids])
            stdev_non_zero_mean = numpy.std([nz_means_per_gene[orf_id]         for orf_id in relevent_orf_ids])
            debug and print(
                "#  %s: genes=%s, meanSat=%s, stdSat=%s, meanNZmean=%s, stdNZmean=%s"
                % (
                    each_state,
                    len(relevent_orf_ids),
                    round(mean_sat, 3),
                    round(stdev_sat, 3),
                    round(mean_non_zero_mean, 1),
                    round(stdev_non_zero_mean, 1),
                )
            )
            mean_sats[each_state]            = mean_sat
            mean_non_zero_means[each_state]  = mean_non_zero_mean
            stdev_sats[each_state]           = stdev_sat # NOTE: original was being assigned to mean_sat (pretty sure it was a bug)
            stdev_non_zero_means[each_state] = stdev_non_zero_mean
        
        return mean_non_zero_means, mean_sats, stdev_sats, stdev_non_zero_means

    @staticmethod
    def second_pass(rows, mean_non_zero_means, mean_sats, stdev_sats, stdev_non_zero_means):
        row_extensions = []
        for orf, gene_name, description, total_sites, es_count, gd_count, ne_count, ga_count, mean_insertions, mean_reads, state_call in rows:
            if total_sites == 0:
                row_extensions.append(tuple())
                continue
            
            probabilities = [
                HmmConfidenceHelper.calc_probability(
                    sat=mean_insertions,
                    non_zero_mean=mean_reads,
                    mean_sat=mean_sats[each_state],
                    stdev_sat=stdev_sats[each_state],
                    mean_non_zero_mean=mean_non_zero_means[each_state],
                    stdev_non_zero_mean=stdev_non_zero_means[each_state],
                )
                    for each_state in HmmConfidenceHelper.states
            ]
            to_t_prob = float(sum(probabilities))
            relative_probabilites = [each / to_t_prob for each in probabilities]
            confidence = relative_probabilites[HmmConfidenceHelper.states.index(state_call)]
            flag = ""
            if confidence < 0.5:
                flag = "low-confidence"
            if max(relative_probabilites) < 0.7:
                if (state_call == "ES" or state_call == "GD") and (
                    relative_probabilites[0] > 0.25 and relative_probabilites[1] > 0.25
                ):
                    flag = "ambiguous"
                if (state_call == "GD" or state_call == "NE") and (
                    relative_probabilites[1] > 0.25 and relative_probabilites[2] > 0.25
                ):
                    flag = "ambiguous"
                if (state_call == "NE" or state_call == "GA") and (
                    relative_probabilites[2] > 0.25 and relative_probabilites[3] > 0.25
                ):
                    flag = "ambiguous"
            
            consistency = max([ es_count, gd_count, ne_count, ga_count ]) / float(total_sites)
            rounded_relative_probabilites = [ round(each, 6) for each in relative_probabilites ]
            
            row_extensions.append((
                round(consistency, 3),
                *rounded_relative_probabilites,
                round(confidence, 4),
                flag,
            ))
    
    @staticmethod
    def calc_probability(sat, non_zero_mean, mean_sat, stdev_sat, mean_non_zero_mean, stdev_non_zero_mean,):
        a = scipy.stats.norm.pdf(sat, loc=mean_sat, scale=stdev_sat)
        b = scipy.stats.norm.pdf(non_zero_mean, loc=mean_non_zero_mean, scale=stdev_non_zero_mean)
        return a * b
            
    @staticmethod
    def compute_row_extensions(rows):
        mean_non_zero_means, mean_sats, stdev_sats, stdev_non_zero_means = HmmConfidenceHelper.first_pass(
            rows
        )
        row_extensions = HmmConfidenceHelper.second_pass(
            rows,
            mean_non_zero_means,
            mean_sats,
            stdev_sats,
            stdev_non_zero_means
        )
        return row_extensions

row_extensions = HmmConfidenceHelper.compute_row_extensions(get_rows(genes_file_path=sys.argv[1]))