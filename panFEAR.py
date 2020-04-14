import os
import argparse
import scipy.stats as stats
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_path = r"C:\Users\aless\Documents\UniTn\SecondSemester\MicrobialGenomics\Project\Group_1" \
            r"\roary_out_i95_cd90_e_flag\gene_presence_absence.csv"

filter_hypothetical_proteins = True
correction = 'bh'

gene_p_a = pd.read_csv(file_path, header=0, usecols=(0, 2, 3))
pangenome_size = gene_p_a.shape[0]
functions_db = dict()
n_hypothetical = 0
group_genome = []


def check_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-gpa', '--gene_presence_absence', help='Roary output file "gene_presence_absence.csv" from which information is retrieved')
    parser.add_argument('-fhp', '--filter_hypothetical_proteins', action='store_true', help='Flag for excluding "hypothetical protein" from the analysed functions')
    parser.add_argument('-b', '--n_bins', type=int, help='Number of bins in the pangenome')
    parser.add_argument('-up', '--upper_percent', default=100, help='Upper percent among the pangenome of the group of which panFEAR analyses functional enrichment')
    parser.add_argument('-lp', '--lower_percent', default=0, help='Lower percent among the pangenome of the group of which panFEAR analyses functional enrichment')
    parser.add_argument('-c', '--correction', help='P-value correction method from multiple tests bias')
    parser.add_argument('-o', '--output_folder', help='ABSOLUTE! path where to store the produced outputs')  # TODO: try to generalize asap
    args = parser.parse_args()
    return args


def database_building(gene_p_a, n_bins, filter_hypothetical_proteins=True, upper_percent=100, lower_percent=0, functions_db=None):
    # Build the necessary structures
    min_strains_in_core = round(float(lower_percent) * n_bins / 100)
    max_strains_in_core = round(float(upper_percent) * n_bins / 100)
    pangenome_size = gene_p_a.shape[0]
    group_genome = []
    if not functions_db:
        functions_db = dict()
        n_hypothetical = 0
        for i, row in gene_p_a.iterrows():
            if filter_hypothetical_proteins:
                if row[1] != 'hypothetical protein':
                    if not row[1] in functions_db:
                        functions_db[row[1]] = [row[0]]
                    else:
                        functions_db[row[1]].append(row[0])
                    if min_strains_in_core <= int(row[2]) <= max_strains_in_core:
                        group_genome.append(row[0])
                else:
                    n_hypothetical += 1
            else:
                print('Chosen not to remove hypothetical proteins from FEA, possible unbalances during analysis')
                print(f'Pangenome size is: {pangenome_size}')
                if not row[1] in functions_db:
                    functions_db[row[1]] = [row[0]]
                else:
                    functions_db[row[1]].append(row[0])
                if min_strains_in_core <= int(row[2]) <= max_strains_in_core:
                    group_genome.append(row[0])
        if filter_hypothetical_proteins:
            print('Chosen to neglect hypothetical proteins from FEA')
            print(f'Pangenome size before hypothetical proteins filtering is: {pangenome_size}')
            pangenome_size -= n_hypothetical
            print(f'New pangenome size is: {pangenome_size}')
        group_genome_size = len(group_genome)
        print(f'Selected group genome size is: {group_genome_size}')
        out_of_group_genome_size = pangenome_size-group_genome_size
    else:
        for i, row in gene_p_a.iterrows():
            if filter_hypothetical_proteins:
                if row[1] != 'hypothetical protein':
                    if min_strains_in_core <= int(row[2]) <= max_strains_in_core:
                        group_genome.append(row[0])
            else:
                if min_strains_in_core <= int(row[2]) <= max_strains_in_core:
                    group_genome.append(row[0])
        group_genome_size = len(group_genome)
        print(f'Selected group genome size is: {group_genome_size}')
        out_of_group_genome_size = pangenome_size - group_genome_size
    return functions_db, group_genome, group_genome_size, out_of_group_genome_size


def fisher_test(functions_db, group_genome, group_genome_size, out_of_group_genome_size):
    # Fisher Exact 2-sided Test
    # Contingency table is:
    #                                 genes in group                   |                      genes not in group
    #     genes in funct          genes_in_funct_in_group              |         len(functions[funct]) - genes_in_funct_in_group
    # genes not in funct  group_genome_size - genes_in_funct_in_group  | accessory_genome_size - (len(functions[funct]) - genes_in_funct_in_group)

    tests_results = []
    for funct in functions_db:
        genes_in_funct_in_core = 0  # UpperLeft
        for gene in functions_db[funct]:
            if gene in group_genome:
                genes_in_funct_in_core += 1
        genes_not_in_funct_in_core = group_genome_size - genes_in_funct_in_core  # LowerLeft
        genes_in_funct_not_in_core = len(functions_db[funct]) - genes_in_funct_in_core  # UpperRight
        genes_not_in_funct_not_in_core = out_of_group_genome_size - genes_in_funct_not_in_core  # LowerRight
        odds_fract, p_value = stats.fisher_exact([[genes_in_funct_in_core, genes_in_funct_not_in_core],
                                                  [genes_not_in_funct_in_core, genes_not_in_funct_not_in_core]])
        tests_results.append((funct, p_value))
    tests_results = sorted(tests_results, key=lambda x: x[1])
    return tests_results


def pval_correct(tests_results, correction):
    # Statistical p-value corrections:
    global adj_p_values
    len_tests = len(tests_results)
    if correction.lower() in ['b', 'bonf', 'bonferroni']:
        print(f'Selected p-values correction is Bonferroni, {len_tests} tests done')
        adj_p_values = list(map(lambda x: x[1]*len_tests, tests_results))
    elif correction.lower() in ['bh', 'fdr']:
        print(f'Selected p-values correction is Benjamini-Hochberg, {len_tests} tests done')
        p_values = np.asarray(list(map(lambda x: x[1], tests_results)))
        Ks = np.arange(1, len_tests+1)/float(len_tests)
        adj_p_values_raw = p_values / Ks
        adj_p_values = np.minimum.accumulate(adj_p_values_raw[::-1])[::-1]
        adj_p_values[adj_p_values > 1] = 1
    tests_results = [(el[0], el[1], adj_p_values[i]) for i, el in enumerate(tests_results)]
    return tests_results


def store_results(tests_results, output_folder):
    # Write resulting dataframe to output, creates a csv
    res_df = pd.DataFrame(tests_results)
    res_df.to_csv(os.path.join(output_folder, 'panFEAR.csv'), header=False, index=False)


def plot_results(test_results):
    p_values = -np.log10(list(map(lambda x: x[-1], test_results))[:15])
    labels = list(map(lambda x: x[0], test_results))[:15]
    y_pos = np.arange(len(labels))
    plt.tight_layout()
    plt.barh(y_pos, p_values)
    plt.yticks(y_pos, labels)
    plt.show()


def main(args):
    gene_p_a_file = args.gene_presence_absence
    filter_hypothetical_proteins = args.filter_hypothetical_proteins
    n_bins = args.n_bins
    up = args.upper_percent
    lp = args.lower_percent
    correction = args.correction
    output_folder = args.output_folder
    gene_p_a = pd.read_csv(gene_p_a_file, header=0, usecols=(0, 2, 3))
    functions, group_genome, group_genome_size, out_of_group_genome_size = database_building(gene_p_a, n_bins, filter_hypothetical_proteins, up, lp)
    tests_results = fisher_test(functions, group_genome, group_genome_size, out_of_group_genome_size)
    if correction:
        tests_results = pval_correct(tests_results, correction)
    store_results(tests_results, output_folder)
    plot_results(tests_results)


if __name__ == '__main__':
    args = check_args()
    main(args)

