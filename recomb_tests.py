"""
TODO:
- Write the description of the script
- Run it with two (three?) clades from Avachinskiy volcano, and each clade separately
- Run it with two clades from Khibine Mts, and each clade separately
"""

import matplotlib.pyplot as plt
import statistics
import numpy as np
import scipy.stats as st
from datetime import datetime
from recombfunctions import combfun

"""
Variables that are intended to be specified by user
"""

input_table = "./input/Khibine_3-gene_combinations-Hc.txt" # TSV table with the following fields: Region, SSU, EF1A_1, EF1A_2, COX1, Clade
gener_num = 1000  # number of F1 generations to be produced

# Variables for testing the influence of sample size on the number of genetic combinations
gen_size_start = 100  # starting generation size
gen_size_end = 5000  # final generation size
gen_size_step = 100  # step of the increment of the generation size

"""
define a log file with current timestamp
"""

now = datetime.now()  # current date and time
timestamp = datetime.timestamp(now)
logfile = f"./output/log_{timestamp}.txt"
print(f"Log and results will be written to {logfile}")

"""
Read the input
"""

individuals_list = []  # list of the original individuals (observed combinations)

# Dictionaries of genetic variants and phylogroups they belong to (Variant: Clade)
ssu_dict = {}
ef1a_dict = {}
coi_dict = {}

log_record = f"Reading the input file {input_table}"
combfun.write_log(logfile, log_record)

with open(input_table, "r") as f:
    lines = f.readlines()

# Put the observed genetic combinations into individuals_list, and genetic variants into dictionaries
for line in lines[1:]:
    if line != "":
        ssu, ef1, ef2, coi, clade = line.strip().split("\t")[1:]
        individuals_list.append([ssu, [ef1, ef2], coi])
        if ssu not in ssu_dict:
            ssu_dict[ssu] = clade
        if ef1 not in ef1a_dict:
            ef1a_dict[ef1] = clade
        if ef2 not in ef1a_dict:
            ef1a_dict[ef2] = clade
        if coi not in coi_dict:
            coi_dict[coi] = clade

log_record = f"The input file contains genetic combinations for {len(individuals_list)} individuals.\n"
combfun.write_log(logfile, log_record)

"""
Generate the specified number of F1 generations of the size == observed sample size
"""

gen_size = len(individuals_list)
gen_list = []

log_record = f"Creating {gener_num} F1 generations, each containing {len(individuals_list)} randomly " \
             f"generated individuals..."
combfun.write_log(logfile, log_record)

for i in range(gener_num):
    new_gen = [combfun.combine(combfun.rand_ind(individuals_list), combfun.rand_ind(individuals_list))
               for j in range(gen_size)]
    gen_list.append(new_gen)

log_record = "Done.\n"
combfun.write_log(logfile, log_record)

"""
Count the average number of unique combinations per generation
"""

comb_counts = []

for gen in gen_list:
    comb_set = set()
    for comb in gen:
        comb[1] = sorted(comb[1])
        comb_set.add(str(comb))
    comb_counts.append(len(comb_set))

"""
Test if the number of unique genetic combinations in a sample of defined size is normally
distributed in simulated generations (D’Agostino’s K^2 test)
"""

normality = st.normaltest(comb_counts)
log_record = f"Testing if the number of unique genetic combinations in a sample of defined size is normally " \
             f"distributed in simulated generations.\n" \
             f"D’Agostino’s K^2 normality test: statistic={normality[0]:.4f}, p-value={normality[1]:.4f}."
combfun.write_log(logfile, log_record)

"""
Compare the average number of genetic combinations in simulated F1 with the number of observed combinations
"""

ind_set = set()

for ind in individuals_list:
    ind_set.add(str(ind))

log_record = f"The average number of combinations per generation of the size {len(individuals_list)} is " \
             f"{sum(comb_counts)/gener_num} (95% confidence interval " \
             f"{st.t.interval(0.95, len(comb_counts)-1, loc=np.mean(comb_counts), scale=st.sem(comb_counts))}), " \
             f"1 stdev: from {np.mean(comb_counts) - statistics.stdev(comb_counts):.4f} to " \
             f"{np.mean(comb_counts) + statistics.stdev(comb_counts):.4f}, " \
             f"2 stdev: from {np.mean(comb_counts) - 2 * statistics.stdev(comb_counts):.4f} to " \
             f"{np.mean(comb_counts) + 2 * statistics.stdev(comb_counts):.4f}, " \
             f"3 stdev: from {np.mean(comb_counts) - 3 * statistics.stdev(comb_counts):.4f} to " \
             f"{np.mean(comb_counts) + 3 * statistics.stdev(comb_counts):.4f}, " \
             f"the observed number of combinations is {len(ind_set)}\n"
combfun.write_log(logfile, log_record)

"""
Write the first five F1 generations for recombination graph
"""
log_record = "Writing the first five F1 generations into TSV tables."
combfun.write_log(logfile, log_record)

for i, gen in enumerate(gen_list[:5]):
    recomb_out = ["SSU\tEF1A_1\tEF1A_2\tCOI"]
    for ind in gen:
        recomb_out.append(f"{ind[0]}\t{ind[1][0]}\t{ind[1][1]}\t{ind[2]}")

    log_record = f"Writing generation {i + 1} into file recombination_generation_{i + 1}.txt"
    combfun.write_log(logfile, log_record)

    with open(f"./output/recombination_generation_{i+1}_{timestamp}.txt", "w") as out:
        out.write("\n".join(recomb_out))

log_record = "Done.\n"
combfun.write_log(logfile, log_record)

"""
Count the number of “mixed” multi-locus genotypes (simulated individuals that have genetic variants coming
from different phylogroups)
"""
mixed_over_generations = []

for gen in gen_list:
    mixed_in_generation = 0
    for ind in gen:
        clade_set = set()
        clade_set.add(ssu_dict[ind[0]])
        clade_set.add(ef1a_dict[ind[1][0]])
        clade_set.add(ef1a_dict[ind[1][1]])
        clade_set.add(coi_dict[ind[2]])
        if len(clade_set) > 1:
            mixed_in_generation += 1
    mixed_over_generations.append(mixed_in_generation / len(individuals_list))

log_record = f"Average proportion of 'mixed' multi-locus genotypes among the simulated generations " \
             f"(simulated individuals that have genetic variants coming from different phylogroups): " \
             f"{sum(mixed_over_generations) / len(mixed_over_generations):.4f}, " \
             f"min {min(mixed_over_generations):.4f}, max {max(mixed_over_generations):.4f}"
combfun.write_log(logfile, log_record)

"""
Plot the dependence between the generation size and the number of unique combinations produced, if every 
generation is sub-sampled down to the observed sample size  
"""
#
# log_record = f"Now testing the influence of the generation size.\n" \
#              f"Generation sizes between {gen_size_start} and {gen_size_end} individuals will be tested " \
#              f"with an increment of {gen_size_step}. For every generation size, {gener_num} replications will be " \
#              f"created, and the average number of unique genetic combinations will be calculated. In addition, " \
#              f"every generation will be sub-sampled down to {len(individuals_list)} (the original sample size), " \
#              f"and the average number of unique genetic combinations will be calculated.\n"
# combfun.write_log(logfile, log_record)
#
# comb_per_gen = []
# comb_per_gen_all = []
# gen_size = []
#
# for k in range(gen_size_start, gen_size_end, gen_size_step):
#     gen_list = []
#     for i in range(gener_num):
#         new_gen = [combfun.combine(combfun.rand_ind(individuals_list), combfun.rand_ind(individuals_list)) for j in
#                    range(k)]
#         gen_list.append(new_gen)
#
#     subcomb_counts = []
#     comb_counts = []
#
#     for gen in gen_list:
#         comb_set = set()
#         for comb in gen:
#             comb_set.add(str(comb))
#         comb_counts.append(len(comb_set))
#
#         sub_comb_set = set()
#         for comb in gen[:len(individuals_list)]:  # sub-sampling to the observed sample size
#             sub_comb_set.add(str(comb))
#         subcomb_counts.append(len(sub_comb_set))
#
#     gen_size.append(k)
#     avg_comb_count = sum(comb_counts) / gener_num
#     subsampled_avg_comb_count = sum(subcomb_counts) / gener_num
#
#     log_record = f"Size {k}: {avg_comb_count} unique combinations; " \
#                  f"Sub-sampled: {subsampled_avg_comb_count} unique combinations"
#     combfun.write_log(logfile, log_record)
#
#     comb_per_gen.append(subsampled_avg_comb_count)
#     comb_per_gen_all.append(avg_comb_count)
#
# log_record = "Done.\n\nShowing the plot of the dependence between generation size and the average number of unique " \
#              "combinations of three independently inherited genes.\nYou can modify it and save manually.\n"
# combfun.write_log(logfile, log_record)
#
# plt.plot(comb_per_gen_all, gen_size)
# plt.ylabel('generation size')
# plt.xlabel('average number of unique combinations')
# plt.show()
#
# log_record = f"Showing the same plot for the generations sub-sampled down to the " \
#              f"observed sample size ({len(individuals_list)}).\nYou can modify it and save manually.\n"
# combfun.write_log(logfile, log_record)
#
# plt.plot(comb_per_gen, gen_size)
# plt.ylabel('generation size (before sub-sampling)')
# plt.xlabel('average number of unique combinations')
# plt.show()

log_record = "The job finished successfully."
combfun.write_log(logfile, log_record)
