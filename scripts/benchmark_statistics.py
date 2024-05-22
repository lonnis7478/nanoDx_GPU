import os
import re
import pandas as pd
import math
import markdown as md

# training_set = snakemake.input["trainingset"]
training_set = "Capper_et_al"

rules = set()
samples = set()

rule_data = {}
sample_data = {}

for file in os.listdir("benchmarks"):
    if "fast5" not in file:
        file = file.replace("." + training_set, "")
        rule = re.findall("^\w*\.(\w+\.*\w*)\.benchmark.txt$", file)
        if len(rule) > 0:
            rules.add(rule[0])
        sample = re.findall("^(\w*)\.\w+\.*\w*\.benchmark.txt$", file)
        if len(sample) > 0:
            samples.add(sample[0])

for rule in rules:

    time_sum = 0.0
    file_counter = 0
    rule_data[rule] = {}
    for file in os.listdir("benchmarks"):
        if rule in file:
            file_counter += 1
            stats = pd.read_csv("benchmarks/" + file, sep="\t")
            time_sum += stats["s"][0]
            sample = re.findall("^(\w*)\.\w+\.*\w*\.benchmark.txt$", file)
            if len(sample) > 0:
                rule_data[rule][sample[0]] = stats["s"][0]

    mean = time_sum / file_counter
    rule_data[rule]["total"] = time_sum
    rule_data[rule]["mean"] = mean


cpg_df = pd.read_csv("benchmarks/CpGs_benchmark.txt", sep="\t", names=["sample", "cpgs"])

for sample in samples:

    time_sum = 0.0
    old_time = 0.0
    file_counter = 0

    for file in os.listdir("benchmarks"):
        if sample+"." in file and "fast5" not in file:
            file_counter += 1
            df = pd.read_csv("benchmarks/" + file, sep="\t")


            if "cuml_tsne" not in file and "get_tSNE_data" not in file and "plot_tSNE_CUDA" not in file:

                old_time += df["s"][0]

            if "plot_tSNE" not in file:
                time_sum += df["s"][0]
    data = {}
    cpg_sample = cpg_df.loc[cpg_df['sample'] == sample]

    if not cpg_sample["cpgs"].empty:
        data["cpgs"] = cpg_sample["cpgs"].iloc[0]
    data["total"] = time_sum
    data["total_old"] = old_time

    sample_data[sample] = data


print(rule_data)
print(sample_data)

html = md.markdown("# Benchmarks")
table = "| Rule | Mean time (s) |\n"
table += " |--------|--------|\n"
for rule in rule_data:
    table += "| " + rule + " |" + str(round(rule_data[rule]["mean"])) + " |\n"

html += md.markdown("## Rule statistics")
html += md.markdown(table, extensions=["tables"])

html += "<br/><br/>"

header = "| Sample | CpGs | Time (s) | Time New (min) | Time - original (min) |"
has_completedHeader = False
header_divider = " |--------|--------|--------|--------|--------|"
table = ""

for sample in sample_data:
    min = math.floor(sample_data[sample]["total"] / 60)
    seconds = round((sample_data[sample]["total"] / 60 - min) * 60)

    min_old = math.floor(sample_data[sample]["total_old"] / 60)
    seconds_old = round((sample_data[sample]["total_old"] / 60 - min) * 60)
    table += "|" + sample + "| " + str(sample_data[sample]["cpgs"]) + " |" + str(round(sample_data[sample]["total"])) + "|" + str(min)+"min " +str(seconds)+"s |" + str(min_old)+"min " +str(seconds_old)+"s |"
    for rule in rule_data:
        if not has_completedHeader:
            header += rule+"|"
            header_divider += "--------|"

        if sample in rule_data[rule]:
            table +=str(rule_data[rule][sample])+"|"
        else:
            table += "---| "
    has_completedHeader = True
    table += "\n"

html += md.markdown("## Sample statistics")

table = header +"\n"+ header_divider + "\n"+table
html += md.markdown(table, extensions=["tables"])

header = "| Sample | TSNE (s)| CUDA TSNE (s) | Speedup (X) |\n"
header_divider = " |--------|--------|--------|--------|\n"
table = ""

for sample in sample_data:
    if sample in rule_data["plot_tSNE"] and sample in rule_data["cuml_tsne"]:
        table += "| " + sample + "| "+str(rule_data["plot_tSNE"][sample])+" | "+str(rule_data["cuml_tsne"][sample])+" | "+str((rule_data["plot_tSNE"][sample] / rule_data["cuml_tsne"][sample]))+" |\n"

html += "<br/><br/>"
html += md.markdown("## TSNE vs. CUDA TSNE")

table = header + header_divider + table
html += md.markdown(table, extensions=["tables"])

html += "<style>table tr td{border:1px solid gray;} table th {border-right:1px solid black;} table{table-collapse:collapse;border:1px solid black;}</style>"
file = open("benchmark_stats.html", "w")
file.write(html)
file.close()
