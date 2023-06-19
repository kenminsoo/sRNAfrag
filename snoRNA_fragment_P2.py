from gtf_scripts.base.gtf_groundtruth import *
from gtf_scripts.base.conversion_tools import *
from gtf_scripts.base.basics import *
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import random
import yaml

## -- Config Variables -- ##

with open("snoRNA_frag_config.yaml", "r") as file:
    config_vars = yaml.safe_load(file)

# Working directory - All intermediate files will be deleted
working_dir = config_vars["dir_locations"]["working_dir"]
out_dir = config_vars["dir_locations"]["out_dir"]

# Reference genome location
build_bool = config_vars["module_options"]["P2"]["alignment2"]["build_masked"]["bool"]
reference_genome = config_vars["module_options"]["P1"]["build_index"]["reference_location"]
indexed = config_vars["module_options"]["P2"]["alignment2"]["index_location"]

# Attribute Choice => What to look for, biotype? transcript IDS?
attribute_choice = config_vars["module_options"]["P2"]["look_for"]

# Full annotation
full_annotation = config_vars["module_options"]["P2"]["annotation_file"]

# Make a lot of figures
make_lots = config_vars["module_options"]["P2"]["plot_every_source"]

# Prefix
prefix = config_vars["module_options"]["P1"]["prefix"]

## -- Config End -- ##

# Read in data and plot

counts_dataset = pd.read_csv("filtered_corrected_counts.csv")

plt.hist(counts_dataset["num"],bins = 50)
plt.title("# Filter Passing Fragment Loci Sources")
plt.xlabel("Number of Sources")
plt.savefig("Filter_Passing_Hist.jpeg", dpi = 500)
plt.clf()

# Create data related to the network creation

sources = list(counts_dataset["sources"])
counts = list(counts_dataset["adjusted_sums"])
lengths = list(counts_dataset["length"])
ids = list(counts_dataset["ID"])

network_building_dict = my_dictionary()

k = 0
for i in sources:
    temp_dict = my_dictionary()
    
    source = sources[k]

    split_source = source.split(sep  = ",")

    split_source2 = [i.split(sep = "_") for i in split_source]

    positions = [i[0].split(sep = ";") for i in split_source2]

    count = counts[k]
    length = lengths[k]
    id = ids[k]

    network_building_dict.add(id, [])

    j = 0 
    for source in split_source2:
        temp_dict.add(source[1], positions[j])
        j += 1

    network_building_dict[id].append(temp_dict)

    network_building_dict[id].append(count)

    network_building_dict[id].append(length)

    k += 1

# Create the graph and generate some figures if desired
connection_graph = nx.Graph()

for key in network_building_dict:
    entry_list = network_building_dict[key][0]

    count = network_building_dict[key][1]

    length = network_building_dict[key][2]

    connection_graph.add_nodes_from([key])

    connection_graph.add_nodes_from(list(entry_list.keys()))

    locations = network_building_dict[key]

    for i in entry_list:
        locations = entry_list[i]
        start = int(locations[0])
        end = int(locations[1])

        connection_graph.add_edge(key, i)

        connection_graph.edges[key, i]["counts"] = count
        connection_graph.edges[key, i]["length"] = length
        connection_graph.edges[key, i]["start"] = start
        connection_graph.edges[key, i]["end"] = end

# extract subgraphs

S = [connection_graph.subgraph(c).copy() for c in nx.connected_components(connection_graph)]

# note: Import the fragmentation prefix

prefix = prefix

to_test_graph = my_dictionary()

for i in S:
    for m,j,k, in i.edges(data = True):
        thing = [m,j,k]

        source_name = j
        ident = m 

        if prefix in m:
            switch = True
        else:
            source_name = m
            ident = j

        if source_name not in to_test_graph:

            to_test_graph.add(source_name, my_dictionary())
            to_test_graph[source_name].add("start", [])
            to_test_graph[source_name].add("end", [])
            to_test_graph[source_name].add("counts", [])
            to_test_graph[source_name].add("orig_ident", [])

        to_test_graph[source_name]["start"].append(k["start"])
        to_test_graph[source_name]["end"].append(k["end"])
        to_test_graph[source_name]["counts"].append(k["counts"])
        to_test_graph[source_name]["orig_ident"].append(ident)

if make_lots == True:
    os.system("mkdir source_peaks")

for key in to_test_graph:

    test_dataframe = pd.DataFrame.from_dict(to_test_graph[key])
    test_dataframe_small = test_dataframe[["start", "end", "counts"]]


    agg_start = test_dataframe_small.groupby("start").agg(sum)
    agg_start = agg_start.reset_index()

    agg_end = test_dataframe_small.groupby("end").agg(sum)
    agg_end = agg_end.reset_index()
    if make_lots == True:
        key_mod = key.replace("/", '.')
        plt.scatter(agg_start["start"], agg_start["counts"], label = "Start")
        plt.scatter(agg_end["end"], agg_end["counts"], label = "End")
        plt.legend()
        plt.xlabel("Position (n.t.)")
        plt.ylabel("Adjusted Counts")
        plt.title(key)
        plt.savefig("source_peaks/" + key_mod + ".jpeg", dpi = 200)
        plt.clf()

if make_lots == True:
    os.system("mv source_peaks " + out_dir)

# Detect Peaks #
x = 1
dis = 0
agg = 0

# Phase one clustering
# Figure creation to show how clusters are made

potential_loci_tracker = my_dictionary()
potential_loci_tracker.add("source_id", [])
potential_loci_tracker.add("loci", [])
potential_loci_tracker.add("type", [])

alt_peak_tracker = []
alt_peak_tracker_num = []
alt_peak_names = []

for key in to_test_graph:

    test_dataframe = pd.DataFrame.from_dict(to_test_graph[key])

    test_dataframe = test_dataframe[["start", "end", "counts"]]

    agg_start = test_dataframe.groupby("start").agg(sum)
    agg_start = agg_start.reset_index()

    start_max = max(agg_start["start"]) + 1

    agg_start_filled = pd.DataFrame(pd.Series(range(0, start_max)))
    agg_start_filled = agg_start_filled.rename(columns = {0:"start"})
    agg_start_filled = agg_start_filled.join(agg_start.set_index("start"), on = "start")
    agg_start_filled = agg_start_filled.fillna(0)

    agg_end = test_dataframe.groupby("end").agg(sum)
    agg_end = agg_end.reset_index()

    end_max = max(agg_end["end"]) + 1

    agg_end_filled = pd.DataFrame(pd.Series(range(0, end_max)))
    agg_end_filled = agg_end_filled.rename(columns = {0:"end"})
    agg_end_filled = agg_end_filled.join(agg_end.set_index("end"), on = "end")
    agg_end_filled = agg_end_filled.fillna(0)

    agg_start_filled["before"] = agg_start_filled["counts"].shift(1)
    agg_start_filled = agg_start_filled.fillna(0)

    agg_start_filled["delta_start"] = agg_start_filled["counts"] - agg_start_filled["before"]

    agg_start_filled["start_delt_pos"] = agg_start_filled["delta_start"] > 0
    agg_start_filled["start_delt_next"] = agg_start_filled["start_delt_pos"].shift(-1)
    agg_start_filled = agg_start_filled.fillna(0)

    agg_start_filled["start_peak"] = (agg_start_filled["start_delt_pos"] == True) & (agg_start_filled["start_delt_next"] == False)
    ## END LOCI

    agg_end_filled["before"] = agg_end_filled["counts"].shift(1)
    agg_end_filled = agg_end_filled.fillna(0)

    agg_end_filled["delta_start"] = agg_end_filled["counts"] - agg_end_filled["before"]

    agg_end_filled["start_delt_pos"] = agg_end_filled["delta_start"] > 0
    agg_end_filled["start_delt_next"] = agg_end_filled["start_delt_pos"].shift(-1)
    agg_end_filled = agg_end_filled.fillna(0)

    agg_end_filled["end_peak"] = (agg_end_filled["start_delt_pos"] == True) & (agg_end_filled["start_delt_next"] == False)

    ## Create mini dfs to join
    agg_start_mini = agg_start_filled[["start", "start_peak"]]

    agg_end_mini = agg_end_filled[["end", "end_peak"]]

    agg_start_mini["type"] = "start"
    agg_start_mini = agg_start_mini.rename(columns={"start":"loci", "start_peak":"peak"})

    agg_end_mini["type"] = "end"
    agg_end_mini = agg_end_mini.rename(columns={"end":"loci", "end_peak":"peak"})

    peaks_dataframe = pd.concat([agg_start_mini, agg_end_mini])
    peaks_dataframe = peaks_dataframe.sort_values("loci")

    clusters = []

    naming_start_cond = True

    current_id = 1

    for row in peaks_dataframe.iterrows():

        row_data = list(row[1])

        type_loci = row_data[2]
        peak_status = row_data[1]

        if peak_status:
            potential_loci_tracker["source_id"].append(key)
            potential_loci_tracker["loci"].append(row_data[0])
            potential_loci_tracker["type"].append(type_loci)

        if type_loci == 'start':
            # first check for the start condition
            if naming_start_cond == False:
                if peak_status:
                    naming_start_cond = True
            
            clusters.append(current_id)

        else:
            if naming_start_cond:
                if peak_status:
                    clusters.append(current_id)
                    current_id += 1
                    naming_start_cond = False
                else:
                    clusters.append(current_id)
            else:
                clusters.append(current_id - 1)

    copy_of_test_graph = to_test_graph.copy()

    peaks_dataframe["cluster"] = clusters

    copy_of_test_graph[key].add("cluster", [])

    # test if end and beginning agree
    inde = 0
    for item in to_test_graph[key]["start"]:
        current_start = item
        current_end = copy_of_test_graph[key]["end"][inde]
        
        current_cluster_start = list(peaks_dataframe[(peaks_dataframe["loci"] == current_start) & (peaks_dataframe["type"] == "start")]["cluster"])
        if len(current_cluster_start) != 1:
            raise ValueError('This should not occur. Submit issue on Github.')
        
        current_cluster_start = current_cluster_start[0]

        current_cluster_end = list(peaks_dataframe[(peaks_dataframe["loci"] == current_end) & (peaks_dataframe["type"] == "end")]["cluster"])
        if len(current_cluster_end) != 1:
            raise ValueError('This should not occur. Submit issue on Github.')
        
        current_cluster_end = current_cluster_end[0]

        if current_cluster_end == current_cluster_start:
            copy_of_test_graph[key]["cluster"].append(current_cluster_start)
            agg += 1

        else:
            moddy_id = copy_of_test_graph[key]["orig_ident"][inde]
            alt_peak_names.append(moddy_id)
            alt_peak_tracker.append(key)
            alt_naming = str(current_cluster_start) + "." + str(current_cluster_end)
            alt_naming = float(alt_naming)
            alt_peak_tracker_num.append(alt_naming)
            copy_of_test_graph[key]["cluster"].append(alt_naming)
            dis += 1

        inde += 1

# Export peaks
potential_peaks = pd.DataFrame.from_dict(potential_loci_tracker)
potential_peaks.to_csv("potential_peaks.csv", index = False)

os.system("mv potential_peaks.csv " + out_dir)

alternate_peaks = {"source_id":alt_peak_tracker, "start_end":alt_peak_tracker_num, "original_id":alt_peak_names}
alt_peaks_df = pd.DataFrame.from_dict(alternate_peaks)
alt_peaks_df.to_csv("alternate_peaks.csv", index = False)
os.system("mv alternate_peaks.csv " + out_dir)

# Cluster Peaks

peak_clustering_graph = nx.Graph()
peak_clustering_graph.add_nodes_from(ids)

for i in peak_clustering_graph.nodes():
    peak_clustering_graph.nodes[i]["number_of_mentions"] = 0
    peak_clustering_graph.nodes[i]["disagreements"] = 0

for i in S:
    for m,j,k, in i.edges(data = True):
        source_name = j
        ident = m 

        if prefix in m:
            switch = True
        else:
            source_name = m
            ident = j

        # pull relevant data

        cluster_values = np.unique(copy_of_test_graph[source_name]["cluster"])

        for value in cluster_values:
            
            indexes = list(np.where(np.array(copy_of_test_graph[source_name]["cluster"]) == value)[0])

            cluster_ids = []

            if "." in str(value):
                for index in indexes:
                    index_associated_id = copy_of_test_graph[source_name]["orig_ident"][index]
                    peak_clustering_graph.nodes[index_associated_id]["number_of_mentions"] += 1
                    peak_clustering_graph.nodes[index_associated_id]["disagreements"] += 1
                    cluster_ids.append(index_associated_id)

                if len(cluster_ids) == 1:
                    continue
                else:
                    kk = 0 
                    for cluster_draw in cluster_ids:
                        if kk == 0:
                            kk += 1
                            continue
                        else:
                            peak_clustering_graph.add_edge(cluster_draw, cluster_ids[kk - 1])

                            if "weight" in peak_clustering_graph.edges[cluster_draw, cluster_ids[kk - 1]]:
                                peak_clustering_graph.edges[cluster_draw, cluster_ids[kk - 1]]["weight"] += 1
                            else:
                                peak_clustering_graph.edges[cluster_draw, cluster_ids[kk - 1]]["weight"] = 1

                            kk += 1

            else:
                for index in indexes:
                    index_associated_id = copy_of_test_graph[source_name]["orig_ident"][index]
                    peak_clustering_graph.nodes[index_associated_id]["number_of_mentions"] += 1
                    cluster_ids.append(index_associated_id)

                if len(cluster_ids) == 1:
                    continue
                else:
                    kk = 0
                    for cluster_draw in cluster_ids:
                        if kk == 0:
                            kk += 1
                            continue
                        else:
                            peak_clustering_graph.add_edge(cluster_draw, cluster_ids[kk - 1])

                            if "weight" in peak_clustering_graph.edges[cluster_draw, cluster_ids[kk - 1]]:
                                peak_clustering_graph.edges[cluster_draw, cluster_ids[kk - 1]]["weight"] += 1
                            else:
                                peak_clustering_graph.edges[cluster_draw, cluster_ids[kk - 1]]["weight"] = 1

                            kk += 1

# Create merged IDS
S_fragments = [peak_clustering_graph.subgraph(c).copy() for c in nx.connected_components(peak_clustering_graph)]

# Calculate disagreement to full ratio
id_counter = 1
new_id = prefix + "_merged_"
number_dis = 0

merged_references = my_dictionary()
merged_references.add("original_id", [])
merged_references.add("new_id", [])
merged_references.add("in_two_peaks", [])

ratio = []
thetotal = 0
for i in S_fragments:
    # obtain node metrics
    for k in i.nodes(data=True):
        if k[1]["disagreements"] != 0:
            ratio.append(k[1]["disagreements"] / k[1]["number_of_mentions"])

            number_dis += 1

            if k[1]["disagreements"] / k[1]["number_of_mentions"] == 1:
                merged_references["in_two_peaks"].append(True)

            else:
                merged_references["in_two_peaks"].append(round(k[1]["disagreements"] / k[1]["number_of_mentions"], 3))
        else:
            merged_references["in_two_peaks"].append(False)

    for k in i.nodes():
        merged_references["original_id"].append(k)
        merged_references["new_id"].append(new_id + str(id_counter))
        thetotal += 1

    id_counter += 1

plt.hist(ratio)
plt.title("Ratio of Fragments Spanning >1 Cluster\n" + str(100 * round(number_dis / (thetotal - 1), 2)) + "% of Fragments")
plt.xlabel("Peaks != / Number of Mentions")
plt.ylabel("Frequency")
plt.savefig("cluster_start_end_disagreement.jpeg", dpi = 500)

ref_table = pd.DataFrame.from_dict(merged_references)
ref_table.to_csv("ref_table.csv", index = False)

os.system("mv ref_table.csv " + out_dir)
os.system("mv cluster_start_end_disagreement.jpeg " + out_dir)

# Merge Tables
joined_data = counts_dataset.join(ref_table.set_index("original_id"), on = "ID")
columns_joined = list(joined_data.columns)
end_index = columns_joined.index("sources")
new_id_index = columns_joined.index("new_id")

counts_table = joined_data.iloc[:,[new_id_index] + list(range(1, end_index))]
sum_table = counts_table.groupby("new_id").agg(sum)
sum_table = sum_table.reset_index()

sum_table.to_csv("merged_counts.csv", index = False)

os.system("mv merged_counts.csv " + out_dir)

## == ## Detection of outside mapping ## == ##
num = list(counts_dataset.columns).index("ID")

tsv_to_fasta(out_dir + "/filtered_corrected_counts.csv", "filtered_sequences.fa",num, 0, delim = ",")

bowtie_align_pipeline(indexed, working_dir, out_dir + "/filtered_sequences.fa")

# Take out unaligned reads
os.system('cd ' + working_dir + "; \
            grep '	0	' lookup_filtered.sam | cut -f 1,2 |  awk -v OFS=',' '{print $1, $2}' | sort | uniq > all_reads.csv;\
            grep ',4' all_reads.csv > nomatches.csv;\
            grep ',0' all_reads.csv | cut -d , -f 1 > matches.txt;\
            grep ',16' all_reads.csv | cut -d , -f 1 >> matches.txt;\
            cat matches.txt | sort | uniq > matches_uniq.txt;\
            mkdir aligned_sams;\
            cat matches_uniq.txt | parallel 'grep @ lookup_filtered.sam > aligned_sams/{}.sam';\
            cat matches_uniq.txt | parallel 'grep {} lookup_filtered.sam >> aligned_sams/{}.sam';\
            featureCounts -a " + full_annotation + " -F 'GTF' -g " + attribute_choice + " -o matches.tsv aligned_sams/*.sam -O -M")

# Find out where the aligned reads mapped to
matches_df = pd.read_csv(working_dir + "/matches.tsv", sep = "\t", skiprows = 1)

matches_df = matches_df.set_index("Geneid")

summary_df = pd.read_csv(working_dir + "/matches.tsv.summary", sep = "\t", nrows = 1)

ids = list(summary_df.columns)

ids_clean = [n.replace("aligned_sams/", "") for n in ids]
ids_clean2 = [n.replace(".sam", "") for n in ids_clean]
ids_clean2.pop(0)

annotation_df_to_join = {"id":ids_clean2}

for row in summary_df.iterrows():
    assigned_list = list(row[1])
    
assigned_list.pop(0)

status = []
anno = []

i = 0

geneids = list(matches_df.index)


for item in assigned_list:
    if item == 0:
        status.append(True)
        
        anno.append("")
        
    else:
        status.append(True)
        
        id_name = ids_clean2[i]
        
        columnwithdata= matches_df["aligned_sams/"+ id_name + ".sam"]
        
        index_list = columnwithdata.to_numpy().nonzero()[0].tolist()
        
        storage_gene = []
        
        for indx in index_list:
            geneid = geneids[indx]
            
            storage_gene.append(geneid)
            
        entry_to_add = ";".join(storage_gene)
        
        anno.append(entry_to_add)
        
    i += 1

annotation_df_to_join.update({"annotation":anno})

converted_df_an = pd.DataFrame.from_dict(annotation_df_to_join)

annotated_ref_table = ref_table.join(converted_df_an.set_index("id"), on = "original_id", how = "left")

annotated_ref_table.to_csv("annotated_ref_table.csv", index = False)

# Make All Figures

if make_lots == True:
    annotated_ref_table
    potential_peaks