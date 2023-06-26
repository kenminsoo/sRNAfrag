from gtf_groundtruth import *
from conversion_tools import *
from basics import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import random
import yaml
import warnings

## -- Config Variables -- ##

with open("sRNA_frag_config.yaml", "r") as file:
    config_vars = yaml.safe_load(file)

# Working directory - All intermediate files will be deleted
working_dir = config_vars["dir_locations"]["working_dir"]
out_dir = config_vars["dir_locations"]["out_dir"]

# Reference genome location
indexed = config_vars["module_options"]["P1"]["built_index_location"]

# Attribute Choice => What to look for, biotype? transcript IDS?
attribute_choice = config_vars["module_options"]["P2"]["look_for"]

# Full annotation
full_annotation = config_vars["module_options"]["P2"]["annotation_file"]

processed_annotation_file = config_vars["module_options"]["P1"]["annotation_options"]["location"]

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
plt.savefig("P2_Filter_Passing_Hist.jpeg", dpi = 500)
plt.clf()

# Create data related to the network creation

sources = list(counts_dataset["sources"])
counts = list(counts_dataset["sums"])
lengths = list(counts_dataset["length"])
ids = list(counts_dataset["ID"])

network_building_dict = my_dictionary()

k = 0
for i in sources:
    # Each row has a unique ID
    # Extract source information that connects to each ID
    temp_dict = my_dictionary()
    
    source = sources[k]

    split_source = source.split(sep  = ">")

    # nested list is generated
    split_source2 = [i.split(sep = "_") for i in split_source]

    positions = [i[0].split(sep = ";") for i in split_source2]

    count = counts[k]
    length = lengths[k]
    id = ids[k]

    # dictionary has id as the key
    network_building_dict.add(id, [])

    j = 0 
    for source in split_source2:

        temp_dict.add(source[1], positions[j])
        j += 1

    # then each id has information regarding count length and sources/positions
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

    connection_graph.add_nodes_from([key], node_type = "fragment")

    connection_graph.add_nodes_from(list(entry_list.keys()), node_type = "source")

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

# extract subgraphs with central node source
# S = [connection_graph.subgraph(c).copy() for c in nx.connected_components(connection_graph)]

S = [connection_graph.subgraph(list(connection_graph.neighbors(c)) + [c]).copy() for c,d in connection_graph.nodes(data = True) if d["node_type"] == "source"]

# note: Import the fragmentation prefix

prefix = prefix

to_test_graph = my_dictionary()

# Oh my, i understand the need to comment now

for i in S:

    # For each star graph loop through the edges
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
        # Original identity refers to the fragmen ID 
        # Source1 => identity1
        # Source2 => identity1 
        # This allows us to track such relationships
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

potential_loci_tracker = my_dictionary()
potential_loci_tracker.add("source_id", [])
potential_loci_tracker.add("loci", [])
potential_loci_tracker.add("type", [])
potential_loci_tracker.add("counts", [])

alt_peak_tracker = []
alt_peak_tracker_num = []
alt_peak_names = []

i_1 = 0

# Loop through each potential source
for key in to_test_graph:

    # Each source has information regarding fragments, loci, and counts
    test_dataframe = pd.DataFrame.from_dict(to_test_graph[key])

    test_dataframe = test_dataframe[["start", "end", "counts"]]

    # Count by start 
    # This part is fairly understandable in the manuscript
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

    ## START LOCI
    agg_start_filled["not_zero_counts"] = agg_start_filled["counts"] > 0
    agg_start_filled["not_zero_before"] = agg_start_filled["not_zero_counts"].shift(1)
    agg_start_filled["not_zero_after"] = agg_start_filled["not_zero_counts"].shift(-1)

    agg_start_filled = agg_start_filled.fillna(0)

    agg_start_filled["zero_smush"] = (agg_start_filled["not_zero_before"] == True) & (agg_start_filled["not_zero_after"] == True) & (agg_start_filled["not_zero_counts"] == False)
    agg_start_filled["next_count"] = agg_start_filled["counts"].shift(-1)
    agg_start_filled = agg_start_filled.fillna(0)
    agg_start_filled["additive_factor"] = agg_start_filled["next_count"] * agg_start_filled["zero_smush"]

    agg_start_filled["counts"] = agg_start_filled["counts"] + agg_start_filled["additive_factor"]

    agg_start_filled["before"] = agg_start_filled["counts"].shift(1)
    agg_start_filled = agg_start_filled.fillna(0)

    agg_start_filled["delta_start"] = (agg_start_filled["counts"] - agg_start_filled["before"])
    agg_start_filled["delta_start_next"] = agg_start_filled["delta_start"].shift(-1)
    agg_start_filled = agg_start_filled.fillna(0)

    agg_start_filled["delta_start_modded"] = agg_start_filled["delta_start"] + ((1/2) * agg_start_filled["delta_start_next"])

    agg_start_filled["start_delt_pos"] = agg_start_filled["delta_start_modded"] > 0
    agg_start_filled["start_delt_next"] = agg_start_filled["start_delt_pos"].shift(-1)
    agg_start_filled = agg_start_filled.fillna(0)

    agg_start_filled["start_peak"] = (agg_start_filled["start_delt_pos"] == True) & (agg_start_filled["start_delt_next"] == False)
    ## END LOCI

    agg_end_filled["not_zero_counts"] = agg_end_filled["counts"] > 0
    agg_end_filled["not_zero_before"] = agg_end_filled["not_zero_counts"].shift(1)
    agg_end_filled["not_zero_after"] = agg_end_filled["not_zero_counts"].shift(-1)

    agg_end_filled = agg_end_filled.fillna(0)

    agg_end_filled["zero_smush"] = (agg_end_filled["not_zero_before"] == True) & (agg_end_filled["not_zero_after"] == True) & (agg_end_filled["not_zero_counts"] == False)
    agg_end_filled["next_count"] = agg_end_filled["counts"].shift(-1)
    agg_end_filled = agg_end_filled.fillna(0)
    agg_end_filled["additive_factor"] = agg_end_filled["next_count"] * agg_end_filled["zero_smush"]

    agg_end_filled["counts"] = agg_end_filled["counts"] + agg_end_filled["additive_factor"]

    agg_end_filled["before"] = agg_end_filled["counts"].shift(1)
    agg_end_filled = agg_end_filled.fillna(0)

    agg_end_filled["delta_start"] = agg_end_filled["counts"] - agg_end_filled["before"]
    agg_end_filled["delta_start_next"] = agg_end_filled["delta_start"].shift(-1)
    agg_end_filled = agg_end_filled.fillna(0)

    agg_end_filled["delta_start_modded"] = agg_end_filled["delta_start"] + ((1/2) * agg_end_filled["delta_start_next"])

    agg_end_filled["start_delt_pos"] = agg_end_filled["delta_start_modded"] > 0
    agg_end_filled["start_delt_next"] = agg_end_filled["start_delt_pos"].shift(-1)
    agg_end_filled = agg_end_filled.fillna(0)

    agg_end_filled["end_peak"] = (agg_end_filled["start_delt_pos"] == True) & (agg_end_filled["start_delt_next"] == False)

    ## Create mini dfs to join
    agg_start_mini = agg_start_filled[["start", "start_peak", "counts"]]

    agg_end_mini = agg_end_filled[["end", "end_peak", "counts"]]

    agg_start_mini["type"] = "start"
    agg_start_mini = agg_start_mini.rename(columns={"start":"loci", "start_peak":"peak"})

    agg_end_mini["type"] = "end"
    agg_end_mini = agg_end_mini.rename(columns={"end":"loci", "end_peak":"peak"})

    peaks_dataframe = pd.concat([agg_start_mini, agg_end_mini])
    peaks_dataframe = peaks_dataframe.sort_values("loci")

    clusters = []

    naming_start_cond = True

    current_id = 1

    # Here we are building the csv file
    # which will track potential peaks
    for row in peaks_dataframe.iterrows():

        row_data = list(row[1])

        type_loci = row_data[3]
        peak_status = row_data[1]
        count_status = row_data[2]

        if peak_status:
            potential_loci_tracker["source_id"].append(key)
            potential_loci_tracker["loci"].append(row_data[0])
            potential_loci_tracker["type"].append(type_loci)
            potential_loci_tracker["counts"].append(count_status)

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
    # Looping through start and end pairs
    # basically each fragments start and end loci
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

        # if the start and end clusters agree
        # we say that they are the same and append a value
        if current_cluster_end == current_cluster_start:
            copy_of_test_graph[key]["cluster"].append(current_cluster_start)
            agg += 1

        # if they do not agree, we input the two clusters that
        # the end and start are apart of
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

    # peaks are tracked...
    if i_1 == 0:

        the_peaker = peaks_dataframe[peaks_dataframe["peak"] == True]

        the_peaker["source"] = key

        i_1 += 1

    else:
        the_peaker_temp = peaks_dataframe[peaks_dataframe["peak"] == True]

        the_peaker_temp["source"] = key

        the_peaker = pd.concat([the_peaker, the_peaker_temp])

# Export peaks
potential_peaks = pd.DataFrame.from_dict(potential_loci_tracker)

# Here we detect if a cluster is a primary one
# Only thing is that ties are not broken
the_peaker["primary_cluster_peak"] = the_peaker.groupby(["source","type","cluster"])["counts"].transform(lambda x: x == x.max())
the_peaker.to_csv("cluster_peak_relationship_table.csv", index = False)
os.system("mv cluster_peak_relationship_table.csv " + out_dir)

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

        # This is kept in for figure generation
        # Performance hit is kind of negligible
        if prefix in m:
            switch = True
        else:
            source_name = m
            ident = j

        # pull relevant data

        cluster_values = np.unique(copy_of_test_graph[source_name]["cluster"])

        for value in cluster_values:
            # We find the indexes such that it is in the current cluster value
            # This allows us to pick out the fragments that associate with each
            # cluster
            # [A,B,C]
            # [1,1,2] => If cluster 1, then we get [A,B] using indices
            # I like this solution. 
            indexes = list(np.where(np.array(copy_of_test_graph[source_name]["cluster"]) == value)[0])

            cluster_ids = []

            # Mark is different start and end peak cluster are used
            # Note: They are still clustered if they have the same start and end clusters
            # i.e. [1.2, 1.2] => Equal
            #      [A, B] => [A,B]
            # We cluster these a little differently
            # The index is taken and fragment id is extracted
            # If there is only one, then the loop ends
            # Otherwise, edges are drawn as follows:
            # [A,B,C] => Fragments
            # [1,1,1]
            # The first node is skipped
            # Then an edge is drawn between the current node (B)
            # And the previous (A)
            # A-B
            # Then the next is done
            # A-B-C

            # If ever the fragment is mentioned again in another source, and we have this set:
            # [A,D,E]
            # [3,3,3]
            # Edges are drawn as follows
            # A - B - C
            #   - D - E
            # And they are now clustered together

            # This is an important functionality because isoforms often have sequences
            # That are very similar
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
            # Actually, the same method is applied here
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
# Now each independent graph is kept
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
plt.savefig("P2_cluster_start_end_disagreement.jpeg", dpi = 500)
plt.clf()

ref_table = pd.DataFrame.from_dict(merged_references)
# No longer need you
#ref_table.to_csv("ref_table.csv", index = False)

os.system("mv ref_table.csv " + out_dir)
os.system("mv P2_cluster_start_end_disagreement.jpeg " + out_dir)

# Merge Tables
joined_data = counts_dataset.join(ref_table.set_index("original_id"), on = "ID")
columns_joined = list(joined_data.columns)
end_index = columns_joined.index("sources")
new_id_index = columns_joined.index("new_id")

counts_table = joined_data.iloc[:,[new_id_index] + list(range(1, end_index))]
sum_table = counts_table.groupby("new_id").agg(sum)
sum_table = sum_table.reset_index()

## == ## Detection of outside mapping ## == ##
num = list(counts_dataset.columns).index("ID")

tsv_to_fasta(out_dir + "/filtered_corrected_counts.csv", "filtered_sequences.fa",num, 0, delim = ",")

bowtie_align_pipeline(indexed, working_dir, out_dir + "/filtered_sequences.fa")

gtf_anti_join(full_annotation, processed_annotation_file, "transcript_id", "anti_joined.gtf")

# move antijoined file
os.system('cd ' + out_dir + ";\
          mv anti_joined.gtf " + working_dir)

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
            featureCounts -a anti_joined.gtf -F 'GTF' -g " + attribute_choice + " -o matches.tsv aligned_sams/*.sam -O -M")

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

# make diagrams showing how peaks were laid out

if make_lots == True:
    os.system("mkdir cluster")

if make_lots == True:
    alt_peaks_df["start_end"] = alt_peaks_df["start_end"].astype("string")

    alt_peaks_df[["start","end"]] = alt_peaks_df["start_end"].str.split(".", expand = True)

    alt_peak_small = alt_peaks_df[["source_id", "start", "end"]]

    for source in list(np.unique(the_peaker["source"])):
        source_specific_df = the_peaker.loc[the_peaker["source"] == source, :]
        source_specific_df = source_specific_df.loc[source_specific_df["primary_cluster_peak"] == True, :]

        source_specific_df = source_specific_df.sort_values("cluster")

        alt_peaks_sources = alt_peak_small.loc[alt_peaks_df["source_id"] == source, :]
        alt_peaks_sources = alt_peaks_sources.drop_duplicates(keep = "first")

        checkpoint = max(list(np.unique(source_specific_df["cluster"])))

        starts_alt = list(alt_peaks_sources["start"])
        ends_alt = list(alt_peaks_sources["end"])

        starts_alt = [int(i) for i in starts_alt]
        ends_alt = [int(i) for i in ends_alt]

        alt_starts = []
        alt_ends = []

        

        for n in list(np.unique(source_specific_df["cluster"])):
            x_constructor = list(source_specific_df.loc[source_specific_df["cluster"] == n, :]["loci"])
            y_constructor = [i] * len(x_constructor)

            plt.plot(x_constructor, y_constructor, label = "Standard")
            plt.title(source)
            plt.xlabel("Position (bp)")
            plt.ylabel("cluster #")
            plt.fill_between(x_constructor, y_constructor, alpha = 0.3)

            if i in starts_alt:
                alt_starts.append(min(x_constructor))
            if i in ends_alt:
                alt_ends.append(max(x_constructor))

            i += 1

        k = 0
        for num in alt_starts:
            plt.plot([num] + [alt_ends[k]], [i] * len([num] + [alt_ends[k]]), label = "Alternate")
            plt.title(source)
            plt.xlabel("Position (bp)")
            plt.ylabel("cluster #")
            plt.fill_between([num] + [alt_ends[k]], [i] * len([num] + [alt_ends[k]]), alpha = 0.1)

            k += 1
            i += 1


        plt.legend()
        plt.savefig("cluster/cluster_loc_" + source + ".jpeg")
        plt.clf()

        i = 1

# The final component of this script is to call sets of license plates
# that associate with each fragment
# this will allow for some cool stuff to be done
# once other species are investigated
# i think the mintmap code can be modified to include
# some degree of similarity
# Also, will include if it is flagged as an out of 
# biotype mapping (just identified i guess)

sum_table
ref_table

col_constructor_orig = []
col_constructor_outside = []
for row in sum_table.iterrows():
    # each row will have a unique value
    row_data = list(row[1])

    # get the merged id
    merged_id = row_data[0]

    # get all rows that have 
    orig_ids = annotated_ref_table.loc[annotated_ref_table["new_id"] == merged_id,"original_id"]
    orig_ids_list = list(orig_ids)

    # outside maps
    outside_maps = annotated_ref_table.loc[annotated_ref_table["new_id"] == merged_id, "annotation"]
    outside_maps = set(outside_maps)

    # split up based on - so we have a list of license plates
    license_plates = [n.split(sep = "-")[-2:] for n in orig_ids_list]

    # make into one string
    license_plates_clean = ["-".join(n) for n in license_plates]

    col_constructor_orig.append(license_plates_clean)
    col_constructor_outside.append(outside_maps)

    # we will add this list to the sum table and export

sum_table_ann = sum_table
sum_table_ann["original_id_set"] = col_constructor_orig
sum_table_ann["outside_maps"] = col_constructor_outside

sum_table_ann.to_csv("merged_counts.csv", index = False)

os.system("mv merged_counts.csv " + out_dir)