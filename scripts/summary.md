# sRNAfrag Summary Report

## P1 Module Outputs

No figures are produced during the first module. Primary function is to align and obtain counts. 

### Directory Structure

Main-Dir

|__out

    |__tables
        |__num_reads_bam.csv // Number of reads for each library. Can be used for CPM/TPM/RPKM normalization.
        |__sRNA_frag_counts.csv // Counts for all fragments. Raw data. 

## S1 Module Outputs

Basic analyses and information regarding transcripts. Prepares for P2 module. 

### Count Distribution and Correction
<img src="figures/S1_correction_density.jpeg" alt="drawing" width="800"/>

Figure 1. Pre and post correction counts with linear fits. Count density is also depicted. 

### Standardized Loci Location
<img src="figures/S1_std_loci.jpeg" alt="drawing" width="800"/>

Figure 2. Start loci are standardized onto a 0-1 scale. If start loci do not agree (due to a fragment possibly being derived from multiple sources) its average is taken and the standard deviation is plotted. The average standardized loci is plotted by frequency. Finally, the count distribution is binned and plotted using a bar graph.

### Motifs by Frequency
<img src="figures/S1_freq_motifs.jpeg" alt="drawing" width="800"/>

Figure 3. Logos of most common motifs by frequency.

### Motifs by Counts
<img src="figures/S1_motifs_counts.jpeg" alt="drawing" width="800"/>

Figure 4. Logos of most common motifs by counts.

### Directory Structure

Main-Dir

|__out

    |__tables
    |   |__num_reads_bam.csv // Number of reads for each library, Can be used for CPM/TPM/RPKM normalization
    |   |__sRNA_frag_counts.csv // Counts for all fragments, Raw data
    |   |__filtered_corrected_counts.csv // Adjusted counts
    |   |__filtered_counts.csv // Unadjusted counts
    |
    |__figures
        |__S1_correction_density.jpeg // Figure depicting adjusted and count distribution
        |__S1_freq_motifs.jpeg // Motifs by frequency
        |__S1_motifs_counts.jpeg // Motifs by counts
        |__S1_std_loci.jpeg // Loci location information
        |__S1_five-mer_3p.html // Most common 3p fiver-mers | Please feel free to modify the S1_figures.R
        |__S1_five-mer_5p.html // Most common 5p fiver-mers | If you'd like to generate different length.

## P2 Module Outputs

### Number of Sources for Each Fragment
<img src="figures/P2_Filter_Passing_Hist.jpeg" alt="drawing" width="800"/>

Figure 7. Distribution of number of potential sources each fragment has. 

### Number of Transcripts that span two clusters
<img src="figures/P2_cluster_start_end_disagreement.jpeg" alt="drawing" width="800"/>

Figure 8. Ratio of disagreements (i.e. start = cluster 1, end = cluster 2) over total mentions. May signify interesting fragmentation pattern. 

### Directory Structure

Main-Dir

|__out

    |__tables
    |   |__num_reads_bam.csv // Number of reads for each library, Can be used for CPM/TPM/RPKM normalization
    |   |__sRNA_frag_counts.csv // Counts for all fragments, Raw data
    |   |__filtered_corrected_counts.csv // Adjusted counts
    |   |__filtered_counts.csv // Unadjusted counts
    |   |__cluster_peak_relationship_table.csv // Associates each start and end loci with a count and if it is a primary peak.
    |   |__alternate_peaks.csv // Marks peaks that span two peaks. Used for plots in cluster folder.
    |   |__annotated_ref_table.csv // Reference table associating fragment ID with merged ID. Includes information about external maps.
    |   |__merged_counts.csv // The final output with sample counts and the set of license plates and external maps. 
    |
    |__figures
    |   |__S1_correction_density.jpeg // Figure depicting adjusted and count distribution
    |   |__S1_freq_motifs.jpeg // Motifs by frequency
    |   |__S1_motifs_counts.jpeg // Motifs by counts
    |   |__S1_std_loci.jpeg // Loci location information
    |   |__S1_five-mer_3p.png // Most common 3p fiver-mers | Please feel free to modify the S1_figures.R
    |   |__S1_five-mer_5p.png // Most common 5p fiver-mers | If you'd like to generate different length.
    |   |__P2_Filter_Passing_Hist.jpeg // The number of sources that fragments could have. Distribution.
    |   |__P2_cluster_start_end_disagreement.jpeg // Distribution of situations where start and end clusters are not equal. Not necessarily bad. Potentially interesting.
    |
    |__sequences
    |   |__filtered_sequences.fa // Sequences in the fasta format. Do Gibbs Sampling to find motif. Or do target prediction.
    |
    |__source_peaks // If makelots is true, each source will have counts plotted against loci position.
    |   |__source1
    |   |__source2
    |   |__...
    |   |__source_n
    |
    |__cluster // If makelots is true, each source will have the cluster zones against loci highlighted, along with disagreements plotted.
    |   |__source1
    |   |__source2
    |   |__...
    |   |__source_n
    |
    |__int_files // Intermediate script files pasted into the directory. Please take a look if you'd like to contribute on the github repo!
        |__scripts...
        |__extra_files...