# NuMap
NuMap - http://www-hsc.usc.edu/~valouev/NuMap/NuMap.html by Anton Valouev

I uploaded this software because I couldn't compile it from the original distribution (there was a header missing), so I thought it would be of interest to others to get pass this problem as well.

I am further uploading the compiled (linux x86_64) binaries as well (they're in `bin`).

## Install
    git clone git@github.com:afrendeiro/NuMap.git
    cd NuMap
    make
    # add to $PATH or copy to /usr/bin/

## Workflow
From the original documentation:

Below is the description of the minimal analysis necessary to perform nucleosome map comparison. Details of other additional modules are provided in the 'Module Description' section:

1. Obtain SAM MNase-Seq files.
2. Calculate align files:

    `sam_2_align_bin genome_table=hg19.gt input_sam_file=MNase_seq.sam analysis_path=./MNase_seq_analysis`
    
    Note:hg19.gt is provided with the distribution. You can use QuEST software to generate genome tables for other genomes.
3. Estimate MNase-Seq library fragment size

    `dist_plots analysis_path=./MNase_seq_analysis`

4. Calculate dyads from MNase-seq reads

    `align_2_dyads analysis_path=./MNase_seq_analysis`

5. Calculate stringency plots
    
    `dyad_stringency analysis_path=./MNase_seq_analysis`
6. Call dyads
    
    `call_dyads analysis_path=./MNase_seq_analysis output_path=././MNase_seq_analysis/dyad_calls/`

7. Match dyads between two experiments
    
    `match_dyads called_dyads_file1=./MNase_seq_analysis1/dyad_calls/dyad_positions.txt called_dyads_file2=././MNase_seq_analysis2/dyad_calls/dyad_positions.txt match_dist=30 genome_table=hg19.gt output_path=./dyad_matching`

8. Calculate nucleosome binding differences at a set of genomic sites
    
    `diff_summaries_at_sites called_dyads_file1=./dyad_matching/dyad_matches unmatched_file1=./dyad_matching/unmatched_dyads1 unmatched_file2=./dyad_matching/unmatched_dyads2 positions_file=sites.txt genome_table=hg19.gt output_file=sites.nucl_diff dist=5000`

## NuMap Modules
From the original documentation:

### sam_2_align_bin
Converts sam files to base-by-base alignment binary profiles accross the genome. Example:

    sam_2_align_bin genome_table=hg19.gt input_sam_file=MNase_Seq.sam analysis_path=./MNase_Seq_analysis

genome table is a tab-delimited file providing coordinate sizes:
chr1   249250621
chr2   243199373
...

### dist_plots
Calculates distograms and phasograms of raw data. Estimates average fragment size of MNase-Seq library. Example:

    dist_plots analysis_path=./MNase_Seq_analysis

### align_2_dyads
Calculate dyad coordinates by adjusting coordinates of reads ends by half the library size. Example:

    align_2_dyads analysis_path=./MNase_Seq_analysis

### simulate_dyads
Simulate dyad coordinates by permuting MNase-Seq inferred dyads within a specified window. Example:

    simulate_dyads analysis_path=./MNase_Seq_analysis window=1000

### dyad_stringency
Calculate dyad stringency using dyad coordinates. Example:

    dyad_stringency analysis_path=./MNase_Seq_analysis mock=no core_size=147 bw=100

- mock can be set to yes if stringency profiles to be generated from the mock dyads
- core_size sets the 'infringement' window
bw is a triweight bandwidth for the stringency formula

### call_dyads
Call dyads using dyad stringency profiles. Example:

    call_dyads analysis_paht=./MNase_Seq_analysis output_path=./MNase_Seq_analysis/dyad_calls mock=no peak_value=mean drop=0.1 drop_window=150 bw=100

- mock value can be set to yes to call dyads from the mock data
- peak_value controls the peak threshold of the dyad stringency profile at which positionined nucleosomes are called (uses mean of the chromosome profile by default)
- drop controls the degree of drop to the left and to the right of the peak (0.1 default)
- drop_window controls the size of the window within the drop should be observed
- bw is the triweight bandwidth value necessary to estimate the smoothed dyad count at the peak position, which is reported along with the dyad call coordinates

This module will produce dyad call coordinates in ./MNase_Seq_analysis/dyad_calls/dyad_positions.txt
The file is tab delimited and has the following fields:

    #chrom  pos     stringency      count_enrichment        count
    chr1    10099   0.3758  0.7594  51.2832

- chrom and pos represent coordinate of the called dyad
- stringency represents nucleosome positioning stringency of that call
- count_enrichment represents the enrichment of the smoothed dyad frequencies based on chromomsome-wide normalization
- count represents the smoothed count using triweight smoothing

This module also generates bed files containing dyad calls which can be visualized using UCSC Genome Browser. The files are located under the same directory as dyad calls and are split by chrom

### nucleosome_coverage
Calculates nucleosome coverage profiles. Example:

    nucleosome_coverage analysis_path=./MNase_Seq_analysis/ mock=no core_size=147 spacing=193

- mock can be set to 'yes' to providing coverage profiles from the mock data
- core_size provides the size of the nucleosome footrpint for the coverage calculations
- spacing provides the inter-nucleosome spacing in order to meaure the relative coverage

### coverage_to_wig
Produces coverage wig files to be viewed using UCSC Genome Browser. Example:

    coverage_to_wig analysis_path=./MNase_Seq_analysis positions_file=sites.txt output_file=sites_coverage.wig max_dist=10000 track_name=dyad_coverage

- positions_file represents sites near which to display the wig plots. The files are tab delimited and have the following fileds:
<crhom> <pos>
- max_dist parameter provides the maximum distance for plotting the wig curves from the specified sites
- track_name parameter can be used to define the UCSC track name value

### stringency_to_wig
Produces stringency wig files to be viewed using UCSC Genome Browser. Example:

    stringency_wig analysis_path=./MNase_Seq_analysis positions_file=sites.txt output_file=sites_stringency.wig max_dist=1000 track_name=dyad_stringency

- positions_file represents sites near which to display the wig plots. The files are tab delimited and have the following fileds:
- max_dist parameter provides the maximum distance for plotting the wig curves from the specified sites
- track_name parameter can be used to define the UCSC track name value

### phasogram_of_sites
Calculates phasogram of sites. Dyad calls can be used as an input to calcualte phasogram of the dyad calls. Example:

    phasogram_of_sites positions_file=./MNase_Seq_analysis/dyad_calls/dyad_positions.txt output_file=./MNase_Seq_analysis/dyad_calls/phasogram max_dist=3000

- max_dist provides the maximum distance of the phasogram
- the output file will have the following format
<dist> <count>
for a distance range [0,max_dist]

### estimate_peaks
Provides smoothing for any plot including phasogram and dyadogram using bandwidth of the specified size. Can be used to estimate nucleosome spacing. Example:

    estimate_peaks dist_file=./MNase_Seq_analysis/dyad_calls/phasogram output_file=./MNase_Seq_analysis/dyad_calls/phasogram.smoothed bw=20 field=1

The output file will contain the smoothed positions in the third field. The output file `./MNase_Seq_analysis/dyad_calls/phasogram.smoothed.peak_pos.txt` will contain peaks of the phasogram along to be used for fitting in R and estimating nucleosome spacing. The file has the following format:

    <peak number> <peak positions>

The file can be opened in R to fit the linear regression:

    w<-read.table(file="./MNase_Seq_analysis/dyad_calls/phasogram.smoothed.peak_pos.txt")
    mlm<-lm(w$V2~w$V1)
    mlm

The following output should be seen:

    Call:
    lm(formula = w$V2 ~ w$V1)
    Coefficients:
    (Intercept)         w$V1  
          193.7        193.6

- bw provides the smoothing bandwidth for the triweight kernel
- field provides the field number of the input file to be smoothed (0-base count)

### called_dyad_organization
Plots a dyadogram of called dyads at the collection of sites (e.g. TFBS). The resulting profiles can be further smoothed using estimate peaks module. Example:

    called_dyad_organization positions_file=sites.txt called_dyads_file=./MNase_Seq_analysis/dyad_calls/dyad_positions.txt genome_table=hg19.gt output_file=sites.dyadogram max_dist=3000

- positions_file provides coordinates of sites for which to plot dyadogram. The file should be tab-delimited with the following fields:

    <chrom>	 <coord>

- called_dyads_file provides coordinates of dyad calls based on which the dyadogram is obtained
- max_dist provides the maximum distance for the dyadogram calculation

### match_dyads
Matches the dyads between two nucleosome mapping experiments. Example:

    match_dyads called_dyads_file1=dyad_calls1.txt called_dyads_file2=dyad_calls2.txt genome_table=hg19.gt match_dist=30 output_path=./dyad_matching/

- called_dyads_file1 and called_dyads_file2 provide files containing dyad calls (see above)
- match_dist provides the matching distance

The program outputs the following files of interest within the ouptut path:

- dyad_matches provides coordinates of matched dyads. The file is tab-delimited with the following fields:

    <chrom> <pos1> <pos2> <str1> <str2>

where pos1 and pos2 provide positoins of matched dyads on the chromosome chrom, and str1 and str2 provide their corresponding positioning strengths
- unmatched_dyads1 and unmatched_dyads2 provide the coordinates of dyads that were not matched.
- matching_dyads.<chrom>.bed provide bed files for visualizing dyad matches using UCSC Genome Browser

### permute_dyads
Permute called dyad coordinates to generate 'random' dyad_positions. Example:

    ./permute_dyads called_dyads_file=./MNase_Seq_analysis/dyad_calls/dyad_positions.txt genome_table=hg19.gt output_file=./MNase_Seq_analysis/dyad_calls/mock.dyad_positions.txt window=1000

- window provides a window size in which reference dyad call positions are permuted to achieve randomization

### diff_summaries_at_sites
Calculate nucleosome binding differences at a collection of sites (e.g. TF binding sites). The resulting file can be then smoothed using estimate_peaks module. Example:

    diff_summaries_at_sites dyad_match_file=./dyad_matching/dyad_matches unmatched_file1=./dyad_matching/unmatched_dyads1 unmatched_file2=./dyad_matching/unmatched_dyads2 positions=sites.txt genome_table=hg19.gt output_file=sites.chrom_diff dist=4000

- dist provides the distanace across which to plot the nucleosome binding differences

The output file is tab-delimited and has the following format:

    <offset> <match counts> <umnatch_counts>

umnatch_counts / (2*match_counts) provides a nucleosome binding difference
