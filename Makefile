CC = g++

include config.mk

default: blah

all: blah

SRC_DIR = ./src/

LOCAL_DEPENDENCIES = params genome_table string_utils sam_alignment psam_alignment_file math_functions file_utilities kernels

LOCAL_TARGETS = sam_2_align_bin dist_plots dyads_vs_TFBS dyad_stringency nucleosome_coverage stringency_to_wig coverage_to_wig align_2_dyads call_dyads simulate_dyads estimate_peaks called_dyad_organization compare_dyad_coordinates match_dyads permute_dyads phasogram_of_sites diff_summaries_at_sites dyad_call_qvalues calculate_tags_for_regions dyad_match_qvalues call_enrichment_peaks

LOCAL_DEP_OBJECTS = $(LOCAL_DEPENDENCIES:=.o)
LOCAL_TARGET_OBJECTS = $(LOCAL_TARGETS:=.o)

local_dep_objects:
	@echo
	@echo =================================
	@echo =  Building common dependencies =
	@echo =================================
	@echo
	cd $(SRC_DIR) && $(MAKE) $(LOCAL_DEP_OBJECTS)

local_target_objects:
	@echo
	@echo ================================
	@echo =     Building target objs     =
	@echo ================================
	@echo
	cd $(SRC_DIR) && $(MAKE) $(LOCAL_TARGET_OBJECTS)

local_targets:
	@echo
	@echo ================================
	@echo =     Building executables     =
	@echo ================================
	@echo
	make $(LOCAL_TARGETS)

default: blah

blah: local_target_objects local_dep_objects local_targets

clean:
	rm -rf $(SRC_DIR)*.o
	rm $(LOCAL_TARGETS)

sam_2_align_bin_deps = string_utils params genome_table sam_alignment psam_alignment_file sam_2_align_bin
sam_2_align_bin_objs = $(addprefix $(SRC_DIR), $(sam_2_align_bin_deps:=.o))

sam_2_align_bin:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(sam_2_align_bin_objs) -o sam_2_align_bin

dist_plots_deps = string_utils params genome_table math_functions dist_plots
dist_plots_objs = $(addprefix $(SRC_DIR), $(dist_plots_deps:=.o))

dist_plots:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(dist_plots_objs) -o dist_plots

dyads_vs_TFBS_deps = string_utils params genome_table math_functions dyads_vs_TFBS
dyads_vs_TFBS_objs = $(addprefix $(SRC_DIR), $(dyads_vs_TFBS_deps:=.o))

dyads_vs_TFBS:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(dyads_vs_TFBS_objs) -o dyads_vs_TFBS

dyad_stringency_deps = string_utils params genome_table math_functions dyad_stringency
dyad_stringency_objs = $(addprefix $(SRC_DIR), $(dyad_stringency_deps:=.o))

dyad_stringency:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(dyad_stringency_objs) -o dyad_stringency

nucleosome_coverage_deps = string_utils params genome_table math_functions file_utilities nucleosome_coverage
nucleosome_coverage_objs = $(addprefix $(SRC_DIR), $(nucleosome_coverage_deps:=.o))

nucleosome_coverage:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(nucleosome_coverage_objs) -o nucleosome_coverage

stringency_to_wig_deps = string_utils params genome_table math_functions stringency_to_wig
stringency_to_wig_objs = $(addprefix $(SRC_DIR), $(stringency_to_wig_deps:=.o))

stringency_to_wig:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(stringency_to_wig_objs) -o stringency_to_wig

coverage_to_wig_deps = string_utils params genome_table math_functions coverage_to_wig
coverage_to_wig_objs = $(addprefix $(SRC_DIR), $(coverage_to_wig_deps:=.o))

coverage_to_wig:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(coverage_to_wig_objs) -o coverage_to_wig

align_2_dyads_deps = string_utils params genome_table math_functions file_utilities align_2_dyads
align_2_dyads_objs = $(addprefix $(SRC_DIR), $(align_2_dyads_deps:=.o))

align_2_dyads:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(align_2_dyads_objs) -o align_2_dyads

call_dyads_deps = string_utils params genome_table math_functions file_utilities call_dyads
call_dyads_objs = $(addprefix $(SRC_DIR), $(call_dyads_deps:=.o))

call_dyads:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(call_dyads_objs) -o call_dyads

simulate_dyads_deps = string_utils params genome_table math_functions file_utilities simulate_dyads
simulate_dyads_objs = $(addprefix $(SRC_DIR), $(simulate_dyads_deps:=.o))

simulate_dyads:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(simulate_dyads_objs) -o simulate_dyads

estimate_peaks_deps = string_utils params genome_table math_functions file_utilities estimate_peaks
estimate_peaks_objs = $(addprefix $(SRC_DIR), $(estimate_peaks_deps:=.o))

estimate_peaks:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(estimate_peaks_objs) -o estimate_peaks

called_dyad_organization_deps = string_utils params genome_table math_functions file_utilities called_dyad_organization
called_dyad_organization_objs = $(addprefix $(SRC_DIR), $(called_dyad_organization_deps:=.o))

called_dyad_organization:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(called_dyad_organization_objs) -o called_dyad_organization

compare_dyad_coordinates_deps = string_utils params genome_table math_functions file_utilities compare_dyad_coordinates
compare_dyad_coordinates_objs = $(addprefix $(SRC_DIR), $(compare_dyad_coordinates_deps:=.o))

compare_dyad_coordinates:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(compare_dyad_coordinates_objs) -o compare_dyad_coordinates

match_dyads_deps = string_utils params genome_table math_functions file_utilities match_dyads
match_dyads_objs = $(addprefix $(SRC_DIR), $(match_dyads_deps:=.o))

match_dyads:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(match_dyads_objs) -o match_dyads

permute_dyads_deps = string_utils params genome_table math_functions file_utilities permute_dyads
permute_dyads_objs = $(addprefix $(SRC_DIR), $(permute_dyads_deps:=.o))

permute_dyads:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(permute_dyads_objs) -o permute_dyads

phasogram_of_sites_deps = string_utils params genome_table math_functions phasogram_of_sites
phasogram_of_sites_objs = $(addprefix $(SRC_DIR), $(phasogram_of_sites_deps:=.o))

phasogram_of_sites:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(phasogram_of_sites_objs) -o phasogram_of_sites

diff_summaries_at_sites_deps = string_utils params genome_table math_functions kernels diff_summaries_at_sites
diff_summaries_at_sites_objs = $(addprefix $(SRC_DIR), $(diff_summaries_at_sites_deps:=.o))

diff_summaries_at_sites:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(diff_summaries_at_sites_objs) -o diff_summaries_at_sites

dyad_call_qvalues_deps = string_utils params genome_table dyad_call_qvalues
dyad_call_qvalues_objs = $(addprefix $(SRC_DIR), $(dyad_call_qvalues_deps:=.o))

dyad_call_qvalues:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(dyad_call_qvalues_objs) -o dyad_call_qvalues

calculate_tags_for_regions_deps = string_utils params genome_table calculate_tags_for_regions
calculate_tags_for_regions_objs = $(addprefix $(SRC_DIR), $(calculate_tags_for_regions_deps:=.o))

calculate_tags_for_regions:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(calculate_tags_for_regions_objs) -o calculate_tags_for_regions

dyad_match_qvalues_deps = string_utils params genome_table dyad_match_qvalues
dyad_match_qvalues_objs = $(addprefix $(SRC_DIR), $(dyad_match_qvalues_deps:=.o))

dyad_match_qvalues:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(dyad_match_qvalues_objs) -o dyad_match_qvalues

call_enrichment_peaks_deps = string_utils params genome_table file_utilities math_functions call_enrichment_peaks
call_enrichment_peaks_objs = $(addprefix $(SRC_DIR), $(call_enrichment_peaks_deps:=.o))

call_enrichment_peaks:
	$(CC) -I$(SRC_DIR) $(CCFLAGS) $(call_enrichment_peaks_objs) -o call_enrichment_peaks

