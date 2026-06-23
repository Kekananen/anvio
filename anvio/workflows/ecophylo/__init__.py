""" Classes to define and work with anvi'o ecophylo workflows. """

import os
import argparse
import pandas as pd

import anvio
import anvio.utils as u
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths

with terminal.SuppressAllOutput():
    import anvio.data.hmm

from anvio.errors import ConfigError
from anvio.workflows import WorkflowSuperClass
from anvio.workflows.read_recruitment import ReadRecruitmentModule
from anvio.genomedescriptions import GenomeDescriptions
from anvio.genomedescriptions import MetagenomeDescriptions
from anvio.artifacts.samples_txt import SamplesTxt

import anvio.db as db
import anvio.constants as constants

__copyright__ = "Copyleft 2015-2024, The Anvi'o Project (http://anvio.org/)"
__credits__ = ['mschecht']
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Matthew S. Schechter"
__email__ = "mschechter@uchicago.edu"


run = terminal.Run()

class EcoPhyloWorkflow(ReadRecruitmentModule, WorkflowSuperClass):
    def __init__(self, args=None, run=terminal.Run(), progress=terminal.Progress()):
        self.init_workflow_super_class(args, workflow_name='ecophylo')

        # initialize the read recruitment module (adds bowtie_build, bowtie, etc.)
        ReadRecruitmentModule.__init__(self)

        # reps contigs DB has short ribosomal protein sequences (< 600 bp).
        # default min-contig-length of 1000 filters them all out.
        self.default_config['anvi_profile']['--min-contig-length'] = 0

        # Snakemake rules
        self.rules.extend(['extract_hmm_hit_seqs',
                           'anvi_run_scg_taxonomy',
                           'cat_hmm_hit_seqs',
                            'hmmsearch_combined',
                            'filter_hmm_hits_combined',
                            'filter_hmm_hits_sample',
                            'process_hmm_hits',
                           'combine_sequence_data',

                           'cluster_X_percent_sim_mmseqs',
                           'anvi_profile_blitz',
                           'subset_AA_seqs_with_mmseqs_reps',
                           'subset_AA_seqs_with_coverage_reps',
                           'align_sequences',
                           'trim_alignment',
                           'remove_sequences_with_X_percent_gaps',
                            'count_num_sequences_filtered',
                            'fasttree',
                            'iqtree',
                            'anvi_summarize',
                            'rename_tree_tips',
                           'make_misc_data',
                           'add_misc_data_to_taxonomy',
                           'make_anvio_state_file',
                            'anvi_import_state',
                            'extract_QCd_sequence_headers',
                            'build_rep_external_gene_calls',
                            'anvi_gen_contigs_database_reps',
                            'anvi_run_scg_taxonomy_reps',
                            'anvi_estimate_scg_taxonomy_reps'
                            ])


    def init(self):
        """This function is called from within the Snakefile to initialize parameters."""

        # Remove inherited keys from parent modules so we control insertion order
        # in 00-07 sequence for get_default_config().
        for k in ["QC_DIR", "MAPPING_DIR", "PROFILE_DIR", "MERGE_DIR"]:
            self.dirs_dict.pop(k, None)

        # Set all 8 canonical keys in 00-07 order
        self.dirs_dict["LOGS_DIR"] = "00_LOGS"
        self.dirs_dict["HMM_HITS_DIR"] = "01_HMM_HITS"
        self.dirs_dict["REPRESENTATIVES_DIR"] = "02_REPRESENTATIVES"
        self.dirs_dict["CONTIGS_DIR"] = "03_CONTIGS"
        self.dirs_dict["PHYLO"] = "04_TREE"
        self.dirs_dict["MAPPING_DIR"] = "05_MAPPING"
        self.dirs_dict["PROFILE_DIR"] = "06_ANVIO_PROFILE"
        self.dirs_dict["MERGE_DIR"] = "07_RESULTS"
        self.dirs_dict["POOLED_HMM_DIR"] = "01_HMM_HITS/POOLED_HMM"

        super().init()
        # QC_DIR is inherited from ReadRecruitmentModule but kept at runtime.
        # It is excluded from the default config via get_default_config() override.
        self.dirs_dict["QC_DIR"] = "05_MAPPING"

        # Make log directories
        if not os.path.exists(self.dirs_dict['LOGS_DIR']):
            os.makedirs(self.dirs_dict['LOGS_DIR'])

        self.names_list = []
        self.contigs_db_name_path_dict = {}
        self.contigs_db_name_bam_dict = {}

        # Load input files
        self.metagenomes = self.get_param_value_from_config(['metagenomes'])
        self.external_genomes = self.get_param_value_from_config(['external_genomes'])
        self.hmm_list_path = self.get_param_value_from_config(['hmm_list'])
        self.samples_txt_file = self.get_param_value_from_config(['samples_txt'])
        self.run_genomes_sanity_check = self.get_param_value_from_config(['run_genomes_sanity_check'])

        if not self.metagenomes and not self.external_genomes:
            raise ConfigError('Please provide at least a metagenomes.txt or external-genomes.txt in your '
                              'EcoPhylo config file.')

        if not self.hmm_list_path:
            raise ConfigError('Please provide a path to an hmm_list.txt')

        self.run_scg_taxonomy = self.get_param_value_from_config(['anvi_run_scg_taxonomy', 'run'])

        self.init_hmm_list_txt()

        gene_caller_to_use = self.get_param_value_from_config(['gene_caller_to_use'])
        if not gene_caller_to_use:
            gene_caller_to_use = constants.default_gene_callers[-1]

        sanity_checked_metagenomes_file = os.path.join(self.dirs_dict['LOGS_DIR'], "sanity_checked_metagenomes.txt")
        sanity_checked_genomes_file = os.path.join(self.dirs_dict['LOGS_DIR'], "sanity_checked_genomes.txt")

        if self.metagenomes:
            filesnpaths.is_file_exists(self.metagenomes)
            self.metagenomes_df = pd.read_csv(self.metagenomes, sep='\t', index_col=False)

            if self.run_genomes_sanity_check:
                if not os.path.exists(sanity_checked_metagenomes_file):
                    args = argparse.Namespace(metagenomes=self.get_param_value_from_config(['metagenomes']), gene_caller = gene_caller_to_use)
                    g = MetagenomeDescriptions(args)
                    g.load_metagenome_descriptions(init=False)
                    self.metagenomes_dict = g.metagenomes_dict
                    self.metagenomes_name_list = list(self.metagenomes_dict.keys())
                    self.metagenomes_path_list = [value['contigs_db_path'] for key,value in self.metagenomes_dict.items()]

                    with open(sanity_checked_metagenomes_file, 'w') as fp:
                        pass
                else:
                    self.run.warning("You have declared run_genomes_sanity_check == false. anvi'o takes no responsibility "
                                     "for any genomes or metagenomes that cause issues downstream in ecophylo.")
                    self.metagenomes_name_list = self.metagenomes_df.name.to_list()
                    self.metagenomes_path_list = self.metagenomes_df.contigs_db_path.to_list()
            else:
                self.metagenomes_name_list = self.metagenomes_df.name.to_list()
                self.metagenomes_path_list = self.metagenomes_df.contigs_db_path.to_list()

            self.contigs_db_name_path_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_path_list)))

            if 'bam' in self.metagenomes_df.columns:
                self.contigs_db_name_bam_dict.update(dict(zip(self.metagenomes_name_list, self.metagenomes_df.bam)))
                self.metagenomes_profiles_list = self.metagenomes_df.bam.to_list()

            self.names_list.extend(self.metagenomes_name_list)

        else:
            self.metagenomes_name_list = []

        if self.external_genomes:
            filesnpaths.is_file_exists(self.external_genomes)
            self.external_genomes_df = pd.read_csv(self.external_genomes, sep='\t', index_col=False)

            if self.run_genomes_sanity_check:
                if not os.path.exists(sanity_checked_genomes_file):
                    # FIXME: metagenomes.txt or external-genomes.txt with multiple gene-callers will break
                    # here. Users should only have one type of gene-caller e.g. "NCBI_PGAP".

                    args = argparse.Namespace(external_genomes=self.external_genomes, gene_caller = gene_caller_to_use)
                    genome_descriptions = GenomeDescriptions(args)
                    genome_descriptions.load_genomes_descriptions(init=False)
                    self.external_genomes_dict = genome_descriptions.external_genomes_dict
                    self.external_genomes_names_list = list(self.external_genomes_dict.keys())
                    self.external_genomes_path_list = [value['contigs_db_path'] for key,value in self.external_genomes_dict.items()]

                    with open(sanity_checked_genomes_file, 'w') as fp:
                        pass
                else:
                    self.external_genomes_names_list = self.external_genomes_df.name.to_list()
                    self.external_genomes_path_list = self.external_genomes_df.contigs_db_path.to_list()
            else:
                self.external_genomes_names_list = self.external_genomes_df.name.to_list()
                self.external_genomes_path_list = self.external_genomes_df.contigs_db_path.to_list()

            self.contigs_db_name_path_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_path_list)))

            if 'bam' in self.external_genomes_df.columns:
                self.contigs_db_name_bam_dict.update(dict(zip(self.external_genomes_names_list, self.external_genomes_df.bam)))
                self.external_genomes_profiles_list = self.external_genomes_df.bam.to_list()

            self.names_list.extend(self.external_genomes_names_list)

        else:
            self.external_genomes_names_list = []

        # Pre-compute which (sample, hmm_source) pairs already have HMMs in their contigs DB.
        # This is done once at init(), before any Snakemake resolution, to avoid confusion
        # if anvi-run-hmms modifies the contigs DB during execution (which would make a
        # Path B sample look like Path A on re-resolution).
        self.hmm_source_presence = {}
        all_sources = set(value['source'] for value in self.hmm_dict.values())
        for sample_name in self.names_list:
            contigs_db_path = self.contigs_db_name_path_dict[sample_name]
            database = db.DB(contigs_db_path, None, ignore_version=True)
            sources_in_db = set(database.get_table_as_dict('hmm_hits_info').keys())
            database.disconnect()
            for source in all_sources:
                self.hmm_source_presence[(sample_name, source)] = source in sources_in_db

        # Pre-compute Path A/B decision for each (sample, source).
        # Path A: HMMs pre-exist in contigs DB AND no non-empty domtblout
        #         from a prior anvi-run-hmms run (i.e. this is a fresh run
        #         on samples where anvi-run-hmms was already run externally).
        # Path B: HMMs must be computed by anvi-run-hmms within the workflow,
        #         or a previous run left a non-empty domtblout (rerun stability).
        #         In Path B, survivors are bare gene caller IDs; in Path A,
        #         survivors have full deflines matching the combined hmmsearch.
        self.path_is_a = {}
        for sample_name in self.names_list:
            for source in all_sources:
                domtblout_path = os.path.join(
                    self.dirs_dict['HMM_HITS_DIR'], sample_name,
                    f"{source}-dom-hmmsearch", "hmm.domtable",
                )
                real_domtblout = os.path.exists(domtblout_path) and os.path.getsize(domtblout_path) > 0
                hmm_present = self.hmm_source_presence.get((sample_name, source), False)
                self.path_is_a[(sample_name, source)] = hmm_present and not real_domtblout

        # Make variables that tells whether we have metagenomes.txt, external-genomes.txt, or both
        if self.metagenomes and not self.external_genomes:
            self.mode = 'metagenomes'
        if not self.metagenomes and self.external_genomes:
            self.mode = 'external_genomes'
        if self.metagenomes and self.external_genomes:
            self.mode = 'both'

        self.AA_mode = self.get_param_value_from_config(['cluster_X_percent_sim_mmseqs', 'AA_mode'])

        if self.samples_txt_file:
            # we initialize the samples.txt to run the sanity check before the workflow reaches the
            # metagenomics workflow rule.
            self.samples_txt = SamplesTxt(self.samples_txt_file, expected_format="free")

            self.sample_names_for_mapping_list = self.samples_txt.samples()

            if self.AA_mode == True:
                raise ConfigError("You provided a samples.txt so you're in profile mode! Please change AA_mode to false.")

        else:
            self.run.warning("Since you did not provide a samples.txt, EcoPhylo will assume you do not want "
                             "to profile the ecology of your proteins and will just be making trees for now!")

        # Populate ReadRecruitmentModule interface attributes
        if self.samples_txt_file:
            self.readsets = [
                {
                    'id': s,
                    'type': 'SR',
                    'reads': self.samples_txt.get_sample(s),
                    'base_sample': s,
                }
                for s in self.sample_names_for_mapping_list
            ]
        else:
            self.readsets = []

        self.group_names = list(set(value['group'] for value in self.hmm_dict.values()))
        self.group_sizes = {group: len(self.readsets) for group in self.group_names}

        self.fasta_information = {}
        for group in self.group_names:
            self.fasta_information[group] = {
                'path': os.path.join(
                    self.dirs_dict['REPRESENTATIVES_DIR'],
                    group,
                    f"{group}-mmseqs_NR_rep_seq.fasta",
                ),
                'external_gene_calls': os.path.join(
                    self.dirs_dict['REPRESENTATIVES_DIR'],
                    group,
                    f"{group}-rep-gene-calls.tsv",
                ),
            }

        self.references_mode = True
        self.remove_short_reads_based_on_references = False
        self.references_for_removal = {}
        self.run_qc = False
        self.set_config_param('all_against_all', True)

        # Pick which tree algorithm
        self.run_iqtree = self.get_param_value_from_config(['iqtree', 'run'])
        self.run_fasttree = self.get_param_value_from_config(['fasttree', 'run'])

        if not self.run_iqtree and not self.run_fasttree:
            raise ConfigError("Please choose either iqtree or fasttree in your config file to run your phylogenetic tree.")

        # Pick clustering method
        self.cluster_representative_method = self.get_param_value_from_config(['cluster_representative_method', 'method'])

        if self.cluster_representative_method not in ['mmseqs', 'cluster_rep_with_coverages']:
            raise ConfigError(f"anvi'o has never heard of this method to pick a cluster representative: {self.cluster_representative_method} "
                              f"Please check your config file {self.config_file} and change cluster_representative_method to one of the following: 'mmseqs' and 'cluster_rep_with_coverages'")

        if self.cluster_representative_method == 'cluster_rep_with_coverages' and len(self.contigs_db_name_bam_dict) == 0:
            raise ConfigError("The EcoPhylo workflow can't use the cluster representative method cluster_rep_with_coverages without BAM files..."
                              "Please edit your metagenomes.txt or external-genomes.txt and add BAM files.")

        if self.cluster_representative_method == 'cluster_rep_with_coverages' and self.AA_mode == True:
            raise ConfigError("The EcoPhylo workflow can't use the cluster representative method cluster_rep_with_coverages in AA_mode")

        # global target files
        self.target_files = self.get_target_files()

    def get_default_config(self):
        """Return default config with output_dirs limited to our 8 canonical keys."""
        c = super().get_default_config()
        c["output_dirs"].pop("QC_DIR", None)
        return c

    def get_target_files(self):
        """This function creates a list of target files for Snakemake

        RETURNS
        =======
        target_files: list
            list of target files for snakemake
        """

        target_files = []

        for hmm, value in self.hmm_dict.items():
            group = value['group']
            hmm_source = value['source']
            target_file = os.path.join(self.dirs_dict['REPRESENTATIVES_DIR'], f"{group}", f"{group}_stats.tsv")
            target_files.append(target_file)

            if not self.samples_txt_file:
                # TREE-MODE
                target_file = os.path.join(self.dirs_dict['PHYLO'], f"{group}", "state_imported_tree.done")
                target_files.append(target_file)

            else:
                # PROFILE-MODE
                target_file = os.path.join(self.dirs_dict['MERGE_DIR'], f"{group}", f"{group}_state_imported_profile.done")
                target_files.append(target_file)

                target_file = os.path.join(self.dirs_dict['MERGE_DIR'], f"{group}", f"{group}_summarize.done")
                target_files.append(target_file)

        return target_files

    def get_target_files_make_anvio_state_file(self):
        """This function creates a list of target files for make_anvio_state_file

        RETURNS
        =======
        target_files: list
            list of target files for snakemake
        """

        target_files = []

        target_file = os.path.join(self.dirs_dict['PROFILE_DIR'], "{group}", "{group}_misc.tsv")
        target_files.append(target_file)

        if self.run_scg_taxonomy and not self.AA_mode:
            target_file = os.path.join(self.dirs_dict['PROFILE_DIR'], "{group}", "anvi_estimate_scg_taxonomy_for_SCGs.done")
            target_files.append(target_file)

        return target_files

    def get_input_files_combine_sequence_data(self, group):
        """This function return a list of input file for the rule combine_sequence_data"""

        input_files = []
        hmm_source_name = []

        # get list of unique hmm sources
        for hmm, value in self.hmm_dict.items():
            if value['group'] == group:
                hmm_source_name.append((value['source'], value['name']))

        # for samples and unique hmm_source, get the input files
        for hmm_source, hmm_name in hmm_source_name:
            input_file = [os.path.join(self.dirs_dict['HMM_HITS_DIR'], sample_name, hmm_source, hmm_name, f"{sample_name}-{hmm_name}-processed.done") for sample_name in self.names_list]
            input_files.extend(input_file)

        return input_files

    def init_hmm_list_txt(self):
        """This function will sanity check hmm-list.txt

        PARAMETERS
        ==========
        self.hmm_list_path : str
            Path to hmm_list.txt

        RETURNS
        =======
        self.hmm_dict : dict
            Dict with hmm (source and name) as primary key and values: hmm_name, hmm_source, PATH, group (optional)
        """
        filesnpaths.is_file_exists(self.hmm_list_path)
        filesnpaths.is_file_tab_delimited(self.hmm_list_path)

        try:
            hmm_df = pd.read_csv(self.hmm_list_path, sep='\t', index_col=False)
        except AttributeError as e:
            raise ConfigError(f"The hmm_list.txt file, {self.hmm_list_path}, does not appear to be properly formatted. "
                              f"This is the error from trying to load it: {self.hmm_list_path}")

        hmm_list_txt_columns = ['name', 'source', 'path']

        for column_name in hmm_list_txt_columns:
            if column_name not in list(hmm_df.columns):
                raise ConfigError(f"Looks like your hmm-list.txt file, {self.hmm_list_path}, is not properly formatted. "
                                  f"We are not sure what's wrong, but we can't find a column with title '{column_name}'."
                                  f"Please make sure you have a tsv with the column names: {hmm_list_txt_columns}")

        if any("-" in s for s in hmm_df['name']):
            raise ConfigError(f"Please do not use '-' in your external hmm names in: "
                              f"{self.hmm_list_path}. It will make our lives "
                              f"easier with Snakemake wildcards :)")

        # create a unique name based on the hmm source and hmm name
        hmm_df['id'] = hmm_df['source'] + '_' + hmm_df['name']

        # the group column is optional and used to combine sequences
        # from multiple HMM source/genes.
        # if no group provided, group name is hmm id (source + name).
        # This "group" will be the main hmm wildcards, similarly to the metagenomics workflow
        if 'group' not in hmm_df:
            hmm_df['group'] = hmm_df['id']

        # to dict
        self.hmm_dict = hmm_df.set_index('id').to_dict('index')

        # FIXME: this line prints the list of hmm_sources to stdout and I don't want that
        self.internal_hmm_sources = list(anvio.data.hmm.sources.keys())

        # make a list of unique hmm source
        self.unique_hmm_source = {}

        # make a list of group id for sanity check
        unique_group = []

        for hmm, value in self.hmm_dict.items():
            hmm_name = value['name']
            hmm_source = value['source']
            hmm_path = value['path']

            if hmm_path == "INTERNAL":
                if hmm_source not in self.internal_hmm_sources:
                    raise ConfigError(f"{hmm_source} is not an 'INTERNAL' hmm source for anvi'o. "
                                      f"Please double check {self.hmm_list_path} to see if you spelled it right or "
                                      f"please checkout the default internal hmms here: https://merenlab.org/software/anvio/help/7/artifacts/hmm-source/#default-hmm-sources")
                if hmm_name not in constants.default_scgs_for_taxonomy and self.run_scg_taxonomy:
                    raise ConfigError(f"You asked EcoPhylo to use anvi-estimate-scg-taxonomy but the HMM {hmm_name} in {hmm_source} is not compatible. "
                                      f"You can either turn off anvi-estimate-scg-taxonomy in the config file, or choose a compatible gene in this set: "
                                      f"{constants.default_scgs_for_taxonomy}")

            if not filesnpaths.is_file_exists(hmm_path, dont_raise=True):
                if hmm_path == 'INTERNAL':
                    pass
                else:
                    raise ConfigError(f"The path to your hmm {hmm_name} does not exist: {hmm_path}. "
                                      f"Please double check the paths in our hmm-list.txt: {self.hmm_list_path} "
                                      f"If the hmm you want to use is in an internal anvi'o hmm collection e.g. Bacteria_71 "
                                      f"please put 'INTERNAL' for the path.")

            if hmm_path != "INTERNAL":
                sources = u.get_HMM_sources_dictionary([hmm_path])

                for source,value in sources.items():
                    gene = value['genes']
                    if hmm_source != source:
                        raise ConfigError(f"In your {self.hmm_list_path}, please change the source for gene {hmm_name} to this: {source}")
                    if len(gene) > 1:
                        raise ConfigError("EcoPhylo can only work with one gene at a time in a hmm directory (at the moment)")
                    if hmm_name != gene[0]:
                        raise ConfigError(f"In your {self.hmm_list_path}, please change the gene name {hmm_name} to this: {gene[0]}")

            if hmm_source not in self.unique_hmm_source:
                self.unique_hmm_source[hmm_source] = hmm_path

            if value['group'] not in unique_group:
                unique_group.append(value['group'])

        # Since we now build a reps contigs DB per group (see clustering_alignment.smk), SCG taxonomy
        # works even when groups combine two or more HMMs. The restriction is lifted.
