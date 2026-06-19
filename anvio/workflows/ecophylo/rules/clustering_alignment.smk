rule cluster_X_percent_sim_mmseqs:
    """Cluster NT or AA fasta file with user defined percent identity"""
    input:
        done=rules.combine_sequence_data.output.done,
    output:
        fasta=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-mmseqs_NR_rep_seq.fasta",
        ),
        mmseqs_cluster_rep_index=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-mmseqs_NR_cluster.tsv",
        ),
        done=touch(
            os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-mmseqs_NR_cluster.done",
            )
        ),
    log:
        rule_log("cluster_X_percent_sim_mmseqs", "cluster_X_mmseqs_{group}"),
    threads: M.T("cluster_X_percent_sim_mmseqs")
    params:
        output_prefix=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-mmseqs_NR"
        ),
        mmseqs_tmp=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-tmp"
        ),
        min_seq_id=M.get_param_value_from_config(
            ["cluster_X_percent_sim_mmseqs", "--min-seq-id"]
        ),
        cov_mode=M.get_param_value_from_config(
            ["cluster_X_percent_sim_mmseqs", "--cov-mode"]
        ),
        additional_params=M.get_param_value_from_config(
            ["cluster_X_percent_sim_mmseqs", "additional_params"]
        ),
    run:
        if M.AA_mode == True:
            fasta = os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                f"{wildcards.group}",
                f"{wildcards.group}-all.faa",
            )
        else:
            fasta = os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                f"{wildcards.group}",
                f"{wildcards.group}-all.fna",
            )
        # Exit workflow if we couldn't find any hmm-hits
        from Bio import SeqIO

        fasta_dict = SeqIO.index(fasta, "fasta")
        if len(fasta_dict) == 0:
            raise ConfigError(
                f"anvi'o and the EcoPhylo workflow are sad to announce that the "
                f"hmm, {wildcards.group}, was not not found in any of your contigs_dbs"
            )
        shell(
            f"mmseqs easy-cluster {fasta} \
                                    {params.output_prefix} \
                                    {params.mmseqs_tmp} \
                                    --threads {threads} \
                                    --min-seq-id {params.min_seq_id} \
                                    --cov-mode {params.cov_mode} \
                                    {params.additional_params} >> {log} 2>&1"
        )


rule cluster_X_percent_sim_mmseqs_OTUs:
    """Cluster extracted proteins within percent identity parameter space provided by user.
This will help identify clustering thresholds for OTU like analyses.
"""
    input:
        done=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{group}-mmseqs_NR_cluster.done",
        ),
    output:
        fasta=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{clustering_threshold}",
            "{group}-{clustering_threshold}-mmseqs_NR_rep_seq.fasta",
        ),
        mmseqs_cluster_rep_index=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{clustering_threshold}",
            "{group}-{clustering_threshold}-mmseqs_NR_cluster.tsv",
        ),
    log:
        rule_log(
            "cluster_X_percent_sim_mmseqs_OTUs",
            "cluster_X_mmseqs_{group}_{clustering_threshold}",
        ),
    threads: M.T("cluster_X_percent_sim_mmseqs")
    params:
        NT_all=rules.combine_sequence_data.output.NT_all,
        min_seq_id=lambda wildcards: M.clustering_threshold_dict[
            wildcards.clustering_threshold
        ],
        output_prefix=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{clustering_threshold}",
            "{group}-{clustering_threshold}-mmseqs_NR",
        ),
        mmseqs_tmp=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
            "{group}",
            "{clustering_threshold}",
            "{group}-{clustering_threshold}-tmp",
        ),
        cov_mode=M.get_param_value_from_config(
            ["cluster_X_percent_sim_mmseqs", "--cov-mode"]
        ),
        additional_params=M.get_param_value_from_config(
            ["cluster_X_percent_sim_mmseqs", "additional_params"]
        ),
    shell:
        "mmseqs easy-cluster {params.NT_all} \
                               {params.output_prefix} \
                               {params.mmseqs_tmp} \
                               --threads {threads} \
                               --min-seq-id {params.min_seq_id} \
                               --cov-mode {params.cov_mode} \
                               {params.additional_params} >> {log} 2>&1"


if M.cluster_representative_method == "cluster_rep_with_coverages":

    rule anvi_profile_blitz:
        """Choose a NT cluster rep based on read recruitment!
        The sequence with the most read recruitment from the input profiled assembly will be chosen as the cluster representative.
        """
        output:
            target=os.path.join(
                dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                "{sample_name}-gene-coverages.txt",
            ),
        log:
            rule_log("anvi_profile_blitz", "anvi_profile_blitz-{sample_name}"),
        threads: M.T("anvi_profile_blitz")
        params:
            contigs_db=lambda wildcards: os.path.join(
                M.contigs_db_name_path_dict[wildcards.sample_name]
            ),
            bam=lambda wildcards: os.path.join(
                M.contigs_db_name_bam_dict[wildcards.sample_name]
            ),
        shell:
            "anvi-profile-blitz {params.bam} -c {params.contigs_db} --gene-mode --report-minimal -o {output} >> {log} 2>&1"

    rule cat_anvi_profile_blitz:
        """Cat gene coverages from anvi-profile-blitz"""

        # log: os.path.join(dirs_dict['LOGS_DIR'], "anvi_profile_blitz-{sample_name}.log")
        input:
            targets=expand(
                os.path.join(
                    dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"],
                    "{sample_name}-gene-coverages.txt",
                ),
                sample_name=M.names_list,
            ),
        output:
            txt=os.path.join(
                dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"], "gene-coverages.txt"
            ),
        threads: M.T("anvi_profile_blitz")
        shell:
            """
            echo -e 'gene_callers_id\tcontig\tsample\tlength\tdetection\tmean_cov' > {output}
            awk 'FNR>1' {input} >> {output}
            """

    rule pick_cluster_rep_with_coverage:
        """Pick a cluster rep with coverage values."""
        input:
            mmseqs_cluster_rep_index=rules.cluster_X_percent_sim_mmseqs.output.mmseqs_cluster_rep_index,
            reformat_report=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-reformat-report-all.txt",
            ),
            coverages=os.path.join(
                dirs_dict["EXTRACTED_RIBO_PROTEINS_DIR"], "gene-coverages.txt"
            ),
        output:
            coverage_reps=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-coverage-headers.txt",
            ),
            coverage_cluster_rep_index=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-coverage_cluster.tsv",
            ),
        log:
            rule_log(
                "pick_cluster_rep_with_coverage",
                "pick_cluster_rep_with_coverage-{group}",
            ),
        run:
            # Bind anvi-profile-blitz data with cluster rep data and group_by cluster rep then find the cluster member with the highest coverage to pick new rep
            pd.set_option("expand_frame_repr", False)
            cluster_rep_index = pd.read_csv(
                input.mmseqs_cluster_rep_index,
                sep="\t",
                index_col=False,
                names=["representative", "cluster_members"],
            )
            reformat_report = pd.read_csv(
                input.reformat_report,
                sep="\t",
                index_col=False,
                names=["new_header", "header"],
            )
            bam = pd.read_csv(input.coverages, sep="\t", index_col=False)
            bam["primary_key"] = (
                bam["contig"] + "_" + bam["gene_callers_id"].astype(str)
            )
            reformat_report["gene_callers_id"] = (
                reformat_report["header"]
                .str.split("gene_callers_id:|\|start:", expand=True)[1]
                .astype(str)
            )
            reformat_report["contig"] = (
                reformat_report["header"]
                .str.split("contig:|\|gene_callers_id:", expand=True)[1]
                .astype(str)
            )
            reformat_report["primary_key"] = (
                reformat_report["contig"]
                + "_"
                + reformat_report["gene_callers_id"].astype(str)
            )
            df = pd.merge(reformat_report, bam, on="primary_key", how="inner")
            df2 = pd.merge(
                cluster_rep_index,
                df,
                left_on="cluster_members",
                right_on="new_header",
                how="inner",
            )[["representative", "cluster_members", "mean_cov"]]

            def get_new_seed(group):
                """Find the cluster member with the highest coverage and set it as the new representative."""
                best_member = group.loc[group["mean_cov"].idxmax(), "cluster_members"]
                group["representative"] = best_member
                return group

            df3 = df2.groupby("representative").apply(get_new_seed)
            # export headers of representatives and cluster rep index
            df3[["representative"]].drop_duplicates().to_csv(
                output.coverage_reps, sep="\t", index=None, header=False
            )
            df3[["representative", "cluster_members"]].to_csv(
                output.coverage_cluster_rep_index, sep="\t", index=None, na_rep="NA"
            )

    rule subset_AA_seqs_with_coverage_reps:
        """Subset AA sequences for the mmseqs cluster representatives"""
        input:
            fa=rules.combine_sequence_data.output.AA_all,
            coverage_reps=rules.pick_cluster_rep_with_coverage.output.coverage_reps,
        output:
            fasta=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-AA_subset.fa",
            ),
        log:
            rule_log(
                "subset_AA_seqs_with_coverage_reps",
                "subset_AA_seqs_with_coverage_reps_{group}",
            ),
        threads: M.T("subset_AA_seqs_with_coverage_reps")
        shell:
            "anvi-script-reformat-fasta {input.fa} -I {input.coverage_reps} -o {output.fasta} >> {log} 2>&1"


if M.cluster_representative_method == "mmseqs":

    rule subset_AA_seqs_with_mmseqs_reps:
        """Subset AA sequences for the mmseqs cluster representatives"""
        input:
            mmseqs_reps=rules.cluster_X_percent_sim_mmseqs.output.fasta,
        output:
            fasta=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-AA_subset.fa",
            ),
        log:
            rule_log(
                "subset_AA_seqs_with_mmseqs_reps",
                "subset_AA_seqs_with_mmseqs_reps_{group}",
            ),
        threads: M.T("subset_AA_seqs_with_mmseqs_reps")
        params:
            fa=rules.combine_sequence_data.output.AA_all,
            headers=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-headers.tmp"
            ),
        shell:
            """
            grep '>' {input.mmseqs_reps} | sed 's/>//g' > {params.headers}
            anvi-script-reformat-fasta {params.fa} -I {params.headers} -o {output.fasta} >> {log} 2>&1
            """


rule align_sequences:
    """MSA of AA sequences subset."""
    input:
        source=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"], "{group}", "{group}-AA_subset.fa"
        ),
    output:
        fasta=os.path.join(dirs_dict["MSA"], "{group}", "{group}-aligned.fa"),
    log:
        rule_log("align_sequences", "align_sequences_{group}"),
    threads: M.T("align_sequences")
    params:
        additional_params=M.get_param_value_from_config(
            ["align_sequences", "additional_params"]
        ),
    shell:
        "muscle -in {input} -out {output} {params.additional_params} -verbose 2> {log}"


rule trim_alignment:
    """Trim MSA alignment"""
    input:
        fasta=rules.align_sequences.output.fasta,
    output:
        fasta=os.path.join(dirs_dict["MSA"], "{group}", "{group}_aligned_trimmed.fa"),
    log:
        rule_log("trim_alignment", "trim_alignment_{group}"),
    threads: M.T("trim_alignment")
    params:
        gt=M.get_param_value_from_config(["trim_alignment", "-gt"]),
        gappyout=M.get_rule_param("trim_alignment", "-gappyout"),
        additional_params=M.get_param_value_from_config(
            ["trim_alignment", "additional_params"]
        ),
    shell:
        "trimal -in {input} -out {output} {params.gappyout} {params.additional_params} 2> {log}"


rule remove_sequences_with_X_percent_gaps:
    """Remove sequences that have X% gaps"""
    input:
        fasta=rules.trim_alignment.output.fasta,
    output:
        fasta=os.path.join(
            dirs_dict["MSA"], "{group}", "{group}_aligned_trimmed_filtered.fa"
        ),
        seq_counts_tsv=os.path.join(dirs_dict["MSA"], "{group}", "{group}_gaps_counts"),
    log:
        rule_log(
            "remove_sequences_with_X_percent_gaps",
            "remove_sequences_with_X_percent_gaps_{group}",
        ),
    threads: M.T("remove_sequences_with_X_percent_gaps")
    params:
        max_percentage_gaps=M.get_param_value_from_config(
            ["remove_sequences_with_X_percent_gaps", "--max-percentage-gaps"]
        ),
    shell:
        "anvi-script-reformat-fasta {input} -o {output.fasta} \
                                      --max-percentage-gaps {params.max_percentage_gaps} \
                                      --export-gap-counts-table {output.seq_counts_tsv} >> {log} 2>&1"


rule extract_QCd_sequence_headers:
    """Extract headers from QCd sequences for misc data generation"""
    input:
        fasta=rules.remove_sequences_with_X_percent_gaps.output.fasta,
    output:
        headers=os.path.join(dirs_dict["MSA"], "{group}", "{group}_headers.tmp"),
    threads: M.T("extract_QCd_sequence_headers")
    shell:
        "grep '^>' {input.fasta} | sed 's/>//g' > {output.headers}"


rule count_num_sequences_filtered:
    """Record the number of sequences filtered at each step of the workflow"""
    input:
        remove_seq_with_gaps=rules.remove_sequences_with_X_percent_gaps.output.fasta,
        clustering_thresholds=expand(
            os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{{group}}",
                "{clustering_threshold}",
                "{{group}}-{clustering_threshold}-mmseqs_NR_rep_seq.fasta",
            ),
            clustering_threshold=M.clustering_param_space_list_strings,
        ),
    output:
        target=os.path.join(
            dirs_dict["RIBOSOMAL_PROTEIN_MSA_STATS"], "{group}", "{group}_stats.tsv"
        ),
    log:
        rule_log("count_num_sequences_filtered", "count_num_sequences_filtered_{group}"),
    threads: M.T("count_num_sequences_filtered")
    params:
        combined_seq=rules.combine_sequence_data.output.NT_all,
        cluster_mmseqs=rules.cluster_X_percent_sim_mmseqs.output.fasta,
    run:
        def count_num_sequences(fasta):
            """This function counts the number of sequences in a fasta file

            Parameters
            ==========
            fasta: fasta

            Returns
            =======
            num_seqs : int
            """

            num_seqs = 0

        for line in fasta:
            if line.startswith(">"):
                num_seqs += 1
        return num_seqs
        input_files_list = [
            params.combined_seq,
            params.cluster_mmseqs,
            input.remove_seq_with_gaps,
        ]
        num_seqs_list = [
            count_num_sequences(open(fasta)) for fasta in input_files_list
        ]
        clustering_threshold_attributes_list = []
        for file in input.clustering_thresholds:
            path = file
            threshold = file.split("/")[3]
            with open(file) as fasta:
                num_seqs = count_num_sequences(fasta)
            clustering_threshold_attributes = [str(threshold), str(num_seqs), path]
            clustering_threshold_attributes_list.append(
                clustering_threshold_attributes
            )
        with open(output.target, "w") as f:
            col_names = ["rule_name", "num_sequences_left", "rel_path"]
            step1 = [
                "combine_sequence_data",
                str(num_seqs_list[0]),
                params.combined_seq,
            ]
            step2 = [
                "cluster_X_percent_sim_mmseqs",
                str(num_seqs_list[1]),
                params.cluster_mmseqs,
            ]
            step3 = [
                "remove_sequences_with_X_percent_gaps",
                str(num_seqs_list[2]),
                input.remove_seq_with_gaps,
            ]
            lines = [
                col_names,
                step1,
                step2,
                step3,
            ] + clustering_threshold_attributes_list
            for line in lines:
                f.write("\t".join(line) + "\n")


if not M.AA_mode:

    rule build_rep_external_gene_calls:
        """Build external gene calls TSV from NT and AA representative sequences"""
        input:
            nt_reps=rules.cluster_X_percent_sim_mmseqs.output.fasta,
            aa_reps=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-AA_subset.fa",
            ),
        output:
            gene_calls=os.path.join(
                dirs_dict["RIBOSOMAL_PROTEIN_FASTAS"],
                "{group}",
                "{group}-rep-gene-calls.tsv",
            ),
        log:
            rule_log("build_rep_external_gene_calls", "build_rep_external_gene_calls_{group}"),
        threads: M.T("build_rep_external_gene_calls")
        run:
            from Bio import SeqIO

            nt_index = SeqIO.index(input.nt_reps, "fasta")
            aa_index = SeqIO.index(input.aa_reps, "fasta")

            records = []
            for i, header in enumerate(nt_index, start=1):
                nt_seq = nt_index[header]
                aa_seq = str(aa_index[header].seq)
                records.append({
                    'gene_callers_id': i,
                    'contig': header,
                    'start': 1,
                    'stop': len(nt_seq.seq),
                    'direction': 'f',
                    'partial': 0,
                    'call_type': 1,
                    'source': 'EcoPhylo',
                    'version': anvio.__version__,
                    'aa_sequence': aa_seq,
                })

            nt_index.close()
            aa_index.close()

            pd.DataFrame(records).to_csv(output.gene_calls, sep="\t", index=False)

    # NOTE: The reps contigs DB is built BEFORE the alignment/trimming/gap-filtering
    # steps. This means the DB contains ALL cluster representatives, even those that
    # may later be removed due to excessive gaps in the MSA. The interactive interface
    # will display all representatives from the contigs DB, not just those that
    # survived the gap filter. If this causes issues (e.g., tree tips with no
    # corresponding entries), check how the interactive interface handles the
    # mismatch and consider either (a) building the DB post-QC, or (b) filtering
    # the DB entries to match the final tree.

    rule anvi_gen_contigs_database_reps:
        """Generate contigs database from representative NT sequences"""
        input:
            fasta=rules.cluster_X_percent_sim_mmseqs.output.fasta,
            gene_calls=rules.build_rep_external_gene_calls.output.gene_calls,
        output:
            db=os.path.join(dirs_dict["CONTIGS_DIR"], "{group}.db"),
        log:
            rule_log("anvi_gen_contigs_database_reps", "anvi_gen_contigs_database_reps_{group}"),
        threads: M.T("anvi_gen_contigs_database_reps")
        params:
            skip_gene_calling="--skip-gene-calling",
        shell:
            "anvi-gen-contigs-database -f {input.fasta} -o {output.db} {params.skip_gene_calling} --external-gene-calls {input.gene_calls} -T {threads} >> {log} 2>&1"

    rule anvi_run_scg_taxonomy_reps:
        """Run SCG taxonomy on the representative sequences contigs database"""
        input:
            db=rules.anvi_gen_contigs_database_reps.output.db,
        output:
            done=touch(os.path.join(dirs_dict["CONTIGS_DIR"], "{group}", "scg_taxonomy.done")),
        log:
            rule_log("anvi_run_scg_taxonomy_reps", "anvi_run_scg_taxonomy_reps_{group}"),
        threads: M.T("anvi_run_scg_taxonomy")
        params:
            additional_params=M.get_param_value_from_config(["anvi_run_scg_taxonomy", "additional_params"]),
        shell:
            "anvi-run-scg-taxonomy -c {input.db} --num-threads {threads} {params.additional_params} >> {log} 2>&1"

    rule anvi_estimate_scg_taxonomy_reps:
        """Export SCG taxonomy for each representative from the reps contigs DB"""
        input:
            db=rules.anvi_run_scg_taxonomy_reps.output.done,
        output:
            done=touch(os.path.join(dirs_dict["MISC_DATA"], "{group}", "anvi_estimate_scg_taxonomy_for_SCGs.done")),
            tax_data_final=os.path.join(dirs_dict["MISC_DATA"], "{group}", "{group}_scg_taxonomy_data.tsv"),
        log:
            rule_log("anvi_estimate_scg_taxonomy_reps", "anvi_scg_taxonomy_reps_{group}"),
        threads: M.T("anvi_estimate_scg_taxonomy")
        params:
            per_scg=os.path.join(dirs_dict["MISC_DATA"], "{group}", "{group}_per_scg.txt"),
        run:
            # Get the HMM name for this group
            hmm_name = ""
            for hmm, value in M.hmm_dict.items():
                if value["group"] == wildcards.group:
                    hmm_name = value["name"]

            reps_db = os.path.join(dirs_dict["CONTIGS_DIR"], f"{wildcards.group}.db")

            shell("anvi-estimate-scg-taxonomy -c {reps_db} \
                                               --metagenome-mode \
                                               --scg-name-for-metagenome-mode {hmm_name} \
                                               --per-scg-output-file {params.per_scg} \
                                               -o /dev/null \
                                               -T {threads} >> {log} 2>&1")

            scg_taxonomy = pd.read_csv(params.per_scg, sep="\t", index_col=False)

            scg_taxonomy["split_name"] = scg_taxonomy["bin_name"].astype(str) + "_split_00001"
            scg_taxonomy = scg_taxonomy.rename(columns={"bin_name": "identifier"})

            expected_columns = ["split_name", "identifier", "percent_identity",
                                "t_domain", "t_phylum", "t_class", "t_order",
                                "t_family", "t_genus", "t_species"]
            available_columns = [c for c in expected_columns if c in scg_taxonomy.columns]
            scg_taxonomy = scg_taxonomy[available_columns]
            scg_taxonomy.to_csv(output.tax_data_final, sep="\t", index=None, na_rep="NA")
