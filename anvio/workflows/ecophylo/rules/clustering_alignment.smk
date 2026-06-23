rule cluster_X_percent_sim_mmseqs:
    """Cluster NT or AA fasta file with user defined percent identity"""
    input:
        done=rules.combine_sequence_data.output.done,
    output:
        fasta=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"],
            "{group}",
            "{group}-mmseqs_NR_rep_seq.fasta",
        ),
        mmseqs_cluster_rep_index=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"],
            "{group}",
            "{group}-mmseqs_NR_cluster.tsv",
        ),
        done=touch(
            os.path.join(
                dirs_dict["REPRESENTATIVES_DIR"],
                "{group}",
                "{group}-mmseqs_NR_cluster.done",
            )
        ),
    log:
        rule_log("cluster_X_percent_sim_mmseqs", "cluster_X_percent_sim_mmseqs_{group}"),
    threads: M.T("cluster_X_percent_sim_mmseqs")
    params:
        output_prefix=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"], "{group}", "{group}-mmseqs_NR"
        ),
        mmseqs_tmp=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"], "{group}", "{group}-tmp"
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
                dirs_dict["HMM_HITS_DIR"],
                f"{wildcards.group}",
                f"{wildcards.group}-all.faa",
            )
        else:
            fasta = os.path.join(
                dirs_dict["HMM_HITS_DIR"],
                f"{wildcards.group}",
                f"{wildcards.group}-all.fna",
            )
        # Exit workflow if we couldn't find any hmm-hits
        fasta_dict = SeqIO.index(fasta, "fasta")
        if len(fasta_dict) == 0:
            raise ConfigError(
                f"anvi'o and the EcoPhylo workflow are sad to announce that the "
                f"hmm, {wildcards.group}, was not found in any of your contigs_dbs"
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




if M.cluster_representative_method == "cluster_rep_with_coverages":

    rule anvi_profile_blitz:
        """Choose a NT cluster rep based on read recruitment!
        The sequence with the most read recruitment from the input profiled assembly will be chosen as the cluster representative.
        """
        output:
            target=os.path.join(
                dirs_dict["HMM_HITS_DIR"],
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
                    dirs_dict["HMM_HITS_DIR"],
                    "{sample_name}-gene-coverages.txt",
                ),
                sample_name=M.names_list,
            ),
        output:
            txt=os.path.join(
                dirs_dict["HMM_HITS_DIR"], "gene-coverages.txt"
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
            egc_all=rules.combine_sequence_data.output.external_gene_calls_all,
            coverages=os.path.join(
                dirs_dict["HMM_HITS_DIR"], "gene-coverages.txt"
            ),
        output:
            coverage_reps=os.path.join(
                dirs_dict["REPRESENTATIVES_DIR"],
                "{group}",
                "{group}-coverage-headers.txt",
            ),
            coverage_cluster_rep_index=os.path.join(
                dirs_dict["REPRESENTATIVES_DIR"],
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
            egc_all = pd.read_csv(input.egc_all, sep="\t", index_col=False)
            bam = pd.read_csv(input.coverages, sep="\t", index_col=False)
            bam["primary_key"] = (
                bam["contig"] + "_" + bam["gene_callers_id"].astype(str)
            )
            egc_all["primary_key"] = (
                egc_all["contig"] + "_" + egc_all["gene_callers_id"].astype(str)
            )
            df = pd.merge(egc_all, bam, on="primary_key", how="inner")
            df2 = pd.merge(
                cluster_rep_index,
                df,
                left_on="cluster_members",
                right_on="header",
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
                dirs_dict["REPRESENTATIVES_DIR"],
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
                dirs_dict["REPRESENTATIVES_DIR"],
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
                dirs_dict["REPRESENTATIVES_DIR"], "{group}", "{group}-headers.tmp"
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
            dirs_dict["REPRESENTATIVES_DIR"], "{group}", "{group}-AA_subset.fa"
        ),
    output:
        fasta=os.path.join(dirs_dict["PHYLO"], "{group}", "{group}-aligned.fa"),
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
        fasta=os.path.join(dirs_dict["PHYLO"], "{group}", "{group}_aligned_trimmed.fa"),
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
            dirs_dict["PHYLO"], "{group}", "{group}_aligned_trimmed_filtered.fa"
        ),
        seq_counts_tsv=os.path.join(dirs_dict["PHYLO"], "{group}", "{group}_gaps_counts.tsv"),
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
        gap_counts_prefix=lambda wildcards: os.path.join(
            dirs_dict["PHYLO"], wildcards.group, f"{wildcards.group}_gaps_counts"
        ),
    shell:
        "anvi-script-reformat-fasta {input} -o {output.fasta} \
                                      --max-percentage-gaps {params.max_percentage_gaps} \
                                       --export-gap-counts-table {params.gap_counts_prefix} >> {log} 2>&1"


rule extract_QCd_sequence_headers:
    """Extract headers from QCd sequences for misc data generation"""
    input:
        fasta=rules.remove_sequences_with_X_percent_gaps.output.fasta,
    output:
        headers=os.path.join(dirs_dict["PHYLO"], "{group}", "{group}_headers.tmp"),
    threads: M.T("extract_QCd_sequence_headers")
    shell:
        "grep '^>' {input.fasta} | sed 's/>//g' > {output.headers}"


rule count_num_sequences_filtered:
    """Record the number of sequences filtered at each step of the workflow"""
    input:
        remove_seq_with_gaps=rules.remove_sequences_with_X_percent_gaps.output.fasta,
    output:
        target=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"], "{group}", "{group}_stats.tsv"
        ),
    log:
        rule_log("count_num_sequences_filtered", "count_num_sequences_filtered_{group}"),
    threads: M.T("count_num_sequences_filtered")
    params:
        combined_seq=rules.combine_sequence_data.output.NT_all,
        cluster_mmseqs=rules.cluster_X_percent_sim_mmseqs.output.fasta,
    run:
        def count_num_sequences(fasta_path):
            return len(SeqIO.index(fasta_path, "fasta"))
        input_files_list = [
            params.combined_seq,
            params.cluster_mmseqs,
            input.remove_seq_with_gaps,
        ]
        num_seqs_list = [
            count_num_sequences(fasta) for fasta in input_files_list
        ]
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
            lines = [col_names, step1, step2, step3]
            for line in lines:
                f.write("\t".join(line) + "\n")


if not M.AA_mode:

    rule build_rep_external_gene_calls:
        """Build external gene calls TSV by subsetting combined EGC to rep sequences"""
        input:
            reps=rules.cluster_X_percent_sim_mmseqs.output.fasta,
            egc_all=rules.combine_sequence_data.output.external_gene_calls_all,
        output:
            gene_calls=os.path.join(
                dirs_dict["REPRESENTATIVES_DIR"],
                "{group}",
                "{group}-rep-gene-calls.tsv",
            ),
        log:
            rule_log("build_rep_external_gene_calls", "build_rep_external_gene_calls_{group}"),
        threads: M.T("build_rep_external_gene_calls")
        run:
            rep_headers = set(rec.id for rec in SeqIO.parse(input.reps, "fasta"))

            egc = pd.read_csv(input.egc_all, sep="\t")
            egc = egc[egc['header'].isin(rep_headers)].copy()

            if egc.empty:
                egc = pd.DataFrame(columns=[
                    "gene_callers_id", "contig", "start", "stop",
                    "direction", "partial", "call_type", "source",
                    "version", "aa_sequence",
                ])
            else:
                egc['contig'] = egc['header']
                egc['gene_callers_id'] = range(1, len(egc) + 1)
                egc = egc[["gene_callers_id", "contig", "start", "stop",
                           "direction", "partial", "call_type", "source",
                           "version", "aa_sequence"]]

            egc.to_csv(output.gene_calls, sep="\t", index=False)

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
        shell:
            "anvi-gen-contigs-database -f {input.fasta} -o {output.db} --external-gene-calls {input.gene_calls} -T {threads} >> {log} 2>&1"

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
        run:
            if M.run_scg_taxonomy:
                shell("anvi-run-hmms -c {input.db} --num-threads {threads} >> {log} 2>&1 && \
                       anvi-run-scg-taxonomy -c {input.db} --num-threads {threads} {params.additional_params} >> {log} 2>&1")

    rule anvi_estimate_scg_taxonomy_reps:
        """Export SCG taxonomy for each representative from the reps contigs DB"""
        input:
            db=rules.anvi_run_scg_taxonomy_reps.output.done,
        output:
            done=touch(os.path.join(dirs_dict["PROFILE_DIR"], "{group}", "anvi_estimate_scg_taxonomy_for_SCGs.done")),
            tax_data_final=os.path.join(dirs_dict["PROFILE_DIR"], "{group}", "{group}_scg_taxonomy_data.tsv"),
        log:
            rule_log("anvi_estimate_scg_taxonomy_reps", "anvi_estimate_scg_taxonomy_reps_{group}"),
        threads: M.T("anvi_estimate_scg_taxonomy")
        run:
            if M.run_scg_taxonomy:
                reps_db = os.path.join(dirs_dict["CONTIGS_DIR"], f"{wildcards.group}.db")

                anvio_db = db.DB(reps_db, None, ignore_version=True)

                scg_taxonomy_dict = anvio_db.get_table_as_dict('scg_taxonomy')
                scg_taxonomy = pd.DataFrame.from_dict(scg_taxonomy_dict, orient='index')

                genes_dict = anvio_db.get_table_as_dict('genes_in_contigs')
                gene_to_contig = {gcid: row['contig'] for gcid, row in genes_dict.items()}
                scg_taxonomy['split_name'] = scg_taxonomy.index.map(
                    lambda gcid: f"{gene_to_contig[gcid]}_split_00001"
                )
                scg_taxonomy['identifier'] = scg_taxonomy.index.astype(str)

                expected_columns = ["split_name", "identifier", "percent_identity",
                                    "t_domain", "t_phylum", "t_class", "t_order",
                                    "t_family", "t_genus", "t_species"]
                available_columns = [c for c in expected_columns if c in scg_taxonomy.columns]
                scg_taxonomy = scg_taxonomy[available_columns]
                scg_taxonomy.to_csv(output.tax_data_final, sep="\t", index=None, na_rep="NA")
            else:
                # Write empty TSV with headers so downstream os.path.isfile() checks pass
                expected_columns = ["split_name", "identifier", "percent_identity",
                                    "t_domain", "t_phylum", "t_class", "t_order",
                                    "t_family", "t_genus", "t_species"]
                pd.DataFrame(columns=expected_columns).to_csv(output.tax_data_final, sep="\t", index=None)
