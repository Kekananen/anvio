def extract_misc_data(mmseqs_cluster_rep_index, final_sequences_headers, output):
    """This function creates a tsv of metadata for mmseqs cluster tsv outfile (*_cluster.tsv)

    Parameters
    ==========
    mmseqs_cluster_rep_index: str
        path to tsv containing mmseqs cluster representatives and member index

    final_sequences_headers: str
        path to tsv with fasta file headers

    output:
        target=str
        path to output tsv, metadata values include ['split_name', 'contigs_db_type', 'genomic_seq_in_cluster', 'cluster_size']
        - 'split_name': ID name in interactive interface
        - 'contigs_db_type': cluster representative origin: MAG, SAG, isolate genome, metagenome.
        - 'genomic_seq_in_cluster': detects if a cluster member came from an external-genome
        - 'cluster_size': number of sequences in cluster
    """

    # Import data
    # ------------
    cluster_rep_index = pd.read_csv(
        mmseqs_cluster_rep_index,
        sep="\t",
        index_col=False,
        names=["representative", "cluster_members"],
    )

    final_sequences_headers = pd.read_csv(
        final_sequences_headers, sep="\t", index_col=False, names=["identifier"]
    )

    # Clean data
    # ------------

    # Detect if there is a genomic reference protein in cluster
    cluster_rep_index_dict = (
        cluster_rep_index.groupby("representative")["cluster_members"]
        .apply(list)
        .to_dict()
    )

    cluster_reps_with_genomic_references_list = []
    for seq in final_sequences_headers.iloc[:, 0].tolist():
        cluster_members_list = cluster_rep_index_dict[seq]
        for external_genome in M.external_genomes_names_list:
            check = any(external_genome in s for s in cluster_members_list)
            if check is True:
                cluster_reps_with_genomic_references_list.append(seq)

    # Count size of clusters
    df = cluster_rep_index.groupby("representative").apply(count_cluster_size)[
        ["cluster_members", "cluster_size"]
    ]

    # Make split names for anvi-interactive
    df["split_name"] = df["cluster_members"].astype(str) + "_split_00001"

    # subset misc data to final set of proteins
    df = pd.merge(
        df,
        final_sequences_headers,
        left_on="cluster_members",
        right_on="identifier",
        how="inner",
    )

    df["genomic_seq_in_cluster"] = np.where(
        df["cluster_members"].isin(cluster_reps_with_genomic_references_list),
        "yes",
        "no",
    )

    # Determine contigs_db type: metagenome or genomes
    # FIXME: This will need to be changed in the future to accomidate SAGs, MAGs, and other genomic sources
    # If metagenome_name is in external_genomes_names_list then it's a genome
    contigs_db_type_dict = {}
    for name in list(df.cluster_members):
        if any(x in name for x in M.external_genomes_names_list):
            contigs_db_type_dict[name] = "genome"
        else:
            contigs_db_type_dict[name] = "metagenome"

    df["contigs_db_type"] = df.cluster_members.apply(
        lambda name: contigs_db_type_dict[name]
    )

    # grab the final columns
    df = df[["split_name", "contigs_db_type", "genomic_seq_in_cluster", "cluster_size"]]

    # Export
    # -------
    df.to_csv(output, sep="\t", index=None, na_rep="NA")


def count_cluster_size(group):
    """Count the number of sequences represented by each cluster."""
    c = group["cluster_members"].count()
    group["cluster_size"] = c

    return group


rule make_misc_data:
    """Make misc data file for clustered sequences"""
    input:
        final_list_of_sequences_for_mapping_headers=rules.extract_QCd_sequence_headers.output.headers,
    output:
        misc_data_final=os.path.join(
            dirs_dict["MISC_DATA"], "{group}", "{group}_misc.tsv"
        ),
    log:
        rule_log("make_misc_data", "add_contigs_db_type_{group}"),
    threads: M.T("add_misc_data_to_taxonomy")
    params:
        mmseqs_cluster_rep_index=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"],
            "{group}",
            "{group}-mmseqs_NR_cluster.tsv",
        ),
        coverage_cluster_rep_index=os.path.join(
            dirs_dict["REPRESENTATIVES_DIR"],
            "{group}",
            "{group}-coverage_cluster.tsv",
        ),
    run:
        """Here we determine the origin of each SCG (which kind of contigs_db): metagenome, isolate genome, etc."""
        if M.cluster_representative_method == "mmseqs":
            extract_misc_data(
                mmseqs_cluster_rep_index=params.mmseqs_cluster_rep_index,
                final_sequences_headers=input.final_list_of_sequences_for_mapping_headers,
                output=output.misc_data_final,
            )
        if M.cluster_representative_method == "cluster_rep_with_coverages":
            extract_misc_data(
                mmseqs_cluster_rep_index=params.coverage_cluster_rep_index,
                final_sequences_headers=input.final_list_of_sequences_for_mapping_headers,
                output=output.misc_data_final,
            )




