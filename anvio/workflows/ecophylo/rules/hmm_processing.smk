import shutil

def get_hmm_threads(wildcards):
    """Select threads for anvi-run-hmms based on if a contigs-db is a metagenome or not."""
    threads = M.T("extract_hmm_hit_seqs")
    max_threads = M.get_param_value_from_config("max_threads")
    if not max_threads:
        max_threads = float("Inf")

    if M.metagenomes:
        if wildcards.sample_name in M.metagenomes_name_list:
            threads = M.get_param_value_from_config(
                ["extract_hmm_hit_seqs", "threads_metagenomes"]
            )
    else:
        threads = M.get_param_value_from_config(
            ["extract_hmm_hit_seqs", "threads_genomes"]
        )

    if threads:
        try:
            if int(threads) > float(max_threads):
                return int(max_threads)
            else:
                return int(threads)
        except:
            raise ConfigError(
                f'"threads" must be an integer number. In your config file you provided "{threads}" for '
                f"the number of threads for rule extract_hmm_hit_seqs"
            )
    else:
        return 1


def locate_hmm_profile(hmm_source):
    """Return path to the .hmm file for a given HMM source."""
    if hmm_source in M.internal_hmm_sources:
        import anvio.data.hmm
        return anvio.data.hmm.sources[hmm_source]['model']
    else:
        return os.path.join(M.unique_hmm_source[hmm_source], 'genes.hmm.gz')


def get_hmm_target(hmm_source):
    """Return target string (e.g. 'AA:GENE', 'RNA:CONTIG') for a given HMM source."""
    if hmm_source in M.internal_hmm_sources:
        import anvio.data.hmm
        return anvio.data.hmm.sources[hmm_source]['target']
    else:
        target_path = os.path.join(M.unique_hmm_source[hmm_source], 'target.txt')
        with open(target_path) as f:
            return f.read().strip()


def get_hmm_hits_txt(contigs_db, out_file):
    """Extract the hmm_hits table from a contigs-db."""
    database = db.DB(contigs_db, None, ignore_version=True)
    tables_in_database = database.get_table_names()
    table_name = "hmm_hits"

    if table_name not in tables_in_database:
        column_names = [
            "entry_id", "source", "gene_unique_identifier",
            "gene_callers_id", "gene_name", "gene_hmm_id", "e_value",
        ]
        df = pd.DataFrame(columns=column_names)
        df.to_csv(out_file, sep="\t", index=False, header=True)
    else:
        table_columns = database.get_table_structure(table_name)
        table_content = database.get_table_as_dataframe(
            args_table, columns_of_interest=table_columns, error_if_no_data=False
        )
        u.store_dataframe_as_TAB_delimited_file(table_content, out_file)
    database.disconnect()


def get_extract_done_files(wildcards):
    """Return list of extract_hmm_hit_seqs .done files for Path A samples only.

    Path A: HMMs pre-exist in the contigs DB and no non-empty domtblout
    from anvi-run-hmms was left behind. These samples use group-level HMM
    search and filtering (cat_hmm_hit_seqs → hmmsearch_combined →
    filter_hmm_hits_combined). Path B samples run per-sample instead.
    """
    return [
        os.path.join(
            dirs_dict["HMM_HITS_DIR"], sample,
            f"{wildcards.hmm_source}-dom-hmmsearch",
            "contigs-hmm-extracted.done",
        )
        for sample in M.names_list
        if M.path_is_a[(sample, wildcards.hmm_source)]
    ]


def get_process_hmm_hits_input(wildcards):
    """Return the survivors path based on whether the sample's contigs DB has HMMs.

    Path A (DB has HMMs, no non-empty domtblout from anvi-run-hmms)
        → group-level survivors from filter_hmm_hits_combined.
    Path B (DB lacks HMMs, or anvi-run-hmms ran in a previous attempt)
        → per-sample survivors from filter_hmm_hits_sample.
    """
    if M.path_is_a[(wildcards.sample_name, wildcards.hmm_source)]:
        hmm_key = f"{wildcards.hmm_source}_{wildcards.hmm_name}"
        group = M.hmm_dict[hmm_key]['group']
        return os.path.join(
            dirs_dict["POOLED_HMM_DIR"],
            f"{group}-{wildcards.hmm_source}_survivor_headers.txt",
        )

    return os.path.join(
        dirs_dict["HMM_HITS_DIR"],
        wildcards.sample_name,
        f"{wildcards.hmm_source}-dom-hmmsearch",
        f"{wildcards.sample_name}_{wildcards.hmm_source}_survivor_headers.txt",
    )


# --------------------------------------------------------------------------------
# Step A: Extract sequences per (sample, source) with unique headers, filter partials
# --------------------------------------------------------------------------------

rule extract_hmm_hit_seqs:
    """Extract AA sequences per (sample, source) with double-underscore-delimited headers, removing partial gene calls."""
    output:
        faa=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "{sample_name}_{hmm_source}_hmm_hits.faa",
        ),
        hmm_hits=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "hmm_hits.txt",
        ),
        domtblout=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "hmm.domtable",
        ),
        done=touch(os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "contigs-hmm-extracted.done",
        )),
    log:
        rule_log("extract_hmm_hit_seqs", "extract_hmm_hit_seqs-{sample_name}-{hmm_source}"),
    threads: get_hmm_threads
    params:
        additional_params=M.get_param_value_from_config(
            ["extract_hmm_hit_seqs", "additional_params"]
        ),
        filter_partial=M.get_rule_param(
            "extract_hmm_hit_seqs", "--filter-out-partial-gene-calls"
        ),
    run:
        contigs_db_path = os.path.join(M.contigs_db_name_path_dict[wildcards.sample_name])
        hmm_source = wildcards.hmm_source
        hmmer_output_dir = os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            wildcards.sample_name,
            f"{hmm_source}-dom-hmmsearch",
        )
        os.makedirs(hmmer_output_dir, exist_ok=True)

        # Ensure HMM source is in the contigs DB
        if not M.hmm_source_presence[(wildcards.sample_name, hmm_source)]:
            if hmm_source in M.internal_hmm_sources:
                shell("anvi-run-hmms -c {contigs_db_path} \
                                     --hmmer-program hmmsearch \
                                     --hmmer-output-dir {hmmer_output_dir} \
                                     --installed-hmm-profile {hmm_source} \
                                     --domain-hits-table \
                                     --just-do-it \
                                     -T {threads} >> {log} 2>&1")
            else:
                hmm_dir = os.path.join(M.unique_hmm_source[hmm_source])
                shell("anvi-run-hmms -c {contigs_db_path} \
                                     --hmmer-program hmmsearch \
                                     --hmm-profile-dir {hmm_dir} \
                                     --hmmer-output-dir {hmmer_output_dir} \
                                     --domain-hits-table \
                                     --just-do-it \
                                     -T {threads} >> {log} 2>&1")

        target = get_hmm_target(hmm_source)
        alphabet = target.split(':')[0]
        defline_fmt = f"{wildcards.sample_name}__{hmm_source}__{{gene_name}}__{{gene_callers_id}}"

        # Narrow to the specific HMM gene name(s) declared in hmm_list.txt
        # instead of extracting the entire collection. All names are from
        # M.hmm_dict entries that share this source.
        gene_names = ' '.join(sorted(set(
            v['name'] for v in M.hmm_dict.values()
            if v['source'] == hmm_source
        )))

        if alphabet == 'AA':
            raw_faa = os.path.join(hmmer_output_dir, "raw_hits.faa")
            shell("anvi-get-sequences-for-hmm-hits -c {contigs_db_path} \
                                                     --hmm-sources {hmm_source} \
                                                     --gene-names {gene_names} \
                                                     --get-aa-sequences \
                                                     -o {raw_faa} \
                                                     --defline-format \"{defline_fmt}\" \
                                                     --just-do-it >> {log} 2>&1")
        else:
            raw_faa = os.path.join(hmmer_output_dir, "raw_hits.fna")
            shell("anvi-get-sequences-for-hmm-hits -c {contigs_db_path} \
                                                     --hmm-sources {hmm_source} \
                                                     --gene-names {gene_names} \
                                                     -o {raw_faa} \
                                                     --defline-format \"{defline_fmt}\" \
                                                     --just-do-it >> {log} 2>&1")

        # Optionally filter out partial gene calls
        if params.filter_partial:
            database = db.DB(contigs_db_path, None, ignore_version=True)
            genes_table = database.get_table_as_dataframe('genes_in_contigs')
            database.disconnect()
            partial_ids = set(genes_table[genes_table['partial'] == 1]['gene_callers_id'].astype(str))

            if partial_ids:
                with open(raw_faa) as f_in, open(output.faa + ".tmp", 'w') as f_out:
                    keep = True
                    for line in f_in:
                        if line.startswith('>'):
                            gid = line.strip().rsplit('__', 1)[-1]
                            keep = gid not in partial_ids
                            if keep:
                                f_out.write(line)
                        elif keep:
                            f_out.write(line)
                shutil.move(output.faa + ".tmp", output.faa)
            else:
                shutil.copy(raw_faa, output.faa)
            os.unlink(raw_faa)
        else:
            shutil.move(raw_faa, output.faa)

        # Get hmm_hits.txt
        get_hmm_hits_txt(contigs_db_path, output.hmm_hits)

        # Ensure domtblout exists (created by anvi-run-hmms for Path B;
        # empty touch for Path A where it won't be consumed downstream)
        if not os.path.exists(output.domtblout):
            shell("touch {output.domtblout}")


# --------------------------------------------------------------------------------
# Step B1: Concatenate per-sample FASTA files per (group, source)
# --------------------------------------------------------------------------------

rule cat_hmm_hit_seqs:
    """Concatenate all per-sample FAA files for a (group, source) into one combined FAA."""
    input:
        done_files=get_extract_done_files,
    output:
        combined=os.path.join(
            dirs_dict["POOLED_HMM_DIR"],
            "{group}-{hmm_source}_combined.faa",
        ),
    log:
        rule_log("cat_hmm_hit_seqs", "cat_hmm_hit_seqs-{group}-{hmm_source}"),
    run:
        faa_list = []
        for sample_name in M.names_list:
            if not M.path_is_a[(sample_name, wildcards.hmm_source)]:
                continue
            faa = os.path.join(
                dirs_dict["HMM_HITS_DIR"],
                sample_name,
                f"{wildcards.hmm_source}-dom-hmmsearch",
                f"{sample_name}_{wildcards.hmm_source}_hmm_hits.faa",
            )
            if os.path.exists(faa):
                faa_list.append(faa)

        if not faa_list:
            shell("touch {output.combined}")
        else:
            import tempfile
            tmp = tempfile.NamedTemporaryFile(mode='w', delete=False)
            for f in faa_list:
                tmp.write(f + '\0')
            tmp.close()
            shell(f"xargs -0 cat < {tmp.name} > {output.combined}")
            os.unlink(tmp.name)


# --------------------------------------------------------------------------------
# Step B2: hmmsearch on combined FAA per (group, source)
# --------------------------------------------------------------------------------

rule hmmsearch_combined:
    """Run hmmsearch/nhmmscan once per (group, source) on the combined FASTA."""
    input:
        combined=rules.cat_hmm_hit_seqs.output.combined,
    output:
        domtblout=os.path.join(
            dirs_dict["POOLED_HMM_DIR"],
            "{group}-{hmm_source}_combined.domtblout",
        ),
    log:
        rule_log("hmmsearch_combined", "hmmsearch_combined-{group}-{hmm_source}"),
    threads: M.T("hmmsearch_combined")
    params:
        additional_params=M.get_param_value_from_config(
            ["hmmsearch_combined", "additional_params"]
        ),
    run:
        fasta_path = input.combined

        # Skip if combined FASTA is empty
        total = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
        if total == 0:
            shell("touch {output.domtblout}")
        else:
            hmm_source = wildcards.hmm_source
            hmm_profile = locate_hmm_profile(hmm_source)
            target = get_hmm_target(hmm_source)
            alphabet = target.split(':')[0]
            hmmer_prog = "nhmmscan" if alphabet in ('RNA', 'DNA') else "hmmsearch"

            shell("{hmmer_prog} --cpu {threads} \
                                 --domtblout {output.domtblout} \
                                 -o /dev/null \
                                 {hmm_profile} {fasta_path} >> {log} 2>&1")


# --------------------------------------------------------------------------------
# Step B3: Filter combined domtblout by model coverage per (group, source)
# --------------------------------------------------------------------------------

rule filter_hmm_hits_combined:
    """Filter combined domtblout by model coverage, write survivor headers."""
    input:
        domtblout=rules.hmmsearch_combined.output.domtblout,
    output:
        survivors=os.path.join(
            dirs_dict["POOLED_HMM_DIR"],
            "{group}-{hmm_source}_survivor_headers.txt",
        ),
    log:
        rule_log("filter_hmm_hits_combined", "filter_hmm_hits_combined-{group}-{hmm_source}"),
    threads: M.T("filter_hmm_hits_combined")
    params:
        min_model_coverage=M.get_param_value_from_config(
            ["filter_hmm_hits_combined", "--min-model-coverage"]
        ),
    run:
        domtblout_path = input.domtblout

        # If domtblout is empty (no input sequences), write empty survivors
        if os.path.getsize(domtblout_path) == 0:
            shell("touch {output.survivors}")
            total = 0
            passing = 0
        else:
            # Only keep survivors matching the HMM names declared for this source
            target_names = set(
                v['name'] for v in M.hmm_dict.values()
                if v['source'] == wildcards.hmm_source
            )

            survivors = set()
            total = 0
            with open(domtblout_path) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) < 17:
                        continue
                    if parts[3] not in target_names:
                        continue
                    total += 1
                    target_name = parts[0]
                    hmm_length = int(parts[5])
                    hmm_start = int(parts[15])
                    hmm_stop = int(parts[16])
                    model_coverage = (hmm_stop - hmm_start) / hmm_length
                    if model_coverage >= params.min_model_coverage:
                        survivors.add(target_name)

            with open(output.survivors, 'w') as out:
                for h in sorted(survivors):
                    out.write(h + '\n')

            passing = len(survivors)

        with open(str(log), 'w') as f:
            f.write(f"Model coverage threshold: {params.min_model_coverage}\n")
            f.write(f"Total hits in domtblout: {total}\n")
            f.write(f"Hits passing filter: {passing}\n")


# --------------------------------------------------------------------------------
# Step B4: Filter per-sample domtblout by model coverage (Path B)
# --------------------------------------------------------------------------------

rule filter_hmm_hits_sample:
    """Filter per-sample domtblout by model coverage, write survivor headers."""
    input:
        domtblout=rules.extract_hmm_hit_seqs.output.domtblout,
    output:
        survivors=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}-dom-hmmsearch",
            "{sample_name}_{hmm_source}_survivor_headers.txt",
        ),
    log:
        rule_log("filter_hmm_hits_sample", "filter_hmm_hits_sample-{sample_name}-{hmm_source}"),
    threads: M.T("filter_hmm_hits_sample")
    params:
        min_model_coverage=M.get_param_value_from_config(
            ["filter_hmm_hits_sample", "--min-model-coverage"]
        ),
    run:
        domtblout_path = input.domtblout

        if os.path.getsize(domtblout_path) == 0:
            shell("touch {output.survivors}")
            total = 0
            passing = 0
        else:
            # Only keep survivors matching the HMM names declared for this source
            target_names = set(
                v['name'] for v in M.hmm_dict.values()
                if v['source'] == wildcards.hmm_source
            )

            survivors = set()
            total = 0
            with open(domtblout_path) as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    parts = line.split()
                    if len(parts) < 17:
                        continue
                    if parts[3] not in target_names:
                        continue
                    total += 1
                    target_name = parts[0]
                    hmm_length = int(parts[5])
                    hmm_start = int(parts[15])
                    hmm_stop = int(parts[16])
                    model_coverage = (hmm_stop - hmm_start) / hmm_length
                    if model_coverage >= params.min_model_coverage:
                        survivors.add(target_name)

            with open(output.survivors, 'w') as out:
                for h in sorted(survivors):
                    out.write(h + '\n')

            passing = len(survivors)

        with open(str(log), 'w') as f:
            f.write(f"Model coverage threshold: {params.min_model_coverage}\n")
            f.write(f"Total hits in domtblout: {total}\n")
            f.write(f"Hits passing filter: {passing}\n")


# --------------------------------------------------------------------------------
# Step C: Per (sample, source, name) — survivor-filtered AA, NT, EGC extraction
# --------------------------------------------------------------------------------

rule process_hmm_hits:
    """Extract AA/NT fastas and external-gene-calls for survivors only."""
    input:
        survivors=get_process_hmm_hits_input,
    output:
        aa_fasta=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}",
            "{hmm_name}",
            "{sample_name}-{hmm_name}-hmm_hits_renamed.faa",
        ),
        nt_fasta=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}",
            "{hmm_name}",
            "{sample_name}-{hmm_name}-hmm_hits_renamed.fna",
        ),
        egc_renamed=os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}",
            "{hmm_name}",
            "{sample_name}-{hmm_name}-external_gene_calls_renamed.tsv",
        ),
        done=touch(os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            "{sample_name}",
            "{hmm_source}",
            "{hmm_name}",
            "{sample_name}-{hmm_name}-processed.done",
        )),
    log:
        rule_log(
            "process_hmm_hits",
            "process_hmm_hits-{sample_name}-{hmm_source}-{hmm_name}",
        ),
    threads: M.T("process_hmm_hits")
    run:
        contigs_db = os.path.join(M.contigs_db_name_path_dict[wildcards.sample_name])
        hmm_source = wildcards.hmm_source
        hmm_name = wildcards.hmm_name
        sample_name = wildcards.sample_name

        fasta_output_dir = os.path.join(
            dirs_dict["HMM_HITS_DIR"],
            sample_name,
            hmm_source,
            hmm_name,
        )
        os.makedirs(fasta_output_dir, exist_ok=True)

        # Determine group from hmm_dict
        group = None
        for value in M.hmm_dict.values():
            if value['source'] == hmm_source and value['name'] == hmm_name:
                group = value['group']
                break
        if group is None:
            raise ConfigError(f"No group found for HMM source '{hmm_source}', name '{hmm_name}'")

        # Read survivor headers from input (resolved by get_process_hmm_hits_input)
        survivor_path = input.survivors

        # Determine which path we're on – must match get_process_hmm_hits_input's logic
        is_path_b = not M.path_is_a[(sample_name, hmm_source)]

        survivors = []
        survivor_gene_callers_ids = []
        if os.path.getsize(survivor_path) == 0:
            pass
        elif is_path_b:
            # Path B: survivors are bare gene caller IDs already filtered
            # to this HMM name by filter_hmm_hits_sample.
            with open(survivor_path) as f:
                for line in f:
                    gid = line.strip()
                    if gid:
                        survivor_gene_callers_ids.append(gid)
                        survivors.append(f"{sample_name}__{hmm_source}__{hmm_name}__{gid}")
        else:
            # Path A: survivors are full double-underscore-delimited headers
            prefix = f"{sample_name}__{hmm_source}__{hmm_name}__"
            with open(survivor_path) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(prefix):
                        gid = line.rsplit('__', 1)[-1]
                        survivors.append(line)
                        survivor_gene_callers_ids.append(gid)

        if not survivors:
            # No survivors: write empty files
            shell("touch {output.aa_fasta} {output.nt_fasta}")
            col_names = ["gene_callers_id", "contig", "start", "stop",
                         "direction", "partial", "call_type", "source",
                         "version", "aa_sequence"]
            pd.DataFrame(columns=col_names).to_csv(output.egc_renamed, sep="\t", index=False)
        else:
            # Step C1: Extract AA sequences with defline format matching survivor headers
            aa_fmt = f"{sample_name}__{hmm_source}__{hmm_name}__{{gene_callers_id}}"
            aa_all = os.path.join(fasta_output_dir, "aa_all.faa")
            shell("anvi-get-sequences-for-hmm-hits -c {contigs_db} \
                                                     --hmm-sources {hmm_source} \
                                                     --gene-names {hmm_name} \
                                                     --get-aa-sequences \
                                                     -o {aa_all} \
                                                     --defline-format \"{aa_fmt}\" \
                                                     --just-do-it >> {log} 2>&1")

            # Filter AA FASTA to survivor headers only
            survivor_set = set(survivors)
            with open(aa_all) as f_in, open(output.aa_fasta, 'w') as f_out:
                write = True
                for line in f_in:
                    if line.startswith('>'):
                        header = line[1:].strip()
                        write = header in survivor_set
                        if write:
                            f_out.write(line)
                    elif write:
                        f_out.write(line)
            os.unlink(aa_all)

            # Step C2: Extract NT sequences and EGC for surviving gene callers
            gid_str = ",".join(survivor_gene_callers_ids)
            raw_nt = os.path.join(fasta_output_dir, "raw_nt.fna")
            raw_egc = os.path.join(fasta_output_dir, "raw_egc.tsv")
            shell("anvi-get-sequences-for-gene-calls -c {contigs_db} \
                                                       --gene-caller-ids {gid_str} \
                                                       --external-gene-calls {raw_egc} \
                                                       -o {raw_nt} >> {log} 2>&1")

            # Post-process NT FASTA headers and EGC with matching numeric identifiers
            egc = pd.read_csv(raw_egc, delim_whitespace=True, index_col=False)
            contig_to_gid = dict(zip(egc['contig'].astype(str),
                                    egc['gene_callers_id'].astype(str)))

            with open(raw_nt) as f_in, open(output.nt_fasta, 'w') as f_out:
                for line in f_in:
                    if line.startswith('>'):
                        raw_id = line[1:].strip()
                        gid = contig_to_gid.get(raw_id, raw_id)
                        f_out.write(f">{sample_name}__{hmm_source}__{hmm_name}__{gid}\n")
                    else:
                        f_out.write(line)
            os.unlink(raw_nt)

            # Post-process EGC: add header column with global identifier,
            # keep original contig column for bam-based coverage joins
            egc['header'] = [f"{sample_name}__{hmm_source}__{hmm_name}__{gid}"
                             for gid in egc['gene_callers_id'].astype(str)]
            egc.to_csv(output.egc_renamed, sep="\t", index=False)
            os.unlink(raw_egc)
