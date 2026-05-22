# profile-mode with read recruitment


rule anvi_summarize:
    """Summarize merged profile for hmm_hits"""
    input:
        profileDB=ancient(os.path.join(dirs_dict["MERGE_DIR"], "{group}", "PROFILE.db")),
    output:
        done=touch(
            os.path.join(dirs_dict["MERGE_DIR"], "{group}", "{group}_summarize.done")
        ),
    log:
        rule_log("anvi_summarize", "anvi_summarize_{group}"),
    threads: M.T("anvi_summarize")
    params:
        contigsDB=ancient(
            os.path.join(
                dirs_dict["HOME"], "METAGENOMICS_WORKFLOW", "03_CONTIGS", "{group}.db"
            )
        ),
        profileDB=os.path.join(
            dirs_dict["HOME"],
            "METAGENOMICS_WORKFLOW",
            "06_MERGED",
            "{group}",
            "PROFILE.db",
        ),
        output_dir=os.path.join(
            dirs_dict["HOME"], "METAGENOMICS_WORKFLOW", "07_SUMMARY", "{group}"
        ),
    shell:
        "anvi-summarize -c {params.contigsDB} -p {params.profileDB} -o {params.output_dir} -C DEFAULT --init-gene-coverages --just-do-it >> {log} 2>&1"


rule make_anvio_state_file:
    """Make a state file customized for EcoPhylo workflow interactive interface"""
    input:
        source=M.get_target_files_make_anvio_state_file(),
    output:
        state_file=os.path.join(
            dirs_dict["MERGE_DIR"],
            "{group}",
            "{group}_ECOPHYLO_WORKFLOW_state.json",
        ),
    log:
        rule_log("make_anvio_state_file", "make_anvio_state_file_{group}"),
    threads: M.T("make_anvio_state_file")
    params:
        tax_data_final=os.path.join(
            dirs_dict["MISC_DATA"], "{group}", "{group}_scg_taxonomy_data.tsv"
        ),
        misc_data_final=os.path.join(
            dirs_dict["MISC_DATA"], "{group}", "{group}_misc.tsv"
        ),
    run:
        # Read in misc data headers for layer_order
        with open(params.misc_data_final) as f:
            lines = f.read()
            first = lines.split("\n", 1)[0]
        misc_layers_list = first.split("\t")
        state_dict = {}
        # basics
        state_dict["version"] = "3"
        state_dict["tree-type"] = "phylogram"
        state_dict["current-view"] = "single"
        # height and width
        # FIXME: It's unclear to me how the interactive interface determines
        # height and width of a tree when the input value is 0. There has to
        # be some kind of calculation to determine the tree shape in the backend
        # of the interface because even after I export a "default" state file
        # the height and width are still "0". However, if you change the height and width
        # values within the interface to "" the tree will disappear. I need to sort this
        # out eventually to have a clean way of changing the tree shape to
        # match the dimensions of the number of SCGs vs metagenomes.
        # num_tree_tips = pd.read_csv(input.num_tree_tips, \
        #                             sep="\t", \
        #                             index_col=None)
        # layer-orders
        first_layers = ["__parent__", "length", "gc_content"]
        metagenomes = []
        for metagenome in M.sample_names_for_mapping_list:
            metagenomes.append(metagenome)
        layer_order = first_layers + metagenomes + misc_layers_list
        # Read in misc data headers for layer_order
        if os.path.isfile(params.tax_data_final):
            with open(params.tax_data_final) as f:
                lines = f.read()
                first = lines.split("\n", 1)[0]
                scg_taxonomy_layers_list = first.split("\t")
            layer_order.extend(scg_taxonomy_layers_list)
        state_dict["layer-order"] = layer_order
        # layers
        layers_dict = {}
        metagenome_layers_dict = {}
        metagenome_attributes = {
            "color": "#000000",
            "height": "180",
            "margin": "15",
            "type": "bar",
            "color-start": "#FFFFFF",
        }
        for metagenome in metagenomes:
            metagenome_layers_dict[str(metagenome)] = metagenome_attributes
        layer_attributes_parent = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF",
        }
        length = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF",
        }
        gc_content = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF",
        }
        identifier = {
            "color": "#000000",
            "height": "0",
            "margin": "15",
            "type": "color",
            "color-start": "#FFFFFF",
        }
        percent_identity = {
            "color": "#000000",
            "height": "180",
            "margin": "15",
            "type": "line",
            "color-start": "#FFFFFF",
        }
        layers_dict.update(metagenome_layers_dict)
        layers_dict["__parent__"] = layer_attributes_parent
        layers_dict["length"] = length
        layers_dict["gc_content"] = gc_content
        layers_dict["identifier"] = identifier
        layers_dict["percent_identity"] = percent_identity
        state_dict["layers"] = layers_dict
        # views
        views_dict = {}
        single_dict = {}
        mean_coverage_dict = {}
        percent_identity = {
            "normalization": "none",
            "min": {"value": "90", "disabled": "false"},
            "max": {"value": "100", "disabled": "false"},
        }
        cluster_size = {"normalization": "none"}
        single_dict["percent_identity"] = percent_identity
        mean_coverage_dict["percent_identity"] = percent_identity
        mean_coverage_dict["cluster_size"] = cluster_size
        views_dict["single"] = single_dict
        views_dict["mean_coverage"] = mean_coverage_dict
        state_dict["views"] = views_dict
        # samples-layer-order
        samples_layers_dict = {
            "default": {
                "num_INDELs_reported": {
                    "height": 0,
                },
                "total_reads_kept": {
                    "height": 0,
                },
                "num_SCVs_reported": {
                    "height": 0,
                },
                "num_SNVs_reported": {
                    "height": 0,
                },
                "total_reads_mapped": {
                    "height": 0,
                },
            }
        }
        state_dict["samples-layers"] = samples_layers_dict
        with open(output.state_file, "w") as outfile:
            json.dump(state_dict, outfile, indent=4)


rule anvi_import_everything_metagenome:
    """Import state file, phylogenetic tree, AND misc data to interactive interface
If samples.txt is NOT provided then we will make an Ad Hoc profileDB for the tree to import misc data

"""
    input:
        tree=rules.rename_tree_tips.output.tree,
        state=rules.make_anvio_state_file.output.state_file,
        profileDB=ancient(os.path.join(dirs_dict["MERGE_DIR"], "{group}", "PROFILE.db")),
    output:
        done=touch(
            os.path.join(
                dirs_dict["MERGE_DIR"],
                "{group}",
                "{group}_state_imported_profile.done",
            )
        ),
    log:
        rule_log("anvi_import_everything_metagenome", "anvi_import_state_{group}"),
    threads: M.T("anvi_import_state")
    params:
        tax_data_final=rules.anvi_estimate_scg_taxonomy.params.tax_data_final,
        profileDB=os.path.join(dirs_dict["MERGE_DIR"], "{group}", "PROFILE.db"),
        tree_profileDB=os.path.join(dirs_dict["TREES"], "{group}", "{group}-PROFILE.db"),
        misc_data=rules.make_misc_data.output.misc_data_final,
    run:
        state = os.path.join(
            dirs_dict["MERGE_DIR"],
            f"{wildcards.group}",
            f"{wildcards.group}_ECOPHYLO_WORKFLOW_state.json",
        )
        shell("echo -e 'Step 1: anvi-import-state:\n' >> {log}")
        shell(
            "anvi-import-state -p {params.profileDB} -s {state} -n default >> {log} 2>&1"
        )
        shell("echo -e '' >> {log}")
        shell("echo -e 'Step 2: anvi-import-items-order:\n' >> {log}")
        shell(
            "anvi-import-items-order -p {params.profileDB} -i {input.tree} --name {wildcards.group}_tree >> {log} 2>&1"
        )
        shell("echo -e '' >> {log}")
        shell("echo -e 'Step 3: anvi-import-misc-data:\n' >> {log}")
        shell(
            "anvi-import-misc-data -p {params.profileDB} --target-data-table items {params.misc_data} --just-do-it >> {log} 2>&1"
        )
        shell("echo -e '' >> {log}")
        if os.path.isfile(params.tax_data_final):
            shell("echo -e 'Step 4: anvi-import-misc-data:\n' >> {log}")
            shell(
                "anvi-import-misc-data -p {params.profileDB} --target-data-table items {params.tax_data_final} --just-do-it >> {log} 2>&1"
            )
            shell("echo -e '' >> {log}")
