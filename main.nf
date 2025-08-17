#!/usr/bin/env nextflow

/*
 * main.nf – Nextflow implementation of the Population Genomics
 * workflow.  This pipeline is inspired by the Huang‑lab
 * PopulationGenomics repository and orchestrates annotation,
 * feature extraction, filtering and summary steps using ANNOVAR,
 * InterVar, bcftools, R and Python.  The code below follows the
 * logical ordering of the original Bash wrappers (s1–s9) and
 * groups related operations into Nextflow processes.  Before
 * running the pipeline you must configure the paths to third‑party
 * tools and data files in `nextflow.config` or on the command line
 * with `--param value` pairs.
 */

// -----------------------------------------------------------------------------
// Input definitions
// -----------------------------------------------------------------------------

// Gather all gzipped VCFs from the input directory.  Each VCF
// represents a cohort or sample to be annotated.  You can adjust
// the glob to match your file naming convention.
Channel.fromPath("${params.input_dir}/*.vcf.gz", checkIfExists: true).set { vcf_files }

// -----------------------------------------------------------------------------
// s1: Annotate variants with ANNOVAR and InterVar
// -----------------------------------------------------------------------------

/*
 * The `annotate_variants` process performs three tasks for each
 * compressed VCF:  (1) decompress the file to a plain VCF, (2)
 * annotate the variants with ANNOVAR using multiple reference
 * databases and (3) interpret the pathogenicity of each variant
 * using InterVar.  The output consists of the decompressed VCF,
 * the ANNOVAR multi‑anno table and the InterVar pathogenicity table.
 *
 * This step corresponds to the original `s1_run_intjob.sh` script
 * from the repository.  Parameters such as
 * genome build, protocols and database paths can be customised via
 * `nextflow.config`.
 */
process annotate_variants {
    label 'large'
    publishDir "${params.outdir}/s1_annotation", mode: 'copy'

    input:
        path vcf_file from vcf_files

    output:
        tuple val(sample_id),
              path("${sample_id}.variants.vcf"),
              path("${sample_id}.${params.genome_build}_multianno.txt"),
              path("${sample_id}.InterVar.tsv")
              into annotated_results_ch

    script:
        // Derive a sample identifier by stripping the `.vcf.gz` suffix.
        def sample_id = vcf_file.baseName.replaceAll(/\.vcf\.gz\$/, '')

        // Compose the command string.  Use a here–document to
        // accommodate multiple steps.  Note that the number of
        // threads for ANNOVAR is determined by the process CPU
        // allocation (task.cpus).
        """
        set -euo pipefail

        # Decompress the input VCF
        zcat '${vcf_file}' > '${sample_id}.variants.vcf'

        # Annotate with ANNOVAR
        perl '${params.annovar_home}/table_annovar.pl' \
            '${sample_id}.variants.vcf' \
            '${params.annovar_db}' \
            -buildver '${params.genome_build}' \
            -out '${sample_id}' \
            -remove \
            -protocol refGene,clinvar_20220320,gnomad30_genome,gnomad211_genome,gnomad30_exome,gnomad211_exome \
            -operation g,f,f,f,f,f \
            -nastring . \
            -vcfinput \
            -thread ${task.cpus}

        # Run InterVar to evaluate variant pathogenicity.  InterVar
        # outputs a file named "<prefix>.InterVar.tsv".  We rename
        # this to a predictable name for downstream processes.  The
        # "-o" flag defines the prefix of the output files.
        python3 '${params.intervar_script}' \
            -i '${sample_id}.${params.genome_build}_multianno.txt' \
            -d '${params.annovar_db}' \
            -b '${params.genome_build}' \
            -o '${sample_id}.InterVar'

        # InterVar writes '${sample_id}.InterVar.tsv' directly; no
        # renaming is required.
        """
}

// Split the annotated results into separate channels for later steps.
annotated_results_ch.map { sample_id, vcf, multianno, intervar -> multianno }  .set { multianno_files }
annotated_results_ch.map { sample_id, vcf, multianno, intervar -> intervar }   .set { intervar_files }
annotated_results_ch.map { sample_id, vcf, multianno, intervar -> vcf }        .set { variants_vcf_files }

// -----------------------------------------------------------------------------
// s2: Merge ANNOVAR and InterVar annotations
// -----------------------------------------------------------------------------

/*
 * After all samples have been annotated, merge the resulting
 * multi‑anno tables with the corresponding InterVar tables.  The
 * supplied R script (`s2_Merge_annotation.R`) scans the working
 * directory for `*.multianno.txt` and `*.InterVar.tsv` files,
 * performs the merge and writes a combined table.  The process is
 * triggered once, after all inputs are available, via the use of
 * `collect()`.
 */
process merge_annotations {
    label 'small'
    publishDir "${params.outdir}/s2_merge", mode: 'copy'

    // Use collect to gather all ANNOVAR and InterVar tables into
    // lists.  The process will run once and wait for all upstream
    // files to arrive.
    input:
        file annovar_list from multianno_files.collect()
        file intervar_list from intervar_files.collect()

    output:
        file "s2_Merge_intervar_annovar_multianno.txt" into merged_annotation

    script:
        """
        set -euo pipefail
        # Copy the input files into the current working directory (Nextflow
        # stages them automatically).  The R script expects to find
        # `*.multianno.txt` and `*.InterVar.tsv` in the working directory.

        # Run the merge script.  It writes the merged table to
        # 's2_Merge_intervar_annovar_multianno.txt'.
        Rscript '${params.scripts_dir}/preprocessing/s2_Merge_annotation.R'
        """
}

// -----------------------------------------------------------------------------
// s3: Extract PLP and PTV variants
// -----------------------------------------------------------------------------

/*
 * The third stage extracts lists of pathogenic/likely pathogenic
 * variants (PLPs) and protein truncating variants (PTVs) from the
 * merged annotation table using the provided R script
 * (`s3_Extract_PLP_PTV_variants.R`).  The
 * script reads the merged table and produces several files,
 * including:
 *
 *   - `s2_Merge_intervar_annovar_multianno_simple.txt` – simplified
 *     annotation used for second‑round VCF filtering.
 *   - `s2_all_PLP_list.txt` – list of PLP variant coordinates.
 *   - `s3_ACMG32_cancer_gene.predis.plp.txt` – list of PLPs in
 *     predefined cancer genes.
 *   - `s3_ACMG32_cancer_gene_truncations.txt` – list of PTVs.
 *
 * These outputs feed into the downstream filtering processes.
 */
process extract_plp_ptv {
    label 'small'
    publishDir "${params.outdir}/s3_extract", mode: 'copy'

    input:
        file merged_file from merged_annotation

    output:
        tuple file("s2_Merge_intervar_annovar_multianno_simple.txt"),
              file("s2_all_PLP_list.txt"),
              file("s3_ACMG32_cancer_gene.predis.plp.txt"),
              file("s3_ACMG32_cancer_gene_truncations_regions.txt"),
              file("s3_ACMG32_cancer_gene_truncations.txt")
              into plp_ptv_outputs

    script:
        """
        set -euo pipefail
        # Run the R script to extract PLP and PTV lists
        Rscript '${params.scripts_dir}/feature_extraction/s3_Extract_PLP_PTV_variants.R'
        """
}

// Unpack the outputs into individual channels
plp_ptv_outputs.map { simple_table, plp_list, predis_plp, trunc_regions, trunc_list -> simple_table } .set { simple_table_ch }
plp_ptv_outputs.map { simple_table, plp_list, predis_plp, trunc_regions, trunc_list -> plp_list }     .set { plp_list_ch }
plp_ptv_outputs.map { simple_table, plp_list, predis_plp, trunc_regions, trunc_list -> trunc_regions } .set { trunc_regions_ch }
plp_ptv_outputs.map { simple_table, plp_list, predis_plp, trunc_regions, trunc_list -> trunc_list }    .set { trunc_list_ch }

// -----------------------------------------------------------------------------
// s4: Filter VCF by PLP variants (two‑stage filtering)
// -----------------------------------------------------------------------------

/*
 * The first filtering step uses bcftools to retain only variants
 * present in the PLP list (`s2_all_PLP_list.txt`) from the original
 * decompressed VCFs.  The second step
 * applies an additional filter using the simplified annotation table
 * produced in s3 via a Python script (`General_2nd_filterVCF.py`)
 *.  Both steps are applied per sample.
 */

process filter_plp_bcftools {
    label 'small'
    publishDir "${params.outdir}/s4_filter_plp", mode: 'copy'

    input:
        tuple val(sample_id), path(variants) from annotated_results_ch
        file plp_list               from plp_list_ch

    output:
        tuple val(sample_id), path("${sample_id}.PLP.raw.vcf") into plp_raw_vcf

    script:
        """
        set -euo pipefail
        # Extract only PLP variants from the sample VCF.  The PLP list
        # contains variant coordinates in the format expected by
        # bcftools (chr:pos).  bcftools view writes a VCF with
        # matching records.
        bcftools view -R '${plp_list}' -Ov -o '${sample_id}.PLP.raw.vcf' '${variants}'
        """
}

process filter_plp_second {
    label 'small'
    publishDir "${params.outdir}/s4_filter_plp", mode: 'copy'

    input:
        tuple val(sample_id), path(plp_vcf)    from plp_raw_vcf
        file simple_table                      from simple_table_ch

    output:
        tuple val(sample_id), path("${sample_id}.PLP.filtered.vcf") into plp_filtered_vcf

    script:
        """
        set -euo pipefail
        # Apply the second filter: remove variants that do not pass
        # additional criteria stored in the simplified annotation
        # table.  The Python script retains entries whose variant
        # key appears in the table (column 'VariantKey').
        python3 '${params.scripts_dir}/feature_extraction/General_2nd_filterVCF.py' \
            -i '${plp_vcf}' \
            -t '${simple_table}' \
            -o '${sample_id}.PLP.filtered.vcf' \
            -k VariantKey
        """
}

// -----------------------------------------------------------------------------
// s5: Retain rare sites by allele frequency
// -----------------------------------------------------------------------------

/*
 * To focus on ultra‑rare variation, we filter each VCF by a
 * maximum allele frequency threshold using bcftools.  The threshold
 * is configurable via `params.max_af` and mirrors the original
 * `--max-af 0.0005` option.
 */
process filter_by_af {
    label 'small'
    publishDir "${params.outdir}/s5_filter_af", mode: 'copy'

    input:
        tuple val(sample_id), path(vcf_file) from plp_filtered_vcf

    output:
        tuple val(sample_id), path("${sample_id}.rare.vcf") into rare_af_vcf

    script:
        """
        set -euo pipefail
        # Filter by allele frequency.  bcftools will drop records
        # whose INFO/AF field exceeds the specified threshold.  The
        # output VCF is sorted to ensure correct ordering.
        bcftools view --max-af ${params.max_af} -Ov -o '${sample_id}.rare.vcf' '${vcf_file}'
        """
}

// -----------------------------------------------------------------------------
// s6: Filter VCF by PTV variants (two‑stage filtering)
// -----------------------------------------------------------------------------

/*
 * This section mirrors s4 but operates on protein truncating
 * variants.  First we extract variants listed in
 * `s3_ACMG32_cancer_gene_truncations.txt`, then apply the same
 * second‑round filtering using the simplified annotation table.
 */
process filter_ptv_bcftools {
    label 'small'
    publishDir "${params.outdir}/s6_filter_ptv", mode: 'copy'

    input:
        tuple val(sample_id), path(variants) from annotated_results_ch
        file trunc_list                         from trunc_list_ch

    output:
        tuple val(sample_id), path("${sample_id}.PTV.raw.vcf") into ptv_raw_vcf

    script:
        """
        set -euo pipefail
        bcftools view -R '${trunc_list}' -Ov -o '${sample_id}.PTV.raw.vcf' '${variants}'
        """
}

process filter_ptv_second {
    label 'small'
    publishDir "${params.outdir}/s6_filter_ptv", mode: 'copy'

    input:
        tuple val(sample_id), path(ptv_vcf)   from ptv_raw_vcf
        file simple_table                     from simple_table_ch

    output:
        tuple val(sample_id), path("${sample_id}.PTV.filtered.vcf") into ptv_filtered_vcf

    script:
        """
        set -euo pipefail
        python3 '${params.scripts_dir}/feature_extraction/General_2nd_filterVCF.py' \
            -i '${ptv_vcf}' \
            -t '${simple_table}' \
            -o '${sample_id}.PTV.filtered.vcf' \
            -k VariantKey
        """
}

// -----------------------------------------------------------------------------
// s7: Obtain rare PTV variants
// -----------------------------------------------------------------------------

/*
 * Combine the allele frequency filter with the PTV filter to
 * isolate ultra‑rare protein truncating variants.  bcftools is used
 * again with `--max-af` and `-R` options.
 */
// Before extracting rare PTVs we join the allele‑frequency filtered
// VCFs with the PTV‑filtered VCFs on the sample identifier.  This
// produces a channel of triples: (sample_id, rare_af_vcf, ptv_vcf).
rare_af_vcf.join(ptv_filtered_vcf, by: 0).map { sid, rareFile, ptvFile -> tuple(sid, rareFile, ptvFile) }.set { rare_ptv_pairs }

process rare_ptv {
    label 'small'
    publishDir "${params.outdir}/s7_rare_ptv", mode: 'copy'

    input:
        tuple val(sample_id), path(rare_af), path(ptv_vcf) from rare_ptv_pairs

    output:
        tuple val(sample_id), path("${sample_id}.PTV.rare.filtered.vcf") into rare_ptv_vcf

    script:
        """
        set -euo pipefail
        # Extract rare PTV variants by intersecting the AF‑filtered
        # VCF with the PTV list.  The `-R` option tells bcftools to
        # keep only positions present in the PTV VCF, while
        # `--max-af` ensures that only variants below the allele
        # frequency threshold are retained.
        bcftools view --max-af ${params.max_af} -R '${ptv_vcf}' -Ov -o '${sample_id}.PTV.rare.filtered.vcf' '${rare_af}'
        """
}

// -----------------------------------------------------------------------------
// s8: Generate summary tables
// -----------------------------------------------------------------------------

/*
 * This stage produces high‑level summary tables for the cohort,
 * combining the PLP and PTV filtered VCFs with gene and patient
 * meta‑data.  It calls the Python script
 * `s8_Generate_summary_df.py` which
 * aggregates variant information at the gene and sample level.  The
 * required inputs include the PLP and PTV filtered VCFs, a gene
 * information table, the gene panel file, patient meta‑data and a
 * BED file specifying genomic regions.  These are supplied via
 * Nextflow parameters.  Multiple summary tables are produced and
 * published to the `s8_summary` directory.
 */
// Join the PLP and PTV filtered VCFs by sample identifier to
// coordinate their inputs into the summary stage.  The join
// produces tuples `(sid, plp_vcf, ptv_vcf)`.
plp_filtered_vcf.join(ptv_filtered_vcf, by: 0).map { sid, plpFile, ptvFile -> tuple(sid, plpFile, ptvFile) }.set { summary_input }

process generate_summary {
    label 'medium'
    publishDir "${params.outdir}/s8_summary", mode: 'copy'

    input:
        tuple val(sample_id), path(plp_vcf), path(ptv_vcf) from summary_input

    output:
        tuple val(sample_id), path("${sample_id}_summary.tsv") into summary_tables

    script:
        """
        set -euo pipefail
        # Run the summary generation script.  It writes multiple
        # summary files; we capture the main summary for the sample
        # as '<sample_id>_summary.tsv'.
        python3 '${params.scripts_dir}/summary/s8_Generate_summary_df.py' \
            --plp_vcf '${plp_vcf}' \
            --ptv_vcf '${ptv_vcf}' \
            --gene_info '${params.gene_info}' \
            --gene_panel '${params.gene_panel}' \
            --patient_info '${params.patient_info}' \
            --region_bed '${params.region_bed}' \
            --sample_id '${sample_id}'

        # Copy the primary summary file to a predictable name.  This
        # assumes the script writes '<sample_id>_summary.tsv'.
        cp '${sample_id}_summary.tsv' '${sample_id}_summary.tsv'
        """
}

// -----------------------------------------------------------------------------
// s9: Calculate carrier frequencies
// -----------------------------------------------------------------------------

/*
 * The final step summarises carrier frequencies across the cohort by
 * ancestry, gene and other groupings using
 * `s9_Cal_carrier_Freq.py`.  It consumes
 * the summary table(s) produced in s8 and the patient meta‑data.  A
 * dictionary of aggregated tables is returned and saved to disk.
 */
process calc_carrier_freq {
    label 'small'
    publishDir "${params.outdir}/s9_carrier_freq", mode: 'copy'

    input:
        tuple val(sample_id), path(summary_table) from summary_tables

    output:
        path("${sample_id}_carrier_frequencies.tsv") into carrier_freq_tables

    script:
        """
        set -euo pipefail
        python3 '${params.scripts_dir}/summary/s9_Cal_carrier_Freq.py' \
            --summary '${summary_table}' \
            --patient_info '${params.patient_info}' \
            --sample_id '${sample_id}'
        cp '${sample_id}_carrier_frequencies.tsv' '${sample_id}_carrier_frequencies.tsv'
        """
}

// -----------------------------------------------------------------------------
// Workflow definition
// -----------------------------------------------------------------------------

/*
 * Define the workflow ordering.  Nextflow implicitly wires the
 * processes together based on the channels defined above.  The
 * ordering here is primarily for readability; the dataflow
 * semantics ensure that each process runs as soon as its inputs
 * become available.
 */
workflow {
    // Stage s1: annotation
    annotated_results_ch

    // Stage s2: merge annotations
    merged_annotation

    // Stage s3: extract PLP/PTV lists
    plp_ptv_outputs

    // Stage s4: PLP filtering
    plp_raw_vcf
    plp_filtered_vcf

    // Stage s5: allele frequency filter
    rare_af_vcf

    // Stage s6–s7: PTV filtering and rare PTV extraction
    ptv_raw_vcf
    ptv_filtered_vcf
    rare_ptv_vcf

    // Stage s8: summary generation
    summary_tables

    // Stage s9: carrier frequency
    carrier_freq_tables
}
