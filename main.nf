include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main.nf'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main.nf'
include { METHYLDACKEL_EXTRACT } from './modules/nf-core/methyldackel/extract/main.nf' 
include { METHYLDACKEL_EXTRACT as METHYLDACKEL_EXTRACT_SIMULATED } from './modules/nf-core/methyldackel/extract/main.nf'

process GENERATE_CPG_TABLE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    path(fasta)

    output:
    tuple val(meta), path('*_cpgtable.txt'), emit: cpgtable

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${projectDir}/scripts/simulate_3gen.py \\
    --bam_file ${bam} \\
    --cpgtable_file ${prefix}_cpgtable.txt \\
    --is_only_cpgtable \\
    --reference_genome ${fasta} \\
    --bed_file ${bed} \\
    --mode regionbed \\
    --output_bam ${prefix}_simulated.bam
    """
}

process SIMULATE_SEQUENCES {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bai), path(bed), path(cpgtable), val(num_turns)
    path(fasta)

    output:
    tuple val(meta), path('*_simulated.bam'), path('*_simulated.bam.bai'), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${projectDir}/scripts/simulate_3gen.py \\
    --bam_file ${bam} \\
    --cpgtable_file ${cpgtable} \\
    --reference_genome ${fasta} \\
    --num_turns ${num_turns} \\
    --bed_file ${bed} \\
    --mode regionbed \\
    --output_bam ${prefix}_simulated.bam
    """
}

workflow {
    // Parse input samplesheet.csv
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.label]
            [meta, file(row.bam)]
        }
        .groupTuple(by: 0)
        .set { ch_bams }

    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.label]
            [meta, file(row.bed)]
        }
        .unique()  // Keep only unique label-bed pairs
        .set { ch_beds }

    Channel
        .value([[:], file(params.reference)])
        .set { ch_fasta }

    def fai_path = "${params.reference}.fai"
    if (file(fai_path).exists()) {
        ch_fai = Channel.value([[:], file(fai_path)])
    } else {
        SAMTOOLS_FAIDX(ch_fasta, [[:], []])
        ch_fai = SAMTOOLS_FAIDX.out.fai
    }

    SAMTOOLS_MERGE(ch_bams, ch_fasta, ch_fai)

    SAMTOOLS_SORT(SAMTOOLS_MERGE.out.bam, ch_fasta)
    ch_merged_bam = SAMTOOLS_SORT.out.bam

    SAMTOOLS_INDEX(ch_merged_bam)
    ch_merged_bam_bai = SAMTOOLS_INDEX.out.bai

    METHYLDACKEL_EXTRACT(
        ch_merged_bam.join(ch_merged_bam_bai).map {meta, bam, bai -> tuple(meta, bam, bai)},
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai },
    )

    ch_merged_bam
        .join(ch_merged_bam_bai)
        .join(ch_beds)
        .set { ch_merged_bams_beds }

    GENERATE_CPG_TABLE(
        ch_merged_bams_beds, 
        ch_fasta.map { meta, fasta -> fasta }
    )

    Channel
        .fromList(params.num_turns_list)
        .combine(Channel.of(1..params.num_replicates))
        .map { num_turns, rep -> 
            def meta = [id: "turns_${num_turns}_rep_${rep}"]
            [meta, num_turns]
        }
        .set { ch_num_turns }

    // Combine channels
    ch_merged_bams_beds
        .join(GENERATE_CPG_TABLE.out.cpgtable)
        .combine(ch_num_turns)
        .map { meta1, bam, bai, bed, cpgtable, meta2, num_turns -> 
            def new_meta = [id: "${meta1.id}_${meta2.id}", dataset: meta1.id, simulate: meta2.id]
            [new_meta, bam, bai, bed, cpgtable, num_turns]
        }
        .set { ch_simulation_input }

    SIMULATE_SEQUENCES(
        ch_simulation_input, 
        ch_fasta.map { meta, fasta -> fasta }
    )

    METHYLDACKEL_EXTRACT_SIMULATED(
        SIMULATE_SEQUENCES.out.bam,
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai },
    )

    // TODO Compare the simulated and real methylation calls
    // METHYLDACKEL_EXTRACT.out.bedgraph
    //     .combine(METHYLDACKEL_EXTRACT_SIMULATED.out.bedgraph)
    //     .set { ch_bedgraphs }

    // EVALUATE(
    //     ch_bedgraphs
    // )
}
