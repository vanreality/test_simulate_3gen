include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main.nf'
include { SAMTOOLS_MERGE } from './modules/nf-core/samtools/merge/main.nf'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main.nf'
include { METHYLDACKEL_EXTRACT } from './modules/nf-core/methyldackel/extract/main.nf' 
include { METHYLDACKEL_EXTRACT as METHYLDACKEL_EXTRACT_SIMULATED } from './modules/nf-core/methyldackel/extract/main.nf'

process GENERATE_CPG_TABLE {
    tag "${meta.id}"

    input:
    tuple val(meta), path(bam), path(bed)
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
    tuple val(meta), path(bam), path(bed), val(num_turns)
    path(cpgtable)
    path(fasta)

    output:
    tuple val(meta), path('*_simulated.bam'), emit: bam

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
        .fromPath('samplesheet.csv')
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.label]
            [meta, file(row.bam)]
        }
        .groupTuple(by: 0)
        .set { ch_bams }

    Channel
        .fromPath('samplesheet.csv')
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

    SAMTOOLS_SORT(SAMTOOLS_MERGE.out.bam)

    METHYLDACKEL_EXTRACT(
        SAMTOOLS_SORT.out.bam,
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai },
    )

    SAMTOOLS_SORT.out.bam
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

    // Join with other channels
    ch_merged_bams_beds
        .combine(ch_num_turns)
        .map { meta1, bam, bed, meta2, num_turns -> 
            def new_meta = meta1 + meta2
            [new_meta, bam, bed, num_turns]
        }
        .set { ch_simulation_input }

    ch_simulation_input.view()

    SIMULATE_SEQUENCES(
        ch_simulation_input, 
        GENERATE_CPG_TABLE.out.cpgtable.map { meta, cpgtable -> cpgtable }, 
        ch_fasta.map { meta, fasta -> fasta }
    )

    METHYLDACKEL_EXTRACT_SIMULATED(
        SIMULATE_SEQUENCES.out.bam,
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai },
    )
}
