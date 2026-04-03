process CONTIG_ABUNDANCE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.idxstats.tsv")           , emit: idxstats
    tuple val(meta), path("*.contig_abundance.tsv")    , emit: abundance
    tuple val(meta), path("*.coverage_per_contig.tsv") , emit: coverage
    tuple val(meta), path("*.flagstat.txt")            , emit: flagstat
    path "versions.yml"                                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set -euo pipefail

    samtools flagstat "${bam}" > "${prefix}.flagstat.txt"

    samtools view -b -F 0x904 -q 10 "${bam}" > filt.bam
    samtools index filt.bam

    samtools idxstats filt.bam > "${prefix}.idxstats.tsv"

    total_primary=\$(samtools view -c -F 0x900 "${bam}")
    mapped_primary=\$(samtools view -c -F 0x904 "${bam}")
    mapped_q10=\$(awk 'BEGIN{sum=0} \$1!="*" {sum+=\$3} END{print sum+0}' "${prefix}.idxstats.tsv")

    lowmapq_primary=\$((mapped_primary - mapped_q10))
    unmapped_primary=\$((total_primary - mapped_primary))
    if (( lowmapq_primary < 0 )); then lowmapq_primary=0; fi
    if (( unmapped_primary < 0 )); then unmapped_primary=0; fi

    {
      printf "Contig\\tLength\\tReads_q10\\tPercent_mapped_q10\\tPercent_primary_reads\\n"
      awk -v TOTAL="\$total_primary" -v SUM="\$mapped_q10" 'BEGIN{OFS="\\t"}
        \$1!="*" {
          pm = (SUM   ? 100*\$3/SUM   : 0);
          pp = (TOTAL ? 100*\$3/TOTAL : 0);
          print \$1, \$2, \$3, sprintf("%.3f%%", pm), sprintf("%.3f%%", pp)
        }' "${prefix}.idxstats.tsv"
      awk -v TOTAL="\$total_primary" -v COUNT="\$lowmapq_primary" 'BEGIN{
          OFS="\\t";
          printf "LOWMAPQ_PRIMARY\\t0\\t%d\\tNA\\t%.3f%%\\n", COUNT, (TOTAL ? 100*COUNT/TOTAL : 0)
        }'
      awk -v TOTAL="\$total_primary" -v COUNT="\$unmapped_primary" 'BEGIN{
          OFS="\\t";
          printf "UNMAPPED_PRIMARY\\t0\\t%d\\tNA\\t%.3f%%\\n", COUNT, (TOTAL ? 100*COUNT/TOTAL : 0)
        }'
    } > "${prefix}.contig_abundance.tsv"

    samtools coverage -A filt.bam > "${prefix}.coverage_per_contig.tsv" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
