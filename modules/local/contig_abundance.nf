process CONTIG_ABUNDANCE {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

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

    total_all=\$(samtools view -c -F 0x900 "${bam}")

    awk -v TOTAL="\$total_all" 'BEGIN{
        OFS="\\t";
        print "Contig","Length","Reads","Percent_mapped","Percent_all_reads"
      }
      \$1!="*" {
        c=\$1; L[c]=\$2; R[c]=\$3; SUM+=\$3
      }
      END{
        for (c in R){
          pm = (SUM   ? 100*R[c]/SUM   : 0);
          pa = (TOTAL ? 100*R[c]/TOTAL : 0);
          print c, L[c], R[c], sprintf("%.3f%%", pm), sprintf("%.3f%%", pa)
        }
        unmapped = (TOTAL - SUM); if (unmapped < 0) unmapped = 0
        print "UNMAPPED", 0, unmapped, sprintf("%.3f%%", 0), sprintf("%.3f%%", (TOTAL ? 100*unmapped/TOTAL : 0))
      }' "${prefix}.idxstats.tsv" \
    | sort -k1,1 > "${prefix}.contig_abundance.tsv"

    samtools coverage -A filt.bam > "${prefix}.coverage_per_contig.tsv" || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | head -1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
