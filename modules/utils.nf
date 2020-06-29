process performHostFilter {
    input:
        tuple(val(sampleId), path(forward), path(reverse))
    output:
        tuple sampleId, path("${sampleId}_hostfiltered_R1.fastq.gz"), path("${sampleId}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.composite_ref} ${forward} ${reverse} | \
            filter_non_human_reads.py -c ${params.viral_contig_name} > ${sampleId}.viral_and_nonmapping_reads.bam
        samtools sort -n ${sampleId}.viral_and_nonmapping_reads.bam | \
             samtools fastq -1 ${sampleId}_hostfiltered_R1.fastq.gz -2 ${sampleId}_hostfiltered_R2.fastq.gz -s ${sampleId}_singletons.fastq.gz -
        """
}
