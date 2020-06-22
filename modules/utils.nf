process performHostFilter {
    input:
        tuple(val(sampleId), path(forward), path(reverse))
    output:
        tuple sampleId, path("${sampleId}_hostfiltered_R1.fastq.gz"), path("${sampleId}_hostfiltered_R2.fastq.gz"), emit: fastqPairs

    script:
        """
        bwa mem -t ${task.cpus} ${params.human_ref} ${forward} ${reverse} | samtools view -q 30 -U ${sampleId}.host_unaligned.bam -Sb > ${sampleId}.host_aligned.bam
        samtools sort -n ${sampleId}.host_unaligned.bam | \
             samtools fastq -1 ${sampleId}_hostfiltered_R1.fastq.gz -2 ${sampleId}_hostfiltered_R2.fastq.gz -s ${sampleId}_singletons.fastq.gz -
        """
}
