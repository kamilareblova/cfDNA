process P1FastqToSam {
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "broadinstitute/picard:3.4.0"
 
	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
        tuple val(name), val(sample), path("${name}.unaligned.bam")

	script:
	"""
	echo FastqToSam $name
        java -jar /usr/picard/picard.jar FastqToSam O=${name}.unaligned.bam F1=$fwd F2=$rev SM=${name} LB=Library1 PU=Unit1 PL=Illumina
	"""
}

process P2ExtractUMI {

        tag "ExtractUMI on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "xreblova/fgbio:1"

        input:
        tuple val(name), val(sample), path("${name}.unaligned.bam")

        output: 
        tuple val(name), val(sample), path("${name}.unaligned_bam_umi_extracted.bam")
 

        script: 
        """
        echo ExtractUMI $name

        fgbio ExtractUmisFromBam --input=${name}.unaligned.bam --output=${name}.unaligned_bam_umi_extracted.bam --read-structure=5M2S+T 5M2S+T --molecular-index-tags=ZA ZB --single-tag=RX
        """
}


process P3BAMtoFASTQ {

        tag "BAMtoFASTQ on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "broadinstitute/picard:3.4.0"

        label "l_cpu"
        label "l_mem"


        input:
        tuple val(name), val(sample), path("${name}.unaligned_bam_umi_extracted.bam")

        output:
        tuple val(name), val(sample), path("${name}.fastq"), path("${name}.unaligned_bam_umi_extracted.bam")

        script:
        """
        echo BAMtoFASTQ $name
        java -jar /usr/picard/picard.jar SamToFastq I=${name}.unaligned_bam_umi_extracted.bam F=${name}.fastq INTERLEAVE=true
        """
}

process P4ALIGN {
     
        tag "ALIGN1 on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        
        label "l_cpu"
        label "l_mem"

        input:
        tuple val(name), val(sample), path("${name}.fastq"), path("${name}.unaligned_bam_umi_extracted.bam")

        output:
        tuple val(name), val(sample), path("${name}.aligned_bam_umi_extracted.bam"), path("${name}.unaligned_bam_umi_extracted.bam") 

        script:
        """
        echo ALIGN1 $name
        source activate bwa

        bwa mem -p -t $task.cpus ${params.refindex}.fa ${name}.fastq | samtools sort -@ 8 -o ${name}.aligned_bam_umi_extracted.bam
        """
}

process P5MERGE {

        tag "5MERGE on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "broadinstitute/picard:3.4.0"

        label "l_cpu"
        label "l_mem"

        input: tuple val(name), val(sample), path("${name}.aligned_bam_umi_extracted.bam"), path("${name}.unaligned_bam_umi_extracted.bam")

        output:
        tuple val(name), val(sample), path("${name}.aligned_tag_umi.bam")

        script:
        """
        echo 5MERGE $name
        java -jar /usr/picard/picard.jar MergeBamAlignment UNMAPPED=${name}.unaligned_bam_umi_extracted.bam UNMAPPED=${name}.unaligned_bam_umi_extracted.bam ALIGNED=${name}.aligned_bam_umi_extracted.bam O=${name}.aligned_tag_umi.bam R=${params.refindex}.fa CLIP_ADAPTERS=false VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true EXPECTED_ORIENTATIONS=FR MAX_GAPS=-1 SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=false 
        """
}

process P6CollectHs {

        tag "6CollectHs on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "broadinstitute/picard:3.4.0"

        label "l_cpu"
        label "l_mem"

        input: 
        tuple val(name), val(sample), path("${name}.aligned_tag_umi.bam")

        output:
        tuple val(name), val(sample), ${name}.picard.txt

        script:
        """
        echo 6CollectHs $name
        java -jar /usr/picard/picard.jar CollectHsMetrics -I ${name}.aligned_tag_umi.bam -O ${name}.picard.txt -R ${params.refindex}.fa --BAIT_INTERVALS ${params.picardintervallist} --TARGET_INTERVALS ${params.varbed1} --PER_TARGET_COVERAGE ${params.varbed1}
        """
}
      
  
process P7GroupbyUMI {

        tag "7GroupbyUMI on $name using $task.cpus CPUs and $task.memory memory"
        publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
        container "xreblova/fgbio:1"

        input:
        tuple val(name), val(sample), path("${name}.aligned_tag_umi.bam")

        output:
        tuple val(name), val(sample), path("${name}.groupedbyumi.bam")


        script:
        """
        echo 7GroupbyUMI $name
        fgbio GroupReadsByUmi --strategy={paired} --input=${name}.aligned_tag_umi.bam --output=${name}.groupedbyumi.bam --raw-tag=RX --min-map-q=10 --edits=1
        """
}


process GATK {
       tag "GATK on $name"
       publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'
        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.vcf")

        script:
        """
        echo GATK $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.ref}.fa -I $bam -L ${params.varbed2}  --dont-use-soft-clipped-bases true -A StrandBiasBySample -minimum-mapping-quality 0 --mapping-quality-threshold-for-genotyping 0 --enable-dynamic-read-disqualification-for-genotyping true --flow-filter-alleles-qual-threshold 0 -O ${name}.vcf
        """
}


process VAFaNORMALIZACE {
        tag "VAFaNORMALIZACE on $name"
        publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path(gatk)

        output:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo VAFaNORMALIZACE $name

        bcftools +fill-tags $gatk -Ob -o ${name}.pom2.bcf -- -t FORMAT/VAF
        bcftools convert -O v -o ${name}.vaf.vcf ${name}.pom2.bcf 
         
        bcftools norm -f ${params.ref}.fa -m -both ${name}.vaf.vcf -o ${name}.norm.vcf
        bgzip ${name}.norm.vcf
        tabix ${name}.norm.vcf.gz

        """
}

process ANOTACE_ACGT {
        tag "ANOTACEACGT on $name"
        // publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.vcf.gz"), path("${name}.norm.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.acgt.vcf.gz"), path("${name}.norm.acgt.vcf.gz.tbi")

        script:
        """
        source activate gatk4610
        echo ANOTACEACGT $name
        gatk --java-options "-Xmx4g"  VariantAnnotator   -V ${name}.norm.vcf.gz -O ${name}.norm.acgt.vcf.gz --resource:ACGT ${params.ACGT} --expression ACGT.AF   --expression ACGT.AC   --expression ACGT.AC_Hom   --expression ACGT.AC_Het   --expression ACGT.AC_Hemi

        tabix ${name}.norm.acgt.vcf.gz
        """
}

process ANOTACE_annovar {
       tag "ANOTACE on $name"
       //publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.acgt.vcf.gz"), path("${name}.norm.acgt.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz.tbi")

        script:
        """
        source activate bcftoolsbgziptabix
        echo ANOTACE $name

        ${params.annovar} -vcfinput ${name}.norm.metarnn.vcf.gz ${params.annovardb}  -buildver hg38 -protocol refGeneWithVer,ensGene,1000g2015aug_all,1000g2015aug_eur,exac03nontcga,avsnp150,clinvar_20250721,dbnsfp41c,gnomad41_exome,gnomad41_genome,cosmic70,revel,GTEx_v8_eQTL \
        -operation gx,g,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo -polish -xreffile ${params.gene_fullxref.txt} -arg '-splicing 50 -exonicsplicing',,,,,,,,,,,, --remove
        bgzip ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf
        tabix ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz
        """
}

process VCF2TXT {
       tag "VCF2TXT on $name"
       // publishDir "${params.outDirectory}/${sample.run}/varianty/", mode:'copy'

        input:
        tuple val(name), val(sample), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz"), path("${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz.tbi")

        output:
        tuple val(name), val(sample), path("${name}.final.txt")

        script:
        """
        echo VCF2TXT $name
        source activate gatk4610
        gatk --java-options "-Xmx4g" VariantsToTable -R ${params.ref}.fa  --show-filtered  -V ${name}.norm.metarnn.vcf.gz.hg38_multianno.vcf.gz -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -GF SB -GF VAF -F dedicnostAR -F dedicnostAD -F dedicnostXlinked -F dedicnostYlinked -F fenotyp  -F ACGT.AF -F ACGT.AC -F ACGT.AC_Hom -F ACGT.AC_Het -F ACGT.AC_Hemi -F Func.refGeneWithVer -F Gene.refGeneWithVer -F GeneDetail.refGeneWithVer -F ExonicFunc.refGeneWithVer -F AAChange.refGeneWithVer -F 1000g2015aug_all -F 1000g2015aug_eur  -F gnomad41_exome_AF -F gnomad41_exome_AF_nfe -F gnomad41_genome_AF -F gnomad41_genome_AF_nfe -F avsnp150 -F CLNSIG -F REVEL -F MetaRNN.Varsome -F SIFT_pred -F MutationTaster_pred -F Gene_full_name.refGeneWithVer -F FATHMM_pred -F PROVEAN_pred -F Function_description.refGeneWithVer -F Disease_description.refGeneWithVer -F Tissue_specificityUniprot.refGeneWithVer -F Expression-egenetics.refGeneWithVer --output ${name}.final.txt
        """
}


process COVERAGE1 {
          tag "COVERAGE1 on $name"
       publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
         container "staphb/samtools:1.20"

        input:
        tuple val(name), val(sample), path(bam), path(bai)

        output:
        tuple val(name), val(sample), path("${name}.coveragefin.txt")

        script:
        """
        echo COVERAGE1 $name
        samtools bedcov ${params.varbed1} $bam -d 20 > ${name}.COV
        awk '{print \$5/(\$3-\$2)}'  ${name}.COV >  ${name}.COV-mean
        awk '{print (\$6/(\$3-\$2))*100"%"}' ${name}.COV > ${name}-procento-nad-20
        paste ${name}.COV-mean ${name}-procento-nad-20 > vysledek
        echo "chr" "start" "stop" "name" ${name}.COV-mean ${name}-procento-nad-20 > hlavicka
        sed -i 's/ /\t/'g hlavicka
        paste ${params.varbed1} vysledek > coverage
        cat hlavicka coverage > ${name}.coveragefin.txt
        sed -i -e "s/\r//g" ${name}.coveragefin.txt
        """
}

process COMBINECOVERAGEMEAN {
    tag "COMBINECOVERAGEMEAN"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/mapped/" }, mode: 'copy'

    output:
    path "coveragemeanALL"


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$5}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed1} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coveragemeanALL
    """
}

process COMBINECOVERAGEPROCENTA {
    tag "COMBINECOVERAGEPROCENTA"

    input:
    tuple val(run_name), path(coverage_files)

    publishDir { "${params.outDirectory}/${run_name}/mapped/" }, mode: 'copy'

    output:
    path "coverageprocentoALL"


    script:
    """
    i=0
    tmpfiles=""
    for file in ${coverage_files}; do
        awk '{print \$6}' "\$file" > column_\$i.txt
        tmpfiles="\$tmpfiles column_\$i.txt"
        i=\$((i+1))
    done
    echo "chr" "start" "stop" "name" > hlavicka
    sed -i 's/ /\t/'g hlavicka
    cat hlavicka ${params.varbed1} > bedshlavickou
    paste  bedshlavickou \$tmpfiles > coverageprocentoALL
    """
}

workflow {
        rawfastq = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    .splitCsv(header: true)
    .map { row ->
        def baseDir = new File("${params.baseDir}")
        def runDir = baseDir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.endsWith(row.run)
            }
        })[0] //get the real folderName that has prepended date

        def fileR1 = file("${runDir}/processed_fastq/${row.name}_R1.fastq.gz", checkIfExists: true)
        def fileR2 = file("${runDir}/processed_fastq/${row.name}_R2.fastq.gz", checkIfExists: true)

                def meta = [name: row.name, run: row.run]
        [
            meta.name,
            meta,
            fileR1,
            fileR2,
                ]
    }
     . view()

unalignedbam = P1FastqToSam(rawfastq)
//umiextracted = P2ExtractUMI(unalignedbam)
//fastqnew = P3BAMtoFASTQ(umiextracted)
//alignedbam = P4ALIGN(fastqnew)
//merged = P5MERGE(alignedbam)
//metrika = P6CollectHs(merged)


//aligned = ALIGN(---------)
//varcalling = GATK(aligned)
//normalizovany = VAFaNORMALIZACE(varcalling)
//anotovanyacgt = ANOTACE_ACGT(normalizovany)
//anotovany = ANOTACE_annovar(anotovanyacgt)
//anotovanyfin = VCF2TXT(anotovany)

//coverage_results = COVERAGE1(aligned)
//coverage_files_collected = coverage_results
//    .map { name, sample, f -> tuple(sample.run, file(f)) }
//    .groupTuple() // groups by sample.run automatically!
//finalcoverage = COMBINECOVERAGEMEAN(coverage_files_collected)
//finalprocenta = COMBINECOVERAGEPROCENTA(coverage_files_collected)


}
