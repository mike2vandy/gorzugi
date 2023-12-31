
import pandas as pd

indList = pd.read_csv("lists/all.list")
individs = indList['samples']

rule all:
  input:
    expand("indSFS/sfs/{ind}.sfs", ind = individs) 

rule indSAF:
  input:
    "../../bams/{ind}.bam"
  output:
    "indSFS/saf/{ind}.saf.idx"
  params:
    ref = "../../reference/PaintedTurtle.dna.toplevel.fa",
    scaffs = "../../reference/regions.txt",
    outpre = "indSFS/saf/{ind}"
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      angsd -i {input} \
        -ref {params.ref} \
        -anc {params.ref} \
        -out {params.outpre} \
        -rf {params.scaffs} \
        -GL 1 -minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 \
        -setMinDepthInd 2 -setMaxDepthInd 25 \
        -doCounts 1 -baq 2 -C 50 -doMajorMinor 1 -doSaf 1 -P {threads}
    '''
 
rule indSFS:
  input:
    saf = "indSFS/saf/{ind}.saf.idx"
  output:
    sfs = "indSFS/sfs/{ind}.sfs"
  params:
  conda: "env/angsd.yaml"
  threads: 4
  shell:
    '''
      realSFS {input.saf} -fold 1 -cores {threads} > {output.sfs}
    '''

 
