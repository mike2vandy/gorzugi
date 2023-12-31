
rule all:
  input:
    "allStr/ngsdist/p_distances.tab",
    "allStr/Ne/all.geno.gz", 
    "allStr/pruned.beagle.gz",
    expand("allStr/ngs_admixes/{K}_{iteration}.log", K = list(range(1,6)), iteration = list(range(1,11))),
    "allStr/pcangsd/all.cov"

rule popGLF:
  input:
    "bamlists/all.bamlist"
  output:
    "allStr/all.beagle.gz"
  params:
    ref = "../../reference/PaintedTurtle.dna.toplevel.fa",
    scaffs = "../../reference/regions.txt",
    outpre = "allStr/all"
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      ind=$(awk 'END {{print int(NR * 0.6)}}' {input})
  
      angsd -bam {input} \
        -ref {params.ref} \
        -out {params.outpre} \
        -rf {params.scaffs} \
        -GL 1 -minMapQ 30 -minQ 20 -skipTriallelic 1 \
        -minInd $ind -setMinDepthInd 2 -setMaxDepthInd 25 \
        -doCounts 1 -doMajorMinor 1 -doGLF 2 -doMaf 1 -minMaf 0.01 -SNP_pval 1e-6 -P {threads}
    '''

rule measureLD:
  input:
    geno = "allStr/all.beagle.gz"
  output:
    "allStr/all.ld"
  params:
    pre = "allStr/",
    ind = 141 
  threads: 8
  shell:
    '''
      zcat {input.geno} |grep -v marker |cut -f 1 |tr '_' '\t' > {params.pre}/ngsld.pos
      nsites=$(zcat {input.geno} |awk 'END {{print NR - 1}}') 

      ~/software/ngsLD/ngsLD \
        --geno {input.geno} \
        --probs \
        --n_ind {params.ind} \
        --n_sites $nsites \
        --pos {params.pre}/ngsld.pos \
        --max_kb_dist 100 \
        --min_maf 0.01 \
        --n_threads {threads} \
        --out {output}
    '''

rule prune:
  input: 
    sites = "allStr/all.ld"
  output:
    unlinked = 'allStr/pruned.sites'
  conda: "env/graphtools.yaml"
  threads: 1
  shell:
    '''
      ~/software/ngsLD/scripts/prune_ngsLD.py \
        --input {input.sites} \
        --max_dist 25000 --min_weight 0.2 \
        --output {output}
    '''

rule pullPruned:
  input:
    sites = 'allStr/pruned.sites',
    beagle = "allStr/all.beagle.gz" 
  output:
    prunedBeagle = "allStr/pruned.beagle.gz"
  threads: 1
  shell:
    '''
      scripts/pullPruned.py {input.sites} {input.beagle} |gzip > {output.prunedBeagle}
    '''

rule NGSadmix:
  input:
    beagle = "allStr/pruned.beagle.gz"
  output:
    "allStr/ngs_admixes/{K}_{iteration}.log"
  params:
    outdir = "allStr/ngs_admixes",
    k = "{K}",
    iterat = "{iteration}"
  conda: "env/angsd.yaml"
  threads: 2
  shell:
    '''
      NGSadmix -likes {input.beagle} \
        -K {params.k} \
        -P {threads} \
        -minMaf 0.01 \
        -o {params.outdir}/{params.k}_{params.iterat}

    '''

rule PCangsd:
  input:
    beagle = "allStr/pruned.beagle.gz"
  output:
    "allStr/pcangsd/all.cov"
  params:
    prefix = "allStr/pcangsd/all"
  conda: "env/pcangsd.yaml"
  threads: 2
  shell:
    '''
      pcangsd --beagle {input.beagle} \
        --n_eig 2 \
        --threads {threads} \
        --out {params.prefix} \
        --maf 0.01 
    '''

rule ngsDist:
  input:
    beagle = "allStr/pruned.beagle.gz"
  output:
    "allStr/ngsdist/p_distances.tab"
  params:
    names = "allStr/sampleNames.txt"
  threads: 4
  shell:
    '''
      nsites=$(zcat {input.beagle} |awk 'END {{print NR - 1}}')
 
      ~/software/ngsDist/ngsDist \
        --geno {input.beagle} \
        --probs \
        --labels {params.names} \
        --n_ind 141 \
        --evol_model 0 \
        --indep_geno \
        --n_sites $nsites \
        --n_threads {threads} \
        --out {output}
    '''

rule genoForNe:
  input:
    "bamlists/all.bamlist"
  output:
    "allStr/Ne/all.geno.gz"
  params:
    sites = "allStr/sites/unprunedSitesForNe.pos",
    pre = "allStr/Ne/all",
    scaffs = "../../reference/regions.txt" 
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      angsd -bam {input} \
        -doGeno 4 \
        -rf {params.scaffs} \
        -sites {params.sites} \
        -doPost 1 -GL 1 -doMajorMinor 1 -doMaf 1 \
        -out {params.pre} -P 8
    '''
