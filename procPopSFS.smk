
import pandas as pd

main = pd.read_csv("2Dpops.txt", sep = '\t').set_index(['pop1','pop2'], drop = False)
popLists = ['br','ps','pn','dr','rg','pecos','riogrande', "filNew", "filOld"]

rule all:
  input:
    expand("popSFS/{pops}/{pops}.win.pestPG", pops = popLists),
    expand("popSFS/2D/{units.pop1}_{units.pop2}/{units.pop1}_{units.pop2}.fst.idx", units = main.itertuples())

rule popSAF:
  input:
    "bamlists/{pops}.bamlist"
  output:
    "popSFS/{pops}/{pops}.saf.idx"
  params:
    ref = "../../reference/PaintedTurtle.dna.toplevel.fa",
    scaffs = "../../reference/regions.txt",
    outpre = "popSFS/{pops}/{pops}"
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      ind=$(awk 'END {{print int(NR * 0.6)}}' {input})
  
      angsd -bam {input} \
        -ref {params.ref} \
        -anc {params.ref} \
        -out {params.outpre} \
        -rf {params.scaffs} \
        -GL 1 -minMapQ 30 -minQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 \
        -minInd $ind -setMinDepthInd 2 -setMaxDepthInd 25 \
        -doCounts 1 -baq 2 -C 50 -doMajorMinor 1 -doSaf 1 -P {threads}
    '''

rule popSFS:
  input:
    saf = "popSFS/{pops}/{pops}.saf.idx"
  output:
    sfs = "popSFS/{pops}/{pops}.sfs"
  params:
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      realSFS {input.saf} -fold 1 -cores {threads} > {output.sfs}
    '''

rule thetas:
  input:
    saf = "popSFS/{pops}/{pops}.saf.idx",
    sfs = "popSFS/{pops}/{pops}.sfs"
  output:
    thetas = "popSFS/{pops}/{pops}.thetas.idx"
  params:
    outpre = "popSFS/{pops}/{pops}"
  conda: "env/angsd.yaml"
  threads: 2
  shell:
    '''
      realSFS saf2theta {input.saf} -sfs {input.sfs} -outname {params.outpre} -cores {threads} -fold 1
    '''

rule windowTheta:
  input:
    "popSFS/{pops}/{pops}.thetas.idx"
  output: 
    "popSFS/{pops}/{pops}.win.pestPG"
  params:
    step = 10000,
    outpre = "popSFS/{pops}/{pops}.win" 
  conda: "env/angsd.yaml"
  threads: 1
  shell:
    '''
      thetaStat do_stat {input} -win {params.step} -step {params.step} -outnames {params.outpre}
    '''

rule makeTwoDSFS:
  input:
    saf1 = "popSFS/{pop1}/{pop1}.saf.idx",
    saf2 = "popSFS/{pop2}/{pop2}.saf.idx"
  output:
    sfs = "popSFS/2D/{pop1}_{pop2}/{pop1}_{pop2}.sfs"
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      realSFS {input.saf1} {input.saf2} -fold 1 -cores {threads} > {output.sfs}
    '''

rule fst:
  input:
    saf1 = "popSFS/{pop1}/{pop1}.saf.idx",
    saf2 = "popSFS/{pop2}/{pop2}.saf.idx",
    sfs = "popSFS/2D/{pop1}_{pop2}/{pop1}_{pop2}.sfs"
  output:
    "popSFS/2D/{pop1}_{pop2}/{pop1}_{pop2}.fst.idx"
  params:
    outpre = "popSFS/2D/{pop1}_{pop2}/{pop1}_{pop2}"
  conda: "env/angsd.yaml"
  threads: 8
  shell:
    '''
      realSFS fst index {input.saf1} {input.saf2} \
        -sfs {input.sfs} \
        -fold 1 \
        -fstOut {params.outpre} \
        -whichFst 1 \
        -cores {threads}

      realSFS fst stats {output} > {params.outpre}.global.fst.txt
    '''

