# mcaller - joint genotyper
[ currently in beta version, I am preparing a quick version in C ]

This program is used to perform joint genotyping of multiple BAM files.
The output is a simple tab separated file (currently not VCF).

## Install 
These are the requirements (listed in requirements.txt)

    a) pysam           
    b) scipy
    c) numpy
    d) Biopython
 
 ### To install deps:
  
    a) Either install globally on machine: 
          pip install -r requirements
          
    b) or create a virtualenv:
          virtualenv venv
          source venv/bin/activate
          pip install -r requirements
          
 ## Usage
 Basic usage:
 
       python mcaller.py -g <genome_fasta> -b <bed_file> -i <BAM1>,<BAM2>,...,<BAMN> -o <output_folder>
       
 The first BAM file is treated as the germline file
 
 ### Options
 
      -h, --help            show this help message and exit
      -g GENOME_FILE, --ref=GENOME_FILE
                            path to fasta (required)
      -b BED, --bed=BED     bed file (required)
      -i BAM, --in=BAM      bam files, semicolon ';' separated (required)
      -o OUT_FOLDER, --out=OUT_FOLDER
                            output folder to store files (required)
      -y, --coord           treat bed as exact coordiantes (quicker if mutations
                            are known)
      -x MIN_ALIGNMENT_SCORE, --aln=MIN_ALIGNMENT_SCORE
                            minimum alignment score
      -s REGION_INTERVAL, --reg=REGION_INTERVAL
                            region size for noise calculation
      -n NOISE, --noise=NOISE
                            noise probability
      -t NORMAL_MUT_THRESH, --ncont=NORMAL_MUT_THRESH
                            normal sample max contamination
      -w MIN_COV, --mincov=MIN_COV
                            minimum coverage during analysis
      -e MAX_READ_EDIT_DISTANCE, --edit=MAX_READ_EDIT_DISTANCE
                            max edit distance
      -c FILTER_COVERAGE, --covfilt=FILTER_COVERAGE
                            filter mutations where coverage was not reached by
                            samples
      -p BACKGROUND_P_THRESH, --pval=BACKGROUND_P_THRESH
                            p-value cutoff
      -q FREQ_FILTER_HIGH_COV, --qfreq=FREQ_FILTER_HIGH_COV
                            mutation frequency filter for high coverage regions
      -f FREQ_THRESH, --freq=FREQ_THRESH
                            mutation frequency filter for non-high coverage
                            regions
      -u MIN_NORM_COV, --ncov=MIN_NORM_COV
                            minimum reads in normal sample
      -m MIN_SOMATIC_MUT, --mut=MIN_SOMATIC_MUT
                            minimum mutant reads in cancer sample
      -a MINIMUM_INDEL_MUT, --indel=MINIMUM_INDEL_MUT
                            minimum mutant reads supporting indel in cancer sample
      -d PROX_INDEL_MUT_COUNT, --maxprox=PROX_INDEL_MUT_COUNT
                            maximum proximal reads with indel before mutation
                            filtering
      -j PROX_INDEL_DIST, --indeldist=PROX_INDEL_DIST
                            distance when looking for proximal indels
      -k, --proxfile        (boolean) filter for proximal indels (0:no / 1:yes)
      -l INDEL_NOISE_FREQ, --indnoise=INDEL_NOISE_FREQ
                            default noise when calculating background
                            
 ## Output example
 
 ### SNP
      chr1	187153	T	C	somatic	1	24:0:0:0:24,1:1:1:1,1:1:1:1:1,1:1:1:1:1	34:0:4:0:30,1:0.018:1:1,0.81:1:1:1:1,1:0.78:1:1:1
Fields

    Chromosome: chr1
    position: 187153
    reference: T
    alteration: C
    mutation status: somatic
    filtering status: 1 (accept)
    germline stats: 24:0:0:0:24,1:1:1:1,1:1:1:1:1,1:1:1:1:1
    tumor mutation stats: 34:0:4:0:30,1:0.018:1:1,1:1:1:1:1,1:1:1:1:1

Explanation of stats (flag)

    stats forma (from tumor):
    coverages: 34:0:4:0:30 [ TOTAL:A:C:G:T ]
    probability: 1:0.018:1:1 [ A:C:G:T ]
    strand bias probability: 0.81:1:1:1:1 [ A:C:G:T ]
    quality bias probability: 1:0.78:1:1:1 [ A:C:G:T ]

Longer explanations:

        a) coverages: 34:0:4:0:30 [ TOTAL:A:C:G:T ]
              Total coverage: 34
              nucl. "A" cov: 0
              nucl. "C" cov: 4
              nucl. "G" cov: 0
              nucl. "T" cov: 30
              
        b) probability: 1:0.018:1:1 [ A:C:G:T ]
              prob. for "A" as mutation: 1
              prob. for "C" as mutation: 0.018 (significant)
              prob. for "G" as mutation: 1
              prob. for "T" as mutation: 1 (this is the ref base)
              
        c) strand bias probability: 1:1:1:1 [ A:C:G:T ]
              prob. for "A" strand bias: 1
              prob. for "C" strand bias: 1
              prob. for "G" strand bias: 1
              prob. for "T" strand bias: 1
              
        d) quality bias probability: 1:0.78:1:1 [ A:C:G:T ]
              prob. for "A" qual. bias: 1
              prob. for "C" qual. bias: 0.78
              prob. for "G" qual. bias: 1
              prob. for "T" qual. bias: 1
