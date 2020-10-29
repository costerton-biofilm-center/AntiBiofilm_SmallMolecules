import argparse
import os
import re
import shlex
import subprocess
import pdb

def parse_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-m", "--metadata",
                            help = "Sample file/metadata list. A tsv file where the sample R1 "
                            "full paths are the first column. Only need to list the R1 files",
                            required = True
                            )
    arg_parser.add_argument("-t", "--threads",
                        help = "Threads for running the analysis",
                        required = True
                        )
    arg_parser.add_argument("-r", "--reference", 
                        help = "Directory to bowtie2 ref genome.",
                        required=True)
    arg_parser.add_argument("-o", "--output_dir", 
                        help = "Base dir for output files.",
                        required=True)
    arg_parser.add_argument("-a", "--annotation", 
                    help = "Path to annotation file.",
                    required=True)

    return arg_parser.parse_args()


def run_sample(sample_info):
    sample_info['path_R1'] = sample_info['File_Path']
    sample_info['path_R2'] = re.sub("_R1", "_R2", sample_info['path_R1'])
    
    analysis = alignment_pipeline(**sample_info)

def alignment_pipeline(**sample_info):
    """ 1. Trim reads
        2. Align to reference genome
        3. Count reads
    """
    #Read trimming with cutadapt

    read1_adapt = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    read2_adapt = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

    trimmed_R1 = f"{sample_info['Sample_Name']}_R1.trimmed.fastq.gz"
    trimmed_R2 = f"{sample_info['Sample_Name']}_R2.trimmed.fastq.gz"


    cutadapt_cmd = f"cutadapt -j 0 -a {read1_adapt}" \
        f" -A {read2_adapt} -m 10 " \
        f"-o {sample_info['output_dir']}{trimmed_R1} " \
        f"-p {sample_info['output_dir']}{trimmed_R2} " \
        f"{sample_info['path_R1']} {sample_info['path_R2']}"
    
    trm = subprocess.Popen(
        shlex.split(cutadapt_cmd), 
        shell=False,
        stdout=open(
            f"{sample_info['output_dir']}{sample_info['Sample_Name']}_cutadaptout.log",'w')
        )  

    trm.communicate()

    ####

    sample_info['aligned_file'] = f"{sample_info['Sample_Name']}.aligned.bam"

    with open(f"{sample_info['output_dir']}{sample_info['aligned_file']}",'bw') as algn_file:

        bt2_cmd = f"bowtie2 -p {sample_info['threads']} " \
            f" -x {sample_info['reference']} " \
            f"-1 {sample_info['output_dir']}{trimmed_R1} " \
            f"-2 {sample_info['output_dir']}{trimmed_R2} "
        
        print(bt2_cmd)

        bt2 = subprocess.Popen(
            shlex.split(bt2_cmd),
            stdout = subprocess.PIPE,
            stderr = open(f"{sample_info['output_dir']}{sample_info['Sample_Name']}.bt2.log",'w')
        )

        samtools_cmd = f"samtools view -Sbh -@ {sample_info['threads']}"

        samtools = subprocess.Popen(
            shlex.split(samtools_cmd),
            stdin=bt2.stdout,
            stdout=algn_file
        )
        
        samtools.communicate()

    #Count reads with featureCounts

    fcount_out = f"{sample_info['Sample_Name']}.counts"

    print(fcount_out)

    fcount_cmd = f"featureCounts -p -a {sample_info['annotation']} " \
                f" -o {sample_info['output_dir']}{fcount_out} -t 'gene' -g 'gene_id' " \
                f" -T {sample_info['threads']}" \
                f" {sample_info['output_dir']}{sample_info['aligned_file']}"
    

    print(shlex.split(fcount_cmd))

    fcount = subprocess.Popen(
        shlex.split(fcount_cmd))
    
    fcount.communicate()


def main():
    cmd_args = parse_args()

    with open(cmd_args.metadata) as samples:
        header = samples.readline()
        header = header.rstrip("\n").split("\t")

        header = [x.strip(' ') for x in header]  # clean headers

        for line in samples: 
            cols = dict(zip(header, line.rstrip('\n').split('\t')))

            run_sample({**cols, 
            'threads' : cmd_args.threads, 
            'output_dir' : cmd_args.output_dir,
            'reference' : cmd_args.reference,
            'annotation' : cmd_args.annotation})

if __name__ == '__main__':
    main()
