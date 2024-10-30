# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:06:39 2023

@author: rjovelin
"""


import os
import subprocess
import argparse
import gzip
import json
from itertools import zip_longest




def group_fastqs(fastqs):
    '''
    (list) -> list
    
    Returns a list of lists of paired fastqs
    
    parameters
    ----------
    - fastqs (list): List of paired fastqs files
    '''
    
    fastqs = sorted(list(fastqs))
    
    L = []
    for i in range(0, len(fastqs), 2):
        L.append([fastqs[i], fastqs[i+1]])
    return L
    


def getread(fastq_file):
    """
    (file) -- > itertools.zip_longest
    :param fastq_file: a fastq file open for reading in plain text mode
    
    Returns an iterator slicing the fastq into 4-line reads.
    Each element of the iterator is a tuple containing read information
    """
    args = [iter(fastq_file)] * 4
    return zip_longest(*args, fillvalue=None)



def downsize_fastqs(r1_file, r2_file, outdir, max_reads):
    """
    (str, str, str, int) -> None
       
    Write new fastqs pairs in outdir keeping the top max_reads in each new file 
    Pre-condition: fastqs are paired, have the same number of reads and files are in sync
        
    Parameters
    ----------
    - r1_file (str): Path to fastq R1
    - r2_file (str): Path to fastq R2
    - outdir (str): Output directory
    - max_reads (int): Maximum number of reads to keep from each fastq
    """
    
    # Read FASTQ in text mode
    r1 = gzip.open(r1_file, "rt")
    r2 = gzip.open(r2_file, "rt")
    
    # Open output fastqs in text mode for writing re-headered reads
    r1_writer = gzip.open(os.path.join(outdir, os.path.basename(r1_file)), "wt")
    r2_writer = gzip.open(os.path.join(outdir, os.path.basename(r2_file)), "wt")
    
    # make a list of fastqs open for reading
    fastqs = [r1, r2]
    
    # make a list of files open for writing
    writers = [r1_writer, r2_writer]

    # create iterator with reads from each file
    I =  zip(getread(fastqs[0]), getread(fastqs[1]))
    
    # count all reads
    total = 0
    
    # loop over iterator with slices of 4 read lines from each line
    for reads in I:
        # extract reads if the maximum number of reads is not reached
        if total <= max_reads:
            for i in range(len(writers)):
                # write new fastqs
                for j in range(len(reads[i])):
                    writers[i].write(reads[i][j])
            # update counter
            total += 1
                
    # close all open files
    for i in writers:
        i.close()
    for i in fastqs:
        i.close()
       
    print("Complete. Output written to {0}".format(outdir))
    
    


def compute_average_reads(fastqs):
    '''
    (dict) -> float
    
    Returns the average number of reads across all fastq files
    
    Parameters
    ----------
    - fastqs (dict): Dictionary of read counts across all fastq files
    '''

    total = 0
    for i in fastqs:
        total += fastqs[i]
    average = total / len(fastqs)
    return average


def compute_max_reads(average, downsize):
    '''
    (float, int) -> int
    
    Returns the maximum number of reads of each fastq files by reducing the 
    average number of reads by a factor of downsize %
    
    Parameters
    ----------
    - average (float): Average number of reads across all samples
    - downsize (int): Downsizing factor, in percent
    '''
       
    cap = average * downsize / 100
    return int(cap)



def link_files(targetdir, fastqs):
    '''
    (str, list) -> None
    
    Link the fastqs files to the target directory
    
    Parameters
    ----------
    - targetdir (str): Directory where the fastas are linked
    - fastqs (list): List of fastqs files
    '''

    for file in fastqs:
        assert os.path.isfile(file)
        filename = os.path.basename(file)
        link = os.path.join(targetdir, filename)
        if os.path.isfile(link) == False:
           os.symlink(file, link)



def extract_fastqs(project, run, fpr):
    '''
    (str, str, str) -> dict
    
    Returns a dictionary of file paths and associated read counts
    
    Parameters
    ----------
    - project (str): Name opf project of interest
    - run (str): Run identifier
    - fpr (str): Path to the File Provenance Reporter
    '''
    
    # open fpr, grab fastqs for the given run, and extract read counts
    counts = {}
    infile = gzip.open(fpr, 'rt', errors='ignore')
    for line in infile:
        if args.run in line:
            line = line.rstrip().split('\t')
            if line[1] == project or line[1] == 'GLCS':
                workflow = line[30]
                if workflow.lower() in ['casava', 'bcl2fastq', 'fileimportforanalysis', 'fileimport']:
                    # keep sequences of the correct run
                    if run == line[18]:
                        file_path = line[46]
                        read_count = line[45]
                        if read_count:
                            read_count = {k.split('=')[0]:k.split('=')[1] for k in read_count.split(';')}
                        if 'read_count' in read_count:
                            read_count = int(read_count['read_count'])
                        else:
                            read_count = -1
                        counts[file_path] = read_count
    infile.close()           

    return counts


def link_fastqs(args):
    
    
    # get working directory
    WorkingDir = os.path.join(args.workingdir, '{0}_{1}'.format(args.run, args.project))

    # exctract fastqs and their read counts 
    counts = extract_fastqs(args.project, args.run, args.fpr)
    
    fastqs = sorted(list(counts.keys()))
    print('extracted {0} fastqs from FPR'.format(len(fastqs)))

    # create target directory
    if args.downsize:
        # downsize the fastqs and keep a percent of the average number of reads
        # across samples
        try:
            int(args.downsize)
        except:
            raise ValueError('the downsize factor should be an integer')
        # compute the mean number of reads across samples
        average = compute_average_reads(counts)
        # compute the maximum number of reads
        max_reads = compute_max_reads(average, int(args.downsize))
        # group all the paired fastqs
        paired_fastqs = group_fastqs(fastqs)
        # downsize by keeping the top max reads from each file
        outdir = os.path.join(WorkingDir, 'fastqs')
        os.makedirs(outdir, exist_ok=True)
        for pair in paired_fastqs:
            r1_file, r2_file = pair
            downsize_fastqs(r1_file, r2_file, outdir, max_reads)
    else:
        # link the original fastqs 
        targetdir = os.path.join(WorkingDir, 'fastqs')
        os.makedirs(targetdir, exist_ok=True)
        link_files(targetdir, fastqs)

    # write json with read counts
    count_file = os.path.join(WorkingDir, '{0}.{1}.read_counts.json'.format(args.project, args.run))
    newfile = open(count_file, 'w')
    json.dump(counts, newfile)
    newfile.close()
    
    
    # write file with downsizing factor, max reads to include on the report
    if args.downsize:
        outputfile = os.path.join(WorkingDir, '{0}.{1}.downsize_summary.txt'.format(args.project, args.run))
        newfile = open(outputfile, 'w')
        summary = 'The mean number of reads across samples ({0}) has been reduced by {1}%.'.format(round(average, 3), args.downsize)
        summary = summary + ' ' + 'Only the top {0} reads have been kept for each sample\n'.format(max_reads)
        newfile.write(summary)
        newfile.close()
    
    
    
    
    

def generate_qsubs(args):
    
    # get working directory
    WorkingDir = os.path.join(args.workingdir, '{0}_{1}'.format(args.run, args.project))

    FastqDir = os.path.join(WorkingDir, 'fastqs')
    analysisdir = os.path.join(WorkingDir, 'smmips_analysis')
    qsubdir = os.path.join(WorkingDir, 'qsubs')
    donedir = os.path.join(qsubdir, 'done')
    logdir = os.path.join(qsubdir, 'logs')
    os.makedirs(analysisdir, exist_ok=True)
    os.makedirs(qsubdir, exist_ok=True)
    os.makedirs(donedir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    # make a list of fastqs to analyse
    Fastqs = [os.path.join(FastqDir, i) for i in os.listdir(FastqDir) if 'fastq' in i]    
    print('# fastq files:', len(Fastqs))

    # map fastqs to samples
    fastq_files = {}
    for i in Fastqs:
        sample = os.path.basename(i)
        sample = sample[:sample.index('_' + args.run)]
        if sample not in fastq_files:
            fastq_files[sample] = []
        fastq_files[sample].append(i)
    for i in fastq_files:
        fastq_files[i].sort()
        assert len(fastq_files[i]) == 2
        assert 'R1' in fastq_files[i][0] and 'R2' in fastq_files[i][1]       

    if os.path.isfile(args.panel) == False:
        raise ValueError("Invalid panel path. {0}".format(args.panel))
        
    bwa='/.mounts/labs/gsi/modulator/sw/Debian8.11/bwa-0.7.12/bin/bwa'
    reference='/.mounts/labs/gsi/modulator/sw/data/hg19-bwa-index-0.7.12/hg19_random.fa'

    # write qsubs for aligning fastqs
    align_cmd = 'module load smmips/1.0.9; smmips align -f1 {0} -f2 {1} -o {2} -r {3} -bwa {4} -pf {5}' 

    for sample in fastq_files:
        BashScript = os.path.join(qsubdir, sample + '.align.sh')
        newfile = open(BashScript, 'w')
        fastq1 = fastq_files[sample][0]
        fastq2 = fastq_files[sample][1]
        outdir = os.path.join(analysisdir, sample)
        mycmd = align_cmd.format(fastq1, fastq2, outdir, reference, bwa, sample)
        newfile.write(mycmd)
        newfile.close()
        # create a qsub script
        qsubscript = os.path.join(qsubdir, sample + '.align.qsub')
        newfile = open(qsubscript, 'w')
        newfile.write('qsub -P gsi -b y -N {0}.smmipQC.alignment.{1} -l h_vmem=30g,h_rt=5:0:0 -e {2} -o {2} \"bash {3}\"'.format(sample, args.run, logdir, BashScript))
        newfile.close()
    

    bed_file = os.path.join(WorkingDir, 'smmipRegions_500.bed')
    infile = open(bed_file)
    regions = infile.read().rstrip().split('\n')
    infile.close()
    for i in range(len(regions)):
        regions[i] = '.'.join(regions[i].split())
        
    # write qsubs for smmip assignments
    assign_cmd = 'module load smmips/1.0.9; smmips assign -pa {0} -o {1} -r {2} -b {3} -pf {4} -s 2 -up 5; mv {5} {6}'

    for sample in fastq_files:
        for region in regions:
            # create a shell script
            BashScript = os.path.join(qsubdir, sample + '.{0}.assign.sh'.format(region))
            qsubscript = os.path.join(qsubdir, sample + '.{0}.assign.qsub'.format(region))
            newfile = open(BashScript, 'w')
            fastq1 = fastq_files[sample][0]
            fastq2 = fastq_files[sample][1]
            outdir = os.path.join(analysisdir, sample)
            sortedbam = os.path.join(WorkingDir, 'smmips_analysis/{0}/out/{0}.sorted.bam'.format(sample))
            mycmd = assign_cmd.format(args.panel, outdir, region, sortedbam, sample, qsubscript, donedir) 
            newfile.write(mycmd)
            newfile.close()
            # create a qsub script
            newfile = open(qsubscript, 'w')
            newfile.write('qsub -P gsi -cwd -b y -N QC.{0}.assign.{1}.{2} -l h_vmem=30g,h_rt=5:0:0 -e {3} -o {3} \"bash {4}\"'.format(args.run, sample, region, logdir, BashScript))
            newfile.close()


    merge_cmd = 'module load smmips/1.0.9; smmips merge -pf {0} -ft {1} -t {2}'


    # make lists of stats files to merge for each sample
    for sample in fastq_files:
        outdir = os.path.join(analysisdir, sample)
        statsdir = os.path.join(outdir, 'stats')
        extraction = [os.path.join(statsdir, '{0}_temp.{1}.extraction_metrics.json'.format(sample, region)) for region in regions]
        counts = [os.path.join(statsdir, '{0}_temp.{1}.smmip_counts.json'.format(sample, region)) for region in regions]
        
        # write qsubs for merging
        BashScript = os.path.join(qsubdir, sample + '.merge.extraction.sh')
        newfile = open(BashScript, 'w')
        prefix = os.path.join(statsdir, sample)
        mycmd = merge_cmd.format(prefix, 'extraction', ' '.join(extraction))
        newfile.write(mycmd)
        newfile.close()
            
        # create a qsub script
        qsubscript = os.path.join(qsubdir, sample + '.merge.extraction.qsub')
        newfile = open(qsubscript, 'w')
        newfile.write('qsub -P gsi -cwd -b y -N {0}.smmipQC.merge.extraction.{1} -l h_vmem=20g,h_rt=5:0:0 -e {2} -o {2} \"bash {3}\"'.format(sample, args.run, logdir, BashScript))
        newfile.close()
        
        # write qsubs for merging
        BashScript = os.path.join(qsubdir, sample + '.merge.counts.sh')
        newfile = open(BashScript, 'w')
        prefix = os.path.join(statsdir, sample)
        mycmd = merge_cmd.format(prefix, 'counts', ' '.join(counts))
        newfile.write(mycmd)
        newfile.close()
            
        # create a qsub script
        qsubscript = os.path.join(qsubdir, sample + '.merge.counts.qsub')
        newfile = open(qsubscript, 'w')
        newfile.write('qsub -P gsi -cwd -b y -N {0}.smmipQC.merge.counts.{1} -l h_vmem=20g,h_rt=5:0:0 -e {2} -o {2} \"bash {3}\"'.format(sample, args.run, logdir, BashScript))
        newfile.close()


def run_alignments(args):
    
    # get working directory
    WorkingDir = os.path.join(args.workingdir, '{0}_{1}'.format(args.run, args.project))

    FastqDir = os.path.join(WorkingDir, 'fastqs')
    analysisdir = os.path.join(WorkingDir, 'smmips_analysis')
    qsubdir = os.path.join(WorkingDir, 'qsubs')
    donedir = os.path.join(qsubdir, 'done')
    logdir = os.path.join(qsubdir, 'logs')
    os.makedirs(analysisdir, exist_ok=True)
    os.makedirs(qsubdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    # make a list of a qsubs
    qsubs = [os.path.join(qsubdir, i) for i in os.listdir(qsubdir) if 'align.qsub' in i]

    while len(qsubs) != 0:
        running_jobs = int(subprocess.check_output('qstat | wc -l', shell=True).decode('utf-8').rstrip())
        #print('running: {0}'.format(running_jobs))
        #print('total jobs left: {0}'.format(len(qsubs)))        
        if running_jobs < int(args.max_jobs):
            jobs_to_submit = int(args.max_jobs) - running_jobs
            if jobs_to_submit > 0:
                #print('jobs to submit', jobs_to_submit)
                current_qsubs = qsubs[:jobs_to_submit]
                for i in current_qsubs:
                    job = subprocess.call('bash {0}'.format(i), shell=True)
                    qsubs.remove(i)
                    
                    
def run_qc(args):

    # get working directory
    WorkingDir = os.path.join(args.workingdir, '{0}_{1}'.format(args.run, args.project))

    FastqDir = os.path.join(WorkingDir, 'fastqs')
    analysisdir = os.path.join(WorkingDir, 'smmips_analysis')
    qsubdir = os.path.join(WorkingDir, 'qsubs')
    donedir = os.path.join(qsubdir, 'done')
    logdir = os.path.join(qsubdir, 'logs')
    os.makedirs(analysisdir, exist_ok=True)
    os.makedirs(qsubdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    # make a list of a qsubs
    qsubs = [os.path.join(qsubdir, i) for i in os.listdir(qsubdir) if 'assign.qsub' in i]

    # count the number of jobs
    
    while len(qsubs) != 0:
        running_jobs = int(subprocess.check_output('qstat | wc -l', shell=True).decode('utf-8').rstrip())
        #print('running: {0}'.format(running_jobs))
        #print('total jobs left: {0}'.format(len(qsubs)))        
        if running_jobs < int(args.max_jobs):
            jobs_to_submit = int(args.max_jobs) - running_jobs
            if jobs_to_submit > 0:
                #print('jobs to submit', jobs_to_submit)
                current_qsubs = qsubs[:jobs_to_submit]
                for i in current_qsubs:
                    job = subprocess.call('bash {0}'.format(i), shell=True)
                    qsubs.remove(i)
                    

def merge_metrics(args):
    
    
    # get working directory
    WorkingDir = os.path.join(args.workingdir, '{0}_{1}'.format(args.run, args.project))
    qsubdir = os.path.join(WorkingDir, 'qsubs')
    
    # make a list of qsubs
    qsubs1 = [os.path.join(qsubdir, i) for i in os.listdir(qsubdir) if '.merge.extraction.qsub' in i]
    qsubs2 = [os.path.join(qsubdir, i) for i in os.listdir(qsubdir) if '.merge.counts.qsub' in i]
    qsubs = qsubs1 + qsubs2
    
    for i in qsubs:
        subprocess.call('bash {0}'.format(i), shell=True)
    

def set_up_analysis(args):
    
    # get working directory
    WorkingDir = os.path.join(args.workingdir, '{0}_{1}'.format(args.run, args.project))
    
    FastqDir = os.path.join(WorkingDir, 'fastqs')
    analysisdir = os.path.join(WorkingDir, 'smmips_analysis')
    qsubdir = os.path.join(WorkingDir, 'qsubs')
    donedir = os.path.join(qsubdir, 'done')
    logdir = os.path.join(qsubdir, 'logs')
    os.makedirs(analysisdir, exist_ok=True)
    os.makedirs(qsubdir, exist_ok=True)
    os.makedirs(logdir, exist_ok=True)

    
    # create qsub script
    qsubscript = os.path.join(WorkingDir, 'Run_qc_analyses.qsub')
    qsubfile = open(qsubscript, 'w')

    # create qsub to link fastqs
    bashscript = os.path.join(qsubdir, '{0}_{1}_link_fastqs.sh'.format(args.project, args.run))
    newfile = open(bashscript, 'w')
    mycmd = 'module load smmip-qc; sleep 60; smmipqc link -fpr {0} -r {1} -pr {2} -w {3}'.format(args.fpr, args.run, args.project, args.workingdir) 
    newfile.write(mycmd)
    newfile.close()
    link_job_name = 'link.{0}.{1}'.format(args.project, args.run)
    qsubfile.write('qsub -P gsi -cwd -b y -N {0} -l h_rt=5:0:0 -e {1} -o {1} \"bash {2}\"'.format(link_job_name, logdir, bashscript))
    qsubfile.write('\n')    
    
    # create qsub to generate bed file
    bashscript = os.path.join(qsubdir, 'Create_regions_bed.sh')
    newfile = open(bashscript, 'w')
    
    if os.path.isfile(args.panel) == False:
        raise ValueError("Invalid panel path. {0}".format(args.panel))
    
    newfile.write('module load smmip-region-finder; smmipRegionFinder -d 500 -p {0} -o {1}'.format(args.panel, WorkingDir))
    newfile.close()
    bed_job_name = 'regions.{0}.{1}'.format(args.project, args.run)
    qsubfile.write("qsub -P gsi -cwd -b y -N {0}  -hold_jid \"{1}\" -l h_vmem=20g,h_rt=5:0:0 -e {2} -o {2} \"bash {3}\"".format(bed_job_name, link_job_name, logdir, bashscript))
    qsubfile.write('\n')
    
    # create qsub to generate alignment, assignment and merge qsubs
    bashscript = os.path.join(qsubdir, 'Create_qsubs.sh')
    newfile = open(bashscript, 'w')
    newfile.write('module load smmip-qc; sleep 60; smmipqc qsubs -r {0} -pr {1} -w {2} -p {3}'.format(args.run, args.project, args.workingdir, args.panel)) 
    newfile.close()
    qsub_job_name= 'script.{0}.{1}'.format(args.project, args.run)
    qsubfile.write("qsub -P gsi -cwd -b y -N {0} -hold_jid \"{1}\" -l h_vmem=20g,h_rt=5:0:0 -e {2} -o {2} \"bash {3}\"".format(qsub_job_name, bed_job_name, logdir, bashscript))
    qsubfile.write('\n')

    # create qsub to run alignments
    bashscript = os.path.join(qsubdir, 'Run_alignments.sh')
    align_job_name = 'align.{0}.{1}'.format(args.project, args.run) 
    newfile = open(bashscript, 'w')
    newfile.write('module load smmip-qc; sleep 60; smmipqc align -r {0} -pr {1} -w {2} -m {3}'.format(args.run, args.project, args.workingdir, args.max_jobs)) 
    newfile.close()
    qsubfile.write("qsub -P gsi -cwd -b y -N {0} -hold_jid \"{1}\" -l h_vmem=20g,h_rt={2}:0:0 -e {3} -o {3} \"bash {4}\"".format(align_job_name, qsub_job_name, args.align_run_time, logdir, bashscript))
    qsubfile.write('\n')
    
    # create qsub for smmip assignment
    bashscript = os.path.join(qsubdir, 'Run_qc.sh')
    qc_job_name = 'qc.{0}.{1}'.format(args.project, args.run) 
    newfile = open(bashscript, 'w')
    newfile.write('module load smmip-qc; sleep 60; smmipqc qc -r {0} -pr {1} -w {2} -m {3}'.format(args.run, args.project, args.workingdir, args.max_jobs)) 
    newfile.close()
    qsubfile.write("qsub -P gsi -cwd -b y -N {0} -hold_jid \"{1}\" -l h_vmem=20g,h_rt={2}:0:0 -e {3} -o {3} \"bash {4}\"".format(qc_job_name, align_job_name, args.qc_run_time, logdir, bashscript))
    qsubfile.write('\n')    
    
    # create qsub for merging
    bashscript = os.path.join(qsubdir, 'Merge_metrics.sh')
    merge_job_name = 'merge.{0}.{1}'.format(args.project, args.run) 
    newfile = open(bashscript, 'w')
    newfile.write('module load smmip-qc; sleep 60; smmipqc merge -r {0} -pr {1} -w {2}'.format(args.run, args.project, args.workingdir)) 
    newfile.close()
    qsubfile.write("qsub -P gsi -cwd -b y -N {0} -hold_jid \"{1}\" -l h_vmem=20g,h_rt=5:0:0 -e {2} -o {2} \"bash {3}\"".format(merge_job_name, qc_job_name, logdir, bashscript))
    qsubfile.write('\n')
        
    qsubfile.close()
    
if __name__ == '__main__':

    # create top-level parser
    parser = argparse.ArgumentParser(prog = 'smmipqc.py', description='A tool to set up manual analyses of smmip QC')
    subparsers = parser.add_subparsers(help='sub-command help', dest='subparser_name')
    
   	# link files in gsi space 
    l_parser = subparsers.add_parser('link', help="Link fastq files to workfing directory")
    l_parser.add_argument('-fpr', '--provenance', dest='fpr', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    l_parser.add_argument('-r', '--run', dest='run', help='Run id', required=True)
    l_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report.', required=True)
    l_parser.add_argument('-w', '--workingdir', dest='workingdir', help='Path to working directory', required=True)
    l_parser.add_argument('-d', '--downsize', dest='downsize', help='Reduce the number of reads by a percent of the average number of reads across samples')
    l_parser.set_defaults(func=link_fastqs)
    
    # generate qsubs 
    q_parser = subparsers.add_parser('qsubs', help="Generate qsubs for alignment, assignment and merging")
    q_parser.add_argument('-r', '--run', dest='run', help='Run id', required=True)
    q_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report.', required=True)
    q_parser.add_argument('-w', '--workingdir', dest='workingdir', help='Path to working directory', required=True)
    q_parser.add_argument('-p', '--panel', dest='panel', help='Smmip panel. Panels are:\
                          /.mounts/labs/gsiprojects/gsi/smmips/smmips_panels/panels/panels_MISO/Myeloid_smMIPS_panel.txt,\
                          /.mounts/labs/gsiprojects/gsi/smmips/smmips_panels/panels/panels_MISO/Myeloid_CML_smMIPS_panel.txt,\
                          /.mounts/labs/gsiprojects/gsi/smmips/smmips_panels/panels/panels_MISO/CB_Paired_B_ALL_smMIPS_panel.txt', required = True)
    q_parser.set_defaults(func=generate_qsubs)


    # run alignments 
    a_parser = subparsers.add_parser('align', help="Launch alignment qsubs")
    a_parser.add_argument('-r', '--run', dest='run', help='Run id', required=True)
    a_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report.', required=True)
    a_parser.add_argument('-w', '--workingdir', dest='workingdir', help='Path to working directory', required=True)
    a_parser.add_argument('-m', '--maxjobs', dest='max_jobs', default=150, help='Maximum number of jobs running in parallel')
    a_parser.set_defaults(func=run_alignments)

    
    # run smmip assignments 
    qc_parser = subparsers.add_parser('qc', help="Launch assignment qsubs")
    qc_parser.add_argument('-r', '--run', dest='run', help='Run id', required=True)
    qc_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report.', required=True)
    qc_parser.add_argument('-w', '--workingdir', dest='workingdir', help='Path to working directory', required=True)
    qc_parser.add_argument('-m', '--maxjobs', dest='max_jobs', default=150, help='Maximum number of jobs running in parallel')
    qc_parser.set_defaults(func=run_qc)

    # merge metrics 
    m_parser = subparsers.add_parser('merge', help="Merge extraction and count metrics")
    m_parser.add_argument('-r', '--run', dest='run', help='Run id', required=True)
    m_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report.', required=True)
    m_parser.add_argument('-w', '--workingdir', dest='workingdir', help='Path to working directory', required=True)
    m_parser.set_defaults(func=merge_metrics)
 
   
   	# set analyses 
    s_parser = subparsers.add_parser('analysis', help="Set qsubs to run smmip QC analyses")
    s_parser.add_argument('-fpr', '--provenance', dest='fpr', default='/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz', help='Path to File Provenance Report. Default is /scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz')
    s_parser.add_argument('-r', '--run', dest='run', help='Run id', required=True)
    s_parser.add_argument('-pr', '--project', dest='project', help='Project name as it appears in File Provenance Report.', required=True)
    s_parser.add_argument('-w', '--workingdir', dest='workingdir', help='Path to working directory', required=True)
    s_parser.add_argument('-p', '--panel', dest='panel', help='Smmip panel. Panels are:\
                          /.mounts/labs/gsiprojects/gsi/smmips/smmips_panels/panels/panels_MISO/Myeloid_smMIPS_panel.txt,\
                          /.mounts/labs/gsiprojects/gsi/smmips/smmips_panels/panels/panels_MISO/Myeloid_CML_smMIPS_panel.txt,\
                          /.mounts/labs/gsiprojects/gsi/smmips/smmips_panels/panels/panels_MISO/CB_Paired_B_ALL_smMIPS_panel.txt', required = True)
    s_parser.add_argument('-m', '--maxjobs', dest='max_jobs', default=150, help='Maximum number of jobs running in parallel')
    s_parser.add_argument('-art', '--align_run_time', dest='align_run_time', default=5, help='Run time in hours for aligning reads. Default is 5 hours')
    s_parser.add_argument('-qrt', '--qc_run_time', dest='qc_run_time', default=10, help='Run time in hours for performing QC. Default is 10 hours')
    s_parser.add_argument('-d', '--downsize', dest='downsize', help='Reduce the number of reads by a percent of the average number of reads across samples')
    s_parser.set_defaults(func=set_up_analysis)
    
    # get arguments from the command line
    args = parser.parse_args()
    # pass the args to the default function
    args.func(args)

