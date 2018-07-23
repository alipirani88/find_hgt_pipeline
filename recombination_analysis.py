from __future__ import division
import os
import argparse
import re
import subprocess
import errno
import readline
from joblib import Parallel, delayed
import multiprocessing
import glob
import csv
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict
from logging_subprocess import *
from log_modules import *
import re
from Bio import SeqIO
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='Recombination/HGT Analysis.\nThe pipeline takes a list of fasta files and aligns All-vs-All using Nucmer.\nExtracts aligned region by parsing nucmer coordinate and gff/bed annotations to extract regions that matches user defined percent identity and minimum aligned length parameters.\nGenerates a preliminary reference database out of extracted aligned regions by deduplicating and removing containments.\nRemoves containment fragments from preliminary database using nucmer.\nPerforms nucmer alignment between pseudomolecule fasta file and final containment removed aligned fragments to generate an alignment matrix score.\n', formatter_class=RawTextHelpFormatter)
required = parser.add_argument_group('Required arguments')
required.add_argument('-filename', action='store', dest="filename", help='This file should contain a list of fasta filenames(one per line) that the user wants to use from argument -dir folder. For Genome coordinate consistency, make sure the fasta files are in a pseudomolecule format', required=True)
required.add_argument('-out', action='store', dest="out", help='Output directory to save the results', required=True)
required.add_argument('-prokka_dir', action='store', dest="prokka_dir", help='Directory containing results of Prokka annotation pipeline or individual sample folders consisting gff and bed file. The folder name should match the fasta file prefix.', required=True)
required.add_argument('-jobrun', action='store', dest="jobrun",
                    help='Type of job to run. Run script on a compute cluster, parallelly on local or on local system(default): cluster, parallel-local, local')
required.add_argument('-dir', action='store', dest="dir", help='Directory containing fasta files specified in -filename list', required=True)
required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique Analysis Name to save results with this prefix', required=True)
optional = parser.add_argument_group('Optional arguments')
#optional.add_argument('-filename_db', action='store', dest="filename_db", help='filenames of nucmer db', required=False)
#optional.add_argument('-dir_db', action='store', dest="dir_db", help='directory of nucmer db fasta files', required=False)
#optional.add_argument('-matrix', action='store', dest="matrix", help='Matrix to parse and remove containments', required=False)
optional.add_argument('-remove_temp', action='store', dest="remove_temp", help='Remove Temporary directories from /tmp/ folder: yes/no', required=False)
optional.add_argument('-steps', action='store', dest="steps",
                    help='Analysis Steps to be performed. Use All or 1,2,3,4,5 to run all steps of pipeline.'
                         '\n1: Align all assembly fasta input file against each other using Nucmer.'
                         '\n2: Parses the Nucmer generated aligned coordinates files, extract individual aligned fragments and their respective annotation for metadata.'
                         '\n3: Generate a database of these extracted aligned regions by deduplicating and removing containments using BBmaps dedupe tool.'
                         '\n4: Remove containments from preliminary database by running nucmer'
                         '\n5: Performs nucmer alignment between input fasta file and final containment removed extracted fragments to generate an alignment score matrix.', required=True)
optional.add_argument('-pbs', action='store', dest="pbs", help='Provide PBS memory resources for individual nucmer jobs. Default: nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00', required=False)

args = parser.parse_args()

# Pipeline Methods
def create_job(jobrun, commands_list, logger):
    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "cluster":
        """
        Supports only PBS clusters for now.
        """
        keep_logging('Running jobs in cluster mode', 'Running jobs in cluster mode', logger, 'info')
        job_directory = args.out + "/"
        make_sure_path_exists(job_directory)
        count = 0
        for i in commands_list:
            job_name = "nucmer_job_command_" + str(count)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, i)
            job_file_name = "%s/%s.pbs" % (job_directory, job_name)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            count += 1
            f1.close()
        pbs_dir = job_directory + "/*/*_parallel_cluster.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: bash %s' % i, 'Running: bash %s' % i, logger, 'info')
            #os.system("bash %s" % i)
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        keep_logging('\nRunning jobs in parallel-local mode', '\nRunning jobs in parallel-local mode', logger, 'info')
        command_file = commands_list
        #print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)


    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            #os.system(i)
            call("%s" % i, logger)

def create_job_parse(jobrun, commands_list, logger):

    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "cluster":
        """
        Supports only PBS clusters for now.
        """
        keep_logging('\nRunning jobs in cluster mode', '\nRunning jobs in cluster mode', logger, 'info')
        job_directory = args.out + "/" + "temp_jobs"
        make_sure_path_exists(job_directory)
        count = 0
        for i in commands_list:
            job_name = "nucmer_coordinates_job_command_" + str(count)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, i)
            job_file_name = "%s/%s.pbs" % (job_directory, job_name)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            count += 1
            f1.close()
        pbs_dir = job_directory + "/nucmer_coordinates_job_command_*.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            #os.system("bash %s" % i)
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        keep_logging('\nRunning jobs in parallel-local mode', '\nRunning jobs in parallel-local mode', logger, 'info')
        command_file = commands_list
        #print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)

    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        keep_logging('\nRunning jobs in local mode', '\nRunning jobs in local mode', logger, 'info')
        command_file = commands_list
        for i in command_file:
            #os.system(i)
            call("%s" % i, logger)

def create_job_containments(jobrun, commands_list, logger):

    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "cluster":
        """
        Supports only PBS clusters for now.
        """
        job_directory = dedup_database_dir + "/" + "temp_jobs"
        make_sure_path_exists(job_directory)
        count = 0
        for i in commands_list:
            job_name = "nucmer_job_command_" + str(count)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, i)
            job_file_name = "%s/%s.pbs" % (job_directory, job_name)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            count += 1
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = job_directory + "/nucmer_job_command_*.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        keep_logging('\nRunning jobs in parallel-local mode', '\nRunning jobs in parallel-local mode', logger, 'info')
        command_file = commands_list
        #print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)


    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            #os.system(i)
            call("%s" % i, logger)

def create_job_containments_final(jobrun, commands_list, logger):

    """
    Based on type of jobrun; generate jobs and run accordingly.
    :param jobrun:
    :param vcf_filenames:
    :return:
    """
    if jobrun == "cluster":
        """
        Supports only PBS clusters for now.
        """
        job_directory = final_containment_removed_dir + "/" + "temp_jobs"
        make_sure_path_exists(job_directory)
        count = 0
        for i in commands_list:
            job_name = "nucmer_job_command_" + str(count)
            job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=1,pmem=4000mb,walltime=6:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, i)
            job_file_name = "%s/%s.pbs" % (job_directory, job_name)
            f1=open(job_file_name, 'w+')
            f1.write(job_print_string)
            count += 1
            f1.close()
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        pbs_dir = job_directory + "/nucmer_job_command_*.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
            #os.system("qsub %s" % i)
            call("qsub %s" % i, logger)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        keep_logging('\nRunning jobs in parallel-local mode', '\nRunning jobs in parallel-local mode', logger, 'info')
        command_file = commands_list
        #print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)


    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            #os.system(i)
            call("%s" % i, logger)

def make_sure_path_exists(out_path):
    """
    Make sure the output folder exists or create at given path
    :param out_path:
    :return:
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name\n', 'Errors in output folder path! please change the output path or analysis name\n', logger, 'info')
            exit()

# Deprecated: Removed before release
# def run_nucmer(temp_cmd):
#     f = open(temp_cmd, 'w+')
#     command_array = []
#     for file in filenames_array:
#         filebase = os.path.basename(file)
#         mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
#         os.system(mkdir_command)
#         for file_again in filenames_array:
#             file_again_base = os.path.basename(file_again)
#             prefix = filebase.replace('.fasta', '') + "_" + file_again_base.replace('.fasta', '')
#             command = "nucmer --maxmatch --prefix=%s %s %s &>/dev/null && show-coords -r -c -l -T %s.delta > %s.coords && rm %s.delta && mv %s.coords %s/%s\n" % (prefix, file, file_again, prefix, prefix, prefix, prefix, args.out, filebase.replace('.fasta', ''))
#             command_array.append(command)
#             f.write(command)
#     f.close()
#     return command_array

def run_nucmer_parallel(temp_cmd, logger):
    command_array_cluster = []
    command_array_local = []

    for file in filenames_array:
        file_command_array = []
        filebase = os.path.basename(file)

        # Generate Output directory for each fasta file
        mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))

        #Set up file names
        filebase_command_cluster = "%s/%s/%s_commands_cluster.sh" % (args.out, filebase.replace('.fasta', ''), filebase)
        filebase_command_local = "%s/%s/%s_commands_local.sh" % (args.out, filebase.replace('.fasta', ''), filebase)
        run_parallel_script_cluster = "%s/%s/%s_parallel_cluster.pbs" % (args.out, filebase.replace('.fasta', ''), filebase)
        run_parallel_script_local = "%s/%s/%s_parallel_local.pbs" % (args.out, filebase.replace('.fasta', ''), filebase)



        #os.system(mkdir_command)
        call("%s" % mkdir_command, logger)

        for file_again in filenames_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + "_" + file_again_base.replace('.fasta', '')
            command_local = "nucmer --maxmatch --prefix=/tmp/%s %s %s &>/dev/null && show-coords -r -c -l -T /tmp/%s.delta > %s/%s/%s.coords && rm /tmp/%s.delta\n" % (prefix, file, file_again, prefix, args.out, filebase.replace('.fasta', ''), prefix, prefix)

            command_cluster = "nucmer --maxmatch --prefix=/tmp/$PBS_JOBID/%s %s %s &>/dev/null && show-coords -r -c -l -T /tmp/$PBS_JOBID/%s.delta > %s/%s/%s.coords && rm /tmp/$PBS_JOBID/%s.delta\n" % (
            prefix, file, file_again, prefix, args.out, filebase.replace('.fasta', ''), prefix, prefix)

            command_array_cluster.append(command_cluster)
            command_array_local.append(command_local)


            #f.write(command)

        f_command = open(filebase_command_local, 'w+')
        for cmd in command_array_local:
            f_command.write(cmd)
        f_command.close()

        f_command = open(filebase_command_cluster, 'w+')
        for cmd in command_array_cluster:
            f_command.write(cmd)
        f_command.close()

        run_parallel_script_file_local = open(run_parallel_script_local, 'w+')
        run_parallel_script_file_cluster = open(run_parallel_script_cluster, 'w+')

        job_name = filebase + "_parallel"

        ## Generate tmp directories
        pbs_tmp_cluster = "/tmp/$PBS_JOBID"
        pbs_tmp_local = "/tmp/"

        #runjobs = "cd %s\n~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/run_commands_parallel.py -command %s" % (args.out, filebase_command)
        runjobs_local = "cd %s\n~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/run_commands_parallel.py -command %s" % (pbs_tmp_local, filebase_command_local)
        runjobs_cluster = "mkdir %s\ncd %s\n~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/run_commands_parallel.py -command %s" % (pbs_tmp_cluster, pbs_tmp_cluster, filebase_command_cluster)

        job_print_string_local = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=8,pmem=4000mb,walltime=240:00:00\n#PBS -q flux\n#PBS -A esnitkin_flux\n#PBS -l qos=flux\n\n%s\n" % (job_name, runjobs_local)
        job_print_string_cluster = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=8,pmem=4000mb,walltime=240:00:00\n#PBS -q flux\n#PBS -A esnitkin_flux\n#PBS -l qos=flux\n\n%s\n" % (job_name, runjobs_cluster)

        run_parallel_script_file_local.write(job_print_string_local)
        run_parallel_script_file_cluster.write(job_print_string_cluster)

        pbs_scripts_local = glob.glob("%s/%s/%s_parallel_local.pbs" % (args.out, filebase.replace('.fasta', ''), filebase))
        pbs_scripts_cluster = glob.glob("%s/%s/%s_parallel_cluster.pbs" % (args.out, filebase.replace('.fasta', ''), filebase))


    keep_logging('Nucmer command list for each file can be found in: %s/*/*.sh' % args.out,
                 'Nucmer command list for each file can be found in: %s/*/*.sh' % args.out, logger, 'info')

    keep_logging('Nucmer jobs for each file can be found in: %s/*/*.pbs' % args.out,
                 'Nucmer jobs for each file can be found in: %s/*/*.pbs' % args.out, logger,
                 'info')


    return command_array_local, command_array_cluster

def run_nucmer_parallel_containment(temp_cmd, logger):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_containments_array:
        file_command_array = []
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (dedup_database_dir, filebase.replace('.fasta', ''))
        filebase_command = "%s/%s/%s_commands.sh" % (dedup_database_dir, filebase.replace('.fasta', ''), filebase)
        run_parallel_script = "%s/%s/%s_parallel.pbs" % (dedup_database_dir, filebase.replace('.fasta', ''), filebase)



        #os.system(mkdir_command)
        call("%s" % mkdir_command, logger)
        for file_again in filenames_containments_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + "_" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=/tmp/%s %s %s &>/dev/null && show-coords -r -c -l -T /tmp/%s.delta > %s/%s/%s.coords && rm /tmp/%s.delta\n" % (
            prefix, file, file_again, prefix, dedup_database_dir, filebase.replace('.fasta', ''), prefix, prefix)

            command_array.append(command)
            file_command_array.append(command)
            f.write(command)
        f_command = open(filebase_command, 'w+')
        for cmd in file_command_array:
            f_command.write(cmd)
        run_parallel_script_file = open(run_parallel_script, 'w+')
        job_name = filebase + "_parallel"
        runjobs = "cd %s\n~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/run_commands_parallel.py -command %s" % (dedup_database_dir, filebase_command)
        job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=8,pmem=4000mb,walltime=240:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, runjobs)
        run_parallel_script_file.write(job_print_string)
        pbs_scripts = glob.glob("%s/%s/%s_parallel.pbs" % (dedup_database_dir, filebase.replace('.fasta', ''), filebase))
        # for i in pbs_scripts:
        #     keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
    f.close()
    keep_logging('Nucmer containment command list to run in parallel-local mode can be found in: %s/*/*.sh' % dedup_database_dir,
                 'Nucmer containment command list to run paralleley can be found in: %s/*/*.sh' % dedup_database_dir, logger,
                 'info')
    keep_logging('Nucmer containment jobs to run in parallel-local mode can be found in: %s/*/*.pbs' % dedup_database_dir,
                 'Nucmer containment jobs to run in parallel-local mode can be found in: %s/*/*.pbs' % dedup_database_dir,
                 logger,
                 'info')
    return command_array

def run_nucmer_db(temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_array:
        file_command_array = []
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
        filebase_command = "%s/%s/%s_commands.sh" % (args.out, filebase.replace('.fasta', ''), filebase)
        run_parallel_script = "%s/%s/%s_parallel.pbs" % (args.out, filebase.replace('.fasta', ''), filebase)

        keep_logging('Nucmer containment commands to run in parallel-local mode can be found in: %s' % filebase_command, 'Nucmer containment commands to run in parallel-local mode can be found in: %s' % filebase_command, logger, 'info')
        keep_logging('Nucmer containment jobs to run in parallel-local mode can be found in: %s' % run_parallel_script, 'Nucmer containment jobs to run in parallel-local mode can be found in: %s' % run_parallel_script, logger, 'info')
        #os.system(mkdir_command)
        call("%s" % mkdir_command, logger)
        for file_again in filenames_db_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + ":" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=/tmp/%s %s %s &>/dev/null && show-coords -r -c -l -T /tmp/%s.delta > %s/%s/%s.coords && rm /tmp/%s.delta\n" % (prefix, file, file_again, prefix, args.out, filebase.replace('.fasta', ''), prefix, prefix)
            command_array.append(command)
            file_command_array.append(command)
            f.write(command)
        f_command = open(filebase_command, 'w+')
        for cmd in file_command_array:
            f_command.write(cmd)
        run_parallel_script_file = open(run_parallel_script, 'w+')
        job_name = filebase + "_parallel"
        runjobs = "cd %s\n~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/run_commands_parallel.py -command %s" % (args.out, filebase_command)
        job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=8,pmem=4000mb,walltime=240:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, runjobs)
        run_parallel_script_file.write(job_print_string)


    f.close()
    return command_array

def run_nucmer_db_final(temp_cmd, logger):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_array:
        file_command_array = []
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (final_containment_removed_dir, filebase.replace('.fasta', ''))
        filebase_command = "%s/%s/%s_commands.sh" % (final_containment_removed_dir, filebase.replace('.fasta', ''), filebase)
        run_parallel_script = "%s/%s/%s_parallel.pbs" % (final_containment_removed_dir, filebase.replace('.fasta', ''), filebase)
        keep_logging('Nucmer containment commands to run in parallel-local mode can be found in: %s' % filebase_command, 'Nucmer containment commands to run in parallel-local mode can be found in: %s' % filebase_command, logger, 'info')
        keep_logging('Nucmer containment job to run in parallel-local mode can be found in: %s' % run_parallel_script, 'Nucmer containment commands to run in parallel-local mode can be found in: %s' % run_parallel_script, logger, 'info')
        #os.system(mkdir_command)
        call("%s" % mkdir_command, logger)
        for file_again in containment_removed_db_filenames:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + ":" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=/tmp/%s %s %s &>/dev/null && show-coords -r -c -l -T /tmp/%s.delta > %s/%s/%s.coords && rm /tmp/%s.delta\n" % (
            prefix, file, file_again, prefix, final_containment_removed_dir, filebase.replace('.fasta', ''), prefix, prefix)

            command_array.append(command)
            file_command_array.append(command)
            f.write(command)
        f_command = open(filebase_command, 'w+')
        for cmd in file_command_array:
            f_command.write(cmd)
        run_parallel_script_file = open(run_parallel_script, 'w+')
        job_name = filebase + "_parallel"
        runjobs = "cd %s\n~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/run_commands_parallel.py -command %s" % (final_containment_removed_dir, filebase_command)
        job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=8,pmem=4000mb,walltime=240:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\n%s\n" % (job_name, runjobs)
        run_parallel_script_file.write(job_print_string)


    f.close()
    return command_array

def run_command(i):
    #os.system(i)
    call("%s" % i, logger)
    done = "done"
    return done

def generate_parse_coord_aggregate_jobs(jobrun, filenames_array, temp_cmd, logger):
    f = open(temp_cmd, 'w+')
    command_array = []
    job_directory = args.out + "/" + "temp_jobs"
    make_sure_path_exists(job_directory)
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python %s/run_nucmer_coordinates.py -folder %s -out %s -dir %s -prokka_dir %s\n" % (os.path.dirname(os.path.abspath(__file__)), folder, args.out, args.dir, args.prokka_dir)
        #print cmd
        command_array.append(cmd)
        f.write(cmd)
    f.close()
    #print len(command_array)
    return command_array

# Deprecated: Removed before release
# def generate_parse_coord_db_aggregate_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd, logger):
#     f = open(temp_cmd, 'w+')
#     command_array = []
#     job_directory = args.out + "/" + "temp_jobs"
#     make_sure_path_exists(job_directory)
#     for folder in filenames_array:
#         cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/run_nucmer_coordinates_db.py -folder %s -out %s -dir %s -prokka_dir %s\n" % (folder, args.out, args.dir, args.prokka_dir)
#         #print cmd
#         command_array.append(cmd)
#         f.write(cmd)
#         #job_name = os.path.basename(folder)
#         #job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.dir, cmd)
#         #job_file_name = "%s/%s_run_nucmer_coordinates.pbs" % (job_directory, job_name)
#         #f1=open(job_file_name, 'w+')
#         #f1.write(job_print_string)
#     #f1.close()
#     f.close()
#     #print len(command_array)
#     return command_array

def generate_parse_containments_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd, logger):
    f = open(temp_cmd, 'w+')
    command_array = []
    job_directory = dedup_database_dir + "/" + "temp_jobs"
    make_sure_path_exists(job_directory)
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/run_nucmer_coordinates_containments.py -folder %s -out %s -dir %s\n" % (folder, dedup_database_dir, dedup_database_dir)
        #print cmd
        command_array.append(cmd)
        f.write(cmd)
        #job_name = os.path.basename(folder)
        #job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.dir, cmd)
        #job_file_name = "%s/%s_run_nucmer_coordinates.pbs" % (job_directory, job_name)
        #f1=open(job_file_name, 'w+')
        #f1.write(job_print_string)
    #f1.close()
    f.close()
    #print len(command_array)
    return command_array

def generate_parse_containments_final_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd, logger):
    f = open(temp_cmd, 'w+')
    command_array = []
    job_directory = final_containment_removed_dir + "/" + "temp_jobs"
    make_sure_path_exists(job_directory)
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/run_nucmer_coordinates_containments.py -folder %s -out %s -dir %s\n" % (folder, final_containment_removed_dir, final_containment_removed_dir)
        #print cmd
        command_array.append(cmd)
        f.write(cmd)
        #job_name = os.path.basename(folder)
        #job_print_string = "#PBS -N %s\n#PBS -M apirani@med.umich.edu\n#PBS -m a\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=76:00:00\n#PBS -q fluxod\n#PBS -A esnitkin_fluxod\n#PBS -l qos=flux\n\ncd %s\n%s\n" % (job_name, args.dir, cmd)
        #job_file_name = "%s/%s_run_nucmer_coordinates.pbs" % (job_directory, job_name)
        #f1=open(job_file_name, 'w+')
        #f1.write(job_print_string)
    #f1.close()
    f.close()
    #print len(command_array)
    return command_array

#keep_logging('', '', logger, 'info')

#Main Steps: Start of pipeline
if __name__ == '__main__':

    """Start Timer"""
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now().strftime('%Y-%m-%d')

    # Generate logger object for logging subprocesses
    global logger
    logger = generate_logger(args.out, args.analysis_name, start_time_2)

    keep_logging('The Script started at: %s\n' % start_time, 'The Script started at: %s\n' % start_time, logger, 'info')

    keep_logging('\nThis pipeline will run the following steps:\n\n1: Align all fasta input file against each other using Nucmer.\n2: Parses the Nucmer generated aligned coordinates files, extract individual aligned fragments and their respective annotation for metadata.\n3: Generate a database of these extracted aligned regions by deduplicating and removing containments.\n4: Performs nucmer alignment between input fasta file and extracted regions to generate a matrix of alignment score.\n', '\nThis pipeline will run the following steps:\n\n1: Align all fasta input file against each other using Nucmer.\n2: Parses the Nucmer generated coordinates files and generate annotated results of aligned fragments.\n3: Generate a preliminary database of these extracted aligned regions by deduplicating and removing containments using BBmaps dedupe tool.\n4: Remove containments from preliminary database by running nucmer\n5: Performs nucmer alignment between input fasta file and final containment removed extracted fragments to generate an alignment score matrix.\n', logger, 'info')

    # GENERATE INPUT ASSEMBLY FASTA FILENAMES ARRAY
    filenames_array = []

    # Copy fasta files to /tmp/ directory.
    temp_fasta_dir = "/tmp/fasta_file"
    make_sure_path_exists(temp_fasta_dir)
    copy_fasta_tmp_cmd = "cp %s/* %s/" % (args.dir, temp_fasta_dir)
    #os.system(copy_fasta_tmp_cmd)
    call("%s" % copy_fasta_tmp_cmd, logger)

    with open(args.filename) as fp:
        for line in fp:
            line = line.strip()
            #line = args.dir + "/" + line
            line = temp_fasta_dir + "/" + line
            filenames_array.append(line)

    # CHECK FOR STEPS ARGUMENTS WITH COMMAND-LINE ARGUMENT STEPS; Run individual pipeline steps
    if args.steps:
        # split values provided with -steps argument and decide the starting point of pipeline
        if args.steps == "All":
            steps_string = "1,2,3,4,5"
            steps_list = steps_string.split(',')
        else:
            steps_list = args.steps.split(',')


        keep_logging('Running pipeline with steps argument: %s\n' % steps_list,
                     'Running pipeline with steps argument: %s\n' % steps_list, logger,
                     'info')

        if "1" in steps_list:
            # Run nucmer on all the input fasta pseudomolecule assembly fasta file. compare All-vs-All.
            keep_logging('Running Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n', 'Running Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n', logger,
                         'info')
            make_sure_path_exists(args.out)
            temp_cmd = "%s/temp_commands" % args.out
            command_array_local, command_array_cluster = run_nucmer_parallel(temp_cmd, logger)
            if args.jobrun:
                if args.jobrun == "parallel-local":
                    create_job(args.jobrun, command_array_local, logger)
                elif args.jobrun == "cluster":
                    create_job(args.jobrun, command_array_cluster, logger)

            keep_logging('End Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n',
                         'End Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n', logger,
                         'info')

        if "2" in steps_list:
            # Step 2
            keep_logging('Running Step 2: Parsing Nucmer Coordinate results generated by Nucmer Step 1 and extracting all the aligned regions.',
                         'Running Step 2: Parsing Nucmer Coordinate results generated by Nucmer Step 1 and extracting all the aligned regions.', logger,
                         'info')
            temp_cmd_parse_coord = "%s/temp_commands_parse_coord" % args.out
            command_array = generate_parse_coord_aggregate_jobs(args.jobrun, filenames_array, temp_cmd_parse_coord, logger)
            if args.jobrun:
                create_job_parse(args.jobrun, command_array, logger)

            # Move Aggregate score results. No use of these results. Remove when done. Pending...
            aggregate_results = args.out + "/aggregate_results/"
            make_sure_path_exists(aggregate_results)
            #os.system("mv %s/*.aligned %s/*.score %s" % (args.out, args.out, aggregate_results))
            call("mv %s/*.aligned %s/*.score %s" % (args.out, args.out, aggregate_results), logger)
            keep_logging(
                'Running: mv %s/*.aligned %s/*.score %s\n' % (args.out, args.out, aggregate_results),
                'Running: mv %s/*.aligned %s/*.score %s\n' % (args.out, args.out, aggregate_results),
                logger,
                'debug')
            keep_logging(
                'End Step 2: Parsing Nucmer Coordinate results generated by Nucmer Step 1 and extracting all the aligned regions.',
                'End Step 2: Parsing Nucmer Coordinate results generated by Nucmer Step 1 and extracting all the aligned regions.',
                logger,
                'info')

            # Add Nursing Home steps here later. The below commented block belongs to Nursing Home project. It should be used to organize the Species-pairs results.

        # if "3" in steps_list:
        #     print "The Results output directory containing species-pair folders is %s" % args.out
        #     list_command = "ls -d %s/*_* | rev | cut -d'/' -f1 | rev" % args.out
        #     species_pair_uniq_array = []
        #     species_array = []
        #     #print list_command
        #     list_of_species_pai_dir = []
        #     list_of_species_pai_dir = subprocess.check_output(list_command, shell=True)
        #     for i in list_of_species_pai_dir.split('\n'):
        #         species_name = (i.split('_'))[0]
        #         if species_name not in species_array and species_name != "":
        #             species_array.append(species_name)
        #     os.system("mkdir %s/redundant_species_pairs" % args.out)
        #     f1=open("%s/move_redundant.sh" % args.out, 'w+')
        #     for file in species_array:
        #         for file_again in species_array:
        #             if file != file_again:
        #                 generate_string = str(file) + "_" + str(file_again)
        #                 reverse_string = str(file_again) + "_" + str(file)
        #                 if reverse_string not in species_pair_uniq_array:
        #                     species_pair_uniq_array.append(generate_string)
        #                     move = "mv %s/%s %s/redundant_species_pairs/" % (args.out, reverse_string, args.out)
        #                     os.system(move)
        #                     f1.write(move + '\n')
        #     f1.close()

        if "3" in steps_list:
            keep_logging('Running Step 3: Generate a preliminary database from the extracted aligned region.',
                         'Running Step 3: Generate a preliminary database from the extracted aligned region.', logger,
                         'info')

            # Set up database directory path to be used for saving dedupe results
            database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"

            # Run dedupe to deduplicate sequences in a preliminary extracted aligned region database.
            dedupe_cmd = "/nfs/esnitkin/bin_group/bbmap/dedupe.sh in=%s/Extracted_aligned_region.fasta out=%s/Extracted_aligned_region_dedup.fasta" % (
            database_directory, database_directory)
            keep_logging('Running dedupe command: %s\n' % dedupe_cmd,
                         'Running dedupe command: %s\n' % dedupe_cmd, logger,
                         'debug')
            #os.system(dedupe_cmd)
            call("%s" % dedupe_cmd, logger)

            # Run dedupe to remove containment sequences from deduplicated preliminary extracted aligned region database results.
            dedupe_containment_cmd = "/nfs/esnitkin/bin_group/bbmap/dedupe.sh in=%s/Extracted_aligned_region_dedup.fasta out=%s/Extracted_aligned_region_dedup_cluster_99.fasta minidentity=99" % (
            database_directory, database_directory)
            #os.system(dedupe_containment_cmd)
            call("%s" % dedupe_containment_cmd, logger)

            keep_logging('Running dedupe containment command: %s\n' % dedupe_containment_cmd,
                         'Running dedupe containment command: %s\n' % dedupe_containment_cmd, logger,
                         'debug')
            # Prepare or extract metadata annotations for sequences that remained after dedupe deduplicate and containment removal steps
            tab = "\\t"
            get_meta_annotations = "for i in `grep '>' %s/Extracted_aligned_region_dedup_cluster_99.fasta | sed 's/>//g'`; do echo \"grep '^$i' %s/Extracted_aligned_region.fasta_meta.tsv | awk -F'\\t' -v OFS='%s' '{print \$1,\$2,\$3,\$4,\$5,\$6}'\"; done > %s/grep_meta.sh" % (
                database_directory, database_directory, tab, database_directory)

            keep_logging('Extracting Metadata annotations for deduplicated and containment removed sequences: %s\n' % get_meta_annotations,
                         'Extracting Metadata annotations for deduplicated and containment removed sequences: %s\n' % get_meta_annotations, logger,
                         'debug')
            keep_logging(
                'Extracting Metadata annotations for deduplicated and containment removed sequences: bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv\n' % (database_directory, database_directory),
                'Extracting Metadata annotations for deduplicated and containment removed sequences: bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv\n' % (database_directory, database_directory),
                logger,
                'debug')

            #os.system(get_meta_annotations)
            #os.system("bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv" % (database_directory, database_directory))

            call("%s" % get_meta_annotations, logger)

            call("bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv" % (database_directory, database_directory), logger)

            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
            make_sure_path_exists(dedup_database_dir)

            get_deduped_extracted_aligned_region_database_files = "for i in `awk -F'\\t' '{print $1}' %s/Extracted_aligned_region_dedup_cluster_99.tsv`; do  echo \"$i\" > %s/$i.txt; done" % (
            database_directory, dedup_database_dir)

            keep_logging(
                'Extract sequence names from deduplicated and containment removed sequences: %s\n' % get_deduped_extracted_aligned_region_database_files,
                'Extract sequence names from deduplicated and containment removed sequences:: %s\n' % get_deduped_extracted_aligned_region_database_files,
                logger,
                'debug')

            #os.system(get_deduped_extracted_aligned_region_database_files)
            call("%s" % get_deduped_extracted_aligned_region_database_files, logger)

            extract_regions_from_file = "for i in `ls %s/*.txt`; do base=`echo $i | sed 's/.txt//g'`; ~/bin/seqtk/seqtk subseq %s/Extracted_aligned_region_dedup_cluster_99.fasta $i > $base.fasta; done" % (
            dedup_database_dir, database_directory)

            keep_logging(
                'Extract individual sequences from deduplicated and containment removed sequence mfa file: %s\n' % extract_regions_from_file,
                'Extract individual sequences from deduplicated and containment removed sequence mfa file: %s\n' % extract_regions_from_file,
                logger,
                'info')
            #os.system(extract_regions_from_file)
            call(extract_regions_from_file, logger)

        if "4" in steps_list:
            keep_logging('Running Step 4: Remove containment sequences from the preliminary database fragments. Part one: Run Nucmer all deduped extracted sequences against each other.',
                         'Running Step 4: Remove containment sequences from the preliminary database fragments. Part one: Run Nucmer all deduped extracted sequences against each other.', logger,
                         'info')

            # Containment sequence: a sequence which is exactly/99% similar to another sequence or a sequence which is contained inside a larger sequence
            # To remove these containment sequences, run nucmer on extracted fragments vs extracted fragments
            database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
            os.chdir(dedup_database_dir)

            get_extracted_aligned_region_dedup_cluster_99_fasta_files = "ls *.fasta | grep -v 'containment' | grep -v 'Extracted' > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % (
                dedup_database_dir)

            keep_logging(
                'Save deduped extracted aligned individual fasta filenames: ls *.fasta > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt\n' % dedup_database_dir,
                'Save deduped extracted aligned individual fasta filenames: ls *.fasta > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt\n' % dedup_database_dir,
                logger,
                'info')
            #os.system(get_extracted_aligned_region_dedup_cluster_99_fasta_files)
            call("%s" % get_extracted_aligned_region_dedup_cluster_99_fasta_files, logger)

            # Run Step 1 nucmer method again but with only Extracted_aligned_region_dedup_cluster_99_fasta_filenames
            output_folder = dedup_database_dir

            # GENERATE FASTA FILENAMES ARRAY
            # Copy fasta files to /tmp/ directory.
            temp_fasta_dir_containment = "/tmp/containment_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment)
            copy_fasta_tmp_cmd = "cp %s/*.fasta %s/" % (dedup_database_dir, temp_fasta_dir_containment)
            #os.system(copy_fasta_tmp_cmd)
            call("%s" % copy_fasta_tmp_cmd, logger)


            filenames_containments_array = []
            with open("%s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % dedup_database_dir) as fp:
                for line in fp:
                    line = line.strip()
                    #line = dedup_database_dir + "/" + line
                    line = temp_fasta_dir_containment + "/" + line
                    filenames_containments_array.append(line)
            fp.close()
            temp_cmd = "%s/temp_commands" % output_folder

            # Run Nucmer All-vs-All for deduped extracted alignment fasta files against each other.
            command_array = run_nucmer_parallel_containment(temp_cmd, logger)
            if args.jobrun:
                create_job_containments(args.jobrun, command_array, logger)

        if "5" in steps_list:
            keep_logging(
                'Running Step 5: Remove containment sequences from the preliminary database fragments. Part two. Parse Nucmer coordinates file of deduped extracted aligned region vs each other file and find contained sequences.\n',
                'Running Step 5: Remove containment sequences from the preliminary database fragments. Part two. Parse Nucmer coordinates file of deduped extracted aligned region vs each other file and find contained sequences.\n',
                logger,
                'info')
            temp_cmd = "%s/temp_commands_parse_containments.sh" % args.out

            # Definition of Containment sequence: a sequence which is exactly/99% similar to another sequence or a smaller sequence which is contained inside a larger sequence.
            # To remove these containment sequences, run nucmer on extracted fragments vs  extracted fragments
            # Parse the nucmer generated coordinates file, generate a matrix of alignment ratios, find out...

            # Set up deduped database directories
            database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"

            keep_logging(
                'Setting up dedupe database directories to store final deduplicated and containment removed fragments: %s\n' % dedup_database_dir,
                'Setting up dedupe database directories to store final deduplicated and containment removed fragments: %s\n' % dedup_database_dir,
                logger,
                'info')
            os.chdir(dedup_database_dir)


            # GENERATE FASTA FILENAMES ARRAY AGAIN IN CASE THIS STEP IS RUN ON ITS OWN.
            # Copy fasta files to /tmp/ directory.
            temp_fasta_dir_containment = "/tmp/containment_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment)
            copy_fasta_tmp_cmd = "cp %s/*.fasta %s/" % (dedup_database_dir, temp_fasta_dir_containment)
            #os.system(copy_fasta_tmp_cmd)
            call("%s" % copy_fasta_tmp_cmd, logger)

            filenames_containments_array = []
            with open("%s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % dedup_database_dir) as fp:
                for line in fp:
                    line = line.strip()
                    # line = dedup_database_dir + "/" + line
                    line = temp_fasta_dir_containment + "/" + line
                    filenames_containments_array.append(line)
            fp.close()

            # Parse Nucmer Coordinates file for extracted aligned regions
            command_array = generate_parse_containments_jobs(args.jobrun, filenames_containments_array, filenames_containments_array, temp_cmd, logger)
            if args.jobrun:
                create_job_parse(args.jobrun, command_array, logger)

            # Read in the Matrix; find alignments with containments
            # os.system("cd %s" % dedup_database_dir)
            # os.system("ls *.fasta | sed 's/.fasta//g' > %s/rownames" % dedup_database_dir)
            # newline = "\n"
            # os.system(
            #     "ls *.fasta | sed 's/.fasta//g' | tr '\\n' '\\t' | sed 's/^/\\t/g' | sed 's/$/\\n/g' > %s/header" % (
            #         dedup_database_dir))
            # os.system("paste %s/rownames %s/*.score > %s/containment_matrix_temp.csv" % (
            # dedup_database_dir, dedup_database_dir, dedup_database_dir))
            # os.system("cat %s/header %s/containment_matrix_temp.csv > %s/containment_matrix.csv" % (
            # dedup_database_dir, dedup_database_dir, dedup_database_dir))

            call("cd %s" % dedup_database_dir, logger)
            call("ls *.fasta | sed 's/.fasta//g' | grep -v 'containment' | grep -v 'Extract' > %s/rownames" % dedup_database_dir, logger)
            newline = "\n"
            call(
                "ls *.fasta | sed 's/.fasta//g' | tr '\\n' '\\t' | sed 's/^/\\t/g' | sed 's/$/\\n/g' > %s/header" % (
                    dedup_database_dir), logger)
            call("paste %s/rownames %s/*.score > %s/containment_matrix_temp.csv" % (
                dedup_database_dir, dedup_database_dir, dedup_database_dir), logger)
            call("cat %s/header %s/containment_matrix_temp.csv > %s/containment_matrix.csv" % (
                dedup_database_dir, dedup_database_dir, dedup_database_dir), logger)

            keep_logging(
                'Reading containment matrix to determine fragments contained in another larger fragment: %s/containment_matrix.csv\n' % dedup_database_dir,
                'Reading containment matrix to determine fragments contained in another larger fragment: %s/containment_matrix.csv\n' % dedup_database_dir,
                logger,
                'info')

            c_reader = csv.reader(open('%s/containment_matrix.csv' % dedup_database_dir, 'r'), delimiter='\t')
            columns = list(zip(*c_reader))
            counts = 1
            end = len(filenames_containments_array) + 1
            #print end
            hit_dict = defaultdict(list)
            # change this array name to something suitable
            v = []

            # Pending: change temporarily assigned variables
            for i in xrange(1, end, 1):
                # Go through each column of containment_matrix.csv;
                # check which fragments have a score > 0.99.
                # Append these fragment names in cluster_rep array.
                # Each column will have its own cluster_rep list.

                # generate cluster_rep for each column.
                cluster_rep = []
                for ind, val in enumerate(columns[i]):
                    if ":" not in val:
                        if float(val) >= 0.99:
                            if columns[i][0] != columns[0][ind]:
                                # print str(columns[i][0]) + "," + str(columns[0][ind])
                                if columns[0][ind] not in cluster_rep:
                                    cluster_rep.append(columns[0][ind])
                #print cluster_rep

                # Find the fragment with longest length in cluster_rep list;
                # This fragment will be considered a representative for this cluster.
                hit_length = 0
                hit_name = ""
                for hits in cluster_rep:
                    hits_split = hits.split(':')
                    if hits_split[3] > hit_length:
                        hit_length = hits_split[3]
                        hit_name = hits
                # hit_name will be the fragment with longest length; this will be the representative for this cluster. All other fragments similar to this will be discarded
                # WTF is happening here?
                # hit_name is the name of representative key.
                for hits in cluster_rep:
                    if hit_name not in hit_dict.values():
                        hit_dict[hit_name].append(hits)
                        if hits != hit_name:
                            v.append(hits)

            #print "length of hit_dict: %s" % len(hit_dict)
            #print "length of hits: %s" % len(set(v))

            # Put all containment fragments key:value pairs in containment_rep.txt file for testing purposes. Note the values are contained in the key, so only key should appear in containment_removed_matrix.csv file.
            f_rep = open("%s/containment_rep.txt" % dedup_database_dir, 'w+')
            for key in hit_dict.keys():
                if key not in v:
                    print_string = "Key:%s\tValues:" % key
                    for containments in hit_dict[key]:
                        #print str(containments)
                        print_string = print_string + "," + str(containments)
                    print_string = print_string + "\n"
                    f_rep.write(print_string)
                    # print str(key) + "\t" + str(hit_dict[key]) + "\t" + str(len(hit_dict[key])) + "\n"
            f_rep.close()

            f = open("%s/containment_removed_matrix.csv" % dedup_database_dir, 'w+')
            f_containments = open("%s/Only_containments_matrix.csv" % dedup_database_dir, 'w+')
            c_reader = csv.reader(open('%s/containment_matrix.csv' % dedup_database_dir, 'r'), delimiter='\t')
            for row in c_reader:
                p_str = ""
                if row[0] not in v:
                    for i in row:
                        p_str = p_str + i + ","
                    f.write(p_str + '\n')
                else:
                    for i in row:
                        p_str = p_str + i + ","
                    f_containments.write(p_str + '\n')

            f.close()
            f_containments.close()

            final_containment_removed_dir = "%s/Final_results/" % dedup_database_dir
            make_sure_path_exists(final_containment_removed_dir)

            f_containment_removed_db_filenames = open(
                "%s/containment_removed_db_filenames.txt" % final_containment_removed_dir, 'w+')
            c_reader_2 = csv.reader(open('%s/containment_removed_matrix.csv' % dedup_database_dir, 'r'), delimiter=',')
            for row in c_reader_2:
                next(c_reader_2, None)
                if row[0] != "":
                    f_containment_removed_db_filenames.write("%s.fasta\n" % str(row[0]))
            f_containment_removed_db_filenames.close()

            # Find linker sequences and break fragments, name individual broken fragments with a suffix
            keep_logging(
                'Searching for linker sequence in containment removed fragments',
                'Searching for linker sequence in containment removed fragments',
                logger,
                'info')

            temp_fasta_dir_containment_removed = "/tmp/containment_removed_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment_removed)
            # os.system("rm %s/*" % temp_fasta_dir_containment_removed)
            call("rm %s/*" % temp_fasta_dir_containment_removed, logger)

            sequence = "NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN"

            containment_removed_db_filenames = []

            with open("%s/containment_removed_db_filenames.txt" % final_containment_removed_dir) as fpp:
                # print "here"
                for line in fpp:
                    line = line.strip()
                    for rec in SeqIO.parse("%s/%s" % (dedup_database_dir, line), "fasta"):
                        if "NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN" in rec.seq:
                            keep_logging('Linker sequence found in: %s' % line, 'Linker sequence found in: %s' % line, logger, 'info')
                            keep_logging('Number of linker sequence found: %s' % rec.seq.count(sequence), 'Number of linker sequence found: %s' % rec.seq.count(sequence), logger, 'info')
                            count = 1
                            split_record = rec.seq.split('NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN')
                            for i in split_record:
                                if len(i) > 1999:
                                    base = (os.path.basename(line)).replace('.fasta', '')
                                    basefile = temp_fasta_dir_containment_removed + "/" + base + "_" + str(
                                        count) + ".fasta"
                                    # print basefile
                                    temp_open = open(basefile, 'w+')
                                    header = ">%s_%s" % (base, count)
                                    fasta_string = "%s\n%s" % (header, i)
                                    temp_open.write(fasta_string)
                                    containment_removed_db_filenames.append(basefile)
                                    # print len(i)
                                    count += 1
                                    temp_open.close()
                        else:
                            copy_fasta_tmp_cmd = "cp %s/%s %s/" % (
                            dedup_database_dir, line, temp_fasta_dir_containment_removed)
                            #os.system(copy_fasta_tmp_cmd)
                            call(copy_fasta_tmp_cmd, logger)
                            line = temp_fasta_dir_containment_removed + "/" + line
                            containment_removed_db_filenames.append(line)
            fpp.close()


            keep_logging(
                'No of deduped fragments before containment removal: %s' % len(filenames_containments_array),
                'No of deduped fragments before containment removal: %s' % len(filenames_containments_array),
                logger,
                'info')

            keep_logging(
                'No of deduped fragments after containment removal: %s' % len(containment_removed_db_filenames),
                'No of deduped fragments after containment removal: %s' % len(containment_removed_db_filenames),
                logger,
                'info')

            # Run Nucmer original fasta vs final fragments database
            temp_cmd = "%s/temp_commands" % final_containment_removed_dir
            command_array = run_nucmer_db_final(temp_cmd, logger)
            if args.jobrun:
                create_job_containments_final(args.jobrun, command_array, logger)

            temp_cmd = "%s/temp_commands_parse_containments.sh" % final_containment_removed_dir

            command_array = generate_parse_containments_final_jobs(args.jobrun, filenames_array, containment_removed_db_filenames, temp_cmd, logger)
            if args.jobrun:
                create_job_parse(args.jobrun, command_array, logger)

            # Read in the Matrix; find alignments with containments
            os.chdir(temp_fasta_dir_containment_removed)
            #os.system("ls *.fasta | sed 's/.fasta//g' > %s/rownames" % final_containment_removed_dir)
            call("ls *.fasta | sed 's/.fasta//g' > %s/rownames" % final_containment_removed_dir, logger)
            newline = "\n"

            # Prepare matrix from alignment score results.
            header = "\t"
            for fasta in filenames_array:
                header = header + fasta.replace('.fasta', '') + "\t"

            header_open = open("%s/header" % final_containment_removed_dir, 'w+')
            header_open.write(header + '\n')
            header_open.close()

            # os.system("paste %s/rownames %s/*.score > %s/Final_HGT_score_matrix_temp.csv" % (
            #     final_containment_removed_dir, final_containment_removed_dir, final_containment_removed_dir))
            # os.system("cat %s/header %s/Final_HGT_score_matrix_temp.csv > %s/Final_HGT_score_matrix.csv" % (
            #     final_containment_removed_dir, final_containment_removed_dir, final_containment_removed_dir))

            call("paste %s/rownames %s/*.score > %s/Final_HGT_score_matrix_temp.csv" % (
                final_containment_removed_dir, final_containment_removed_dir, final_containment_removed_dir), logger)
            call("cat %s/header %s/Final_HGT_score_matrix_temp.csv > %s/Final_HGT_score_matrix.csv" % (
                final_containment_removed_dir, final_containment_removed_dir, final_containment_removed_dir), logger)

            keep_logging(
                'Generating annotation metadata file from %s/Extracted_aligned_region_dedup_cluster_99.tsv.' % final_containment_removed_dir,
                'Generating annotation metadata file from %s/Extracted_aligned_region_dedup_cluster_99.tsv.' % final_containment_removed_dir,
                logger,
                'info')
            aligned_fragments_meta = database_directory + "/Extracted_aligned_region_dedup_cluster_99.tsv"
            aligned_fragments_meta_dict = {}
            with open("%s/Extracted_aligned_region_dedup_cluster_99.tsv" % database_directory) as fpp:
                for lines in fpp:
                    line_split = lines.split('\t')
                    if line_split[0] not in aligned_fragments_meta_dict.keys():
                        aligned_fragments_meta_dict[line_split[0]] = line_split[5]
                    else:
                         keep_logging(
                'Duplicate aligned fragment key found: %s' % line_split[0],
                'Duplicate aligned fragment key found: %s' % line_split[0],
                logger,
                'info')



            grep_meta_file = open(
                "%s/grep_meta.sh" % final_containment_removed_dir, 'w+')
            for frags in containment_removed_db_filenames:
                grep_meta_file.write("grep -w '%s' %s/Extracted_aligned_region_dedup_cluster_99.tsv\n" % (os.path.basename(frags.replace(".fasta", "")), database_directory))

            #os.system("bash %s/grep_meta.sh > %s/Final_HGT_score_matrix_meta.tsv" % (final_containment_removed_dir, final_containment_removed_dir))
            call("bash %s/grep_meta.sh > %s/Final_HGT_score_matrix_meta.tsv" % (
            final_containment_removed_dir, final_containment_removed_dir), logger)
            keep_logging(
                'Running: bash %s/grep_meta.sh > %s/Final_HGT_score_matrix_meta.tsv' % (final_containment_removed_dir, final_containment_removed_dir),
                'Running: bash %s/grep_meta.sh > %s/Final_HGT_score_matrix_meta.tsv' % (final_containment_removed_dir, final_containment_removed_dir),
                logger,
                'info')

            keep_logging(
                'The final HGT score matrix is: %s/Final_HGT_score_matrix.csv' % final_containment_removed_dir,
                'The final HGT score matrix is: %s/Final_HGT_score_matrix.csv' % final_containment_removed_dir,
                logger,
                'info')

        if "6" in steps_list:
            # This step is for testing/debugging purposes only.

            # Set up deduped database directories
            database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"

            final_containment_removed_dir = "%s/Final_results/" % dedup_database_dir
            make_sure_path_exists(final_containment_removed_dir)

            # GENERATE FASTA FILENAMES ARRAY
            temp_fasta_dir_containment_removed = "/tmp/containment_removed_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment_removed)
            #os.system("rm %s/*" % temp_fasta_dir_containment_removed)
            call("rm %s/*" % temp_fasta_dir_containment_removed, logger)


            sequence = "NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN"

            containment_removed_db_filenames = []
            containment_removed_db_filenames_no_linkers = open("containment_removed_db_filenames_no_linkers.txt", 'w+')

            with open("%s/containment_removed_db_filenames.txt" % final_containment_removed_dir) as fpp:
                # print "here"
                for line in fpp:
                    line = line.strip()
                    for rec in SeqIO.parse("%s/%s" % (dedup_database_dir, line), "fasta"):
                        if "NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN" in rec.seq:
                            keep_logging('Linker sequence found in: %s' % line, 'Linker sequence found in: %s' % line, logger, 'info')
                            keep_logging('Number of linker sequence found: %s' % rec.seq.count(sequence), 'Number of linker sequence found: %s' % rec.seq.count(sequence), logger, 'info')
                            count = 1
                            split_record = rec.seq.split('NNNNNCATTCCATTCATTAATTAATTAATGAATGAATGNNNNN')
                            for i in split_record:
                                if len(i) > 1999:
                                    base = (os.path.basename(line)).replace('.fasta', '')
                                    basefile = temp_fasta_dir_containment_removed + "/" + base + "_" + str(count) + ".fasta"
                                    # print basefile
                                    temp_open = open(basefile, 'w+')
                                    header = ">%s_%s" % (base, count)
                                    fasta_string = "%s\n%s" % (header, i)
                                    temp_open.write(fasta_string)
                                    containment_removed_db_filenames.append(basefile)
                                    containment_removed_db_filenames_no_linkers.write('%s\n' % basefile)
                                    # print len(i)
                                    count += 1
                                    temp_open.close()
                        else:
                            copy_fasta_tmp_cmd = "cp %s/%s %s/" % (dedup_database_dir, line, temp_fasta_dir_containment_removed)
                            containment_removed_db_filenames_no_linkers.write('%s\n' % line)
                            os.system(copy_fasta_tmp_cmd)
                            line = temp_fasta_dir_containment_removed + "/" + line
                            containment_removed_db_filenames.append(line)
            fpp.close()

        if args.remove_temp == "yes":
            remove_temp_directories = "rm -rf /tmp/fasta_file /tmp/$PBS_JOBID"
            #os.system(remove_temp_directories)
            keep_logging(
                'Running: %s' % (remove_temp_directories),
                'Running: %s' % (remove_temp_directories),
                logger,
                'debug')
            keep_logging(
                'Removing temporary directories: %s, %s' % ("/tmp/fasta_file", "/tmp/$PBS_JOBID"),
                'Removing temporary directories: %s, %s' % ("/tmp/fasta_file", "/tmp/$PBS_JOBID"),
                logger,
                'info')

    else:
        keep_logging('Exiting. Require -steps argument to run the script.', 'Exiting. Require -steps argument to run the script.', logger, 'exception')