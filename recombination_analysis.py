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
from argparse import RawTextHelpFormatter


parser = argparse.ArgumentParser(description='Pipeline: Recombination/HGT Analysis.\n The pipeline takes a list of fasta files and aligns All-vs-All using Nucmer.\n Parses Nucmer coordinates results and input gff/bed annotation files to extract regions that matches user provided percent identity and aligned region(bp) parameters.\n Makes a reference database of these extracted aligned regions by deduplicating and removing containments.\n Performs nucmer alignment between input fasta file and extracted regions to generate a matrix of aligned bases ratio.\n', formatter_class=RawTextHelpFormatter)
required = parser.add_argument_group('Required arguments')
required.add_argument('-filename', action='store', dest="filename", help='This file contains a list of fasta filenames. For Genome coordinate consistency, make sure the fasta files are in a pseudomolecule assembly format i.e contigs stiched together into one assembly.', required=True)
required.add_argument('-out', action='store', dest="out", help='Output directory to save the results', required=True)
required.add_argument('-prokka_dir', action='store', dest="prokka_dir", help='Directory containing results of Prokka annotation pipeline.', required=True)
required.add_argument('-jobrun', action='store', dest="jobrun",
                    help='Type of job to run. Running on Cluster or Parallelly at local system or locally (default): cluster, parallel-local, local')
required.add_argument('-dir', action='store', dest="dir", help='Directory containing fasta files specified in -filename list', required=True)
required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique Analysis Name to save results with this prefix', required=True)
optional = parser.add_argument_group('Optional arguments')
optional.add_argument('-filename_db', action='store', dest="filename_db", help='filenames of nucmer db', required=False)
optional.add_argument('-dir_db', action='store', dest="dir_db", help='directory of nucmer db fasta files', required=False)
optional.add_argument('-matrix', action='store', dest="matrix", help='Matrix to parse and remove containments', required=False)
optional.add_argument('-steps', action='store', dest="steps",
                    help='Analysis Steps to be performed. Default: All'
                         '\n1: Align all fasta input file against each other using Nucmer.'
                         '\n2: Parses the Nucmer generated coordinates files and generate annotated results of aligned fragments.'
                         '\n3: Generate a database of these extracted aligned regions by deduplicating and removing containments.'
                         '\n4: Performs nucmer alignment between input fasta file and extracted regions to generate a matrix of aligned bases ratio.')
args = parser.parse_args()


def create_job(jobrun, commands_list):

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
        # job_directory = args.out + "/" + "temp_jobs"
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
        #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
        #pbs_dir = job_directory + "/nucmer_job_command_*.pbs"
        pbs_dir = job_directory + "/*/*_parallel_cluster.pbs"
        pbs_scripts = glob.glob(pbs_dir)
        for i in pbs_scripts:
            print "Running: qsub %s" % i
            os.system("bash %s" % i)
            #os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_file = commands_list
        print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)


    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            os.system(i)

def create_job_parse(jobrun, commands_list):

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
            print "Running: qsub %s" % i
            os.system("bash %s" % i)
            #os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
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
            os.system(i)

def create_job_containments(jobrun, commands_list):

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
            print "Running: qsub %s" % i
            #os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_file = commands_list
        print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)


    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            os.system(i)

def create_job_containments_final(jobrun, commands_list):

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
            print "Running: qsub %s" % i
            #os.system("qsub %s" % i)

    elif jobrun == "parallel-local":
        """
        Generate a Command list of each job and run it in parallel on different cores available on local system
        """
        command_file = commands_list
        print len(command_file)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_file)


    else:
        """
        Generate a Command list of each job and run it on local system one at a time
        """
        command_file = commands_list
        for i in command_file:
            os.system(i)

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
            print "Errors in output folder path! please change the output path or analysis name\n"
            exit()

def run_nucmer(temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_array:
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (args.out, filebase.replace('.fasta', ''))
        os.system(mkdir_command)
        for file_again in filenames_array:
            file_again_base = os.path.basename(file_again)
            prefix = filebase.replace('.fasta', '') + "_" + file_again_base.replace('.fasta', '')
            command = "nucmer --maxmatch --prefix=%s %s %s &>/dev/null && show-coords -r -c -l -T %s.delta > %s.coords && rm %s.delta && mv %s.coords %s/%s\n" % (prefix, file, file_again, prefix, prefix, prefix, prefix, args.out, filebase.replace('.fasta', ''))
            command_array.append(command)
            f.write(command)
    f.close()
    return command_array

def run_nucmer_parallel(temp_cmd):
    #f = open(temp_cmd, 'w+')
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

        print filebase_command_cluster
        print filebase_command_local
        print run_parallel_script_cluster
        print run_parallel_script_local

        os.system(mkdir_command)

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

        for i in pbs_scripts_cluster:
            print "Running: qsub %s" % i

    #f.close()

    return command_array_local, command_array_cluster

def run_nucmer_parallel_containment(temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_containments_array:
        file_command_array = []
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (dedup_database_dir, filebase.replace('.fasta', ''))
        filebase_command = "%s/%s/%s_commands.sh" % (dedup_database_dir, filebase.replace('.fasta', ''), filebase)
        run_parallel_script = "%s/%s/%s_parallel.pbs" % (dedup_database_dir, filebase.replace('.fasta', ''), filebase)
        print filebase_command
        print run_parallel_script
        os.system(mkdir_command)
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
        for i in pbs_scripts:
            print "Running: qsub %s" % i
    f.close()

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
        print filebase_command
        print run_parallel_script
        os.system(mkdir_command)
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

def run_nucmer_db_final(temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    for file in filenames_array:
        file_command_array = []
        filebase = os.path.basename(file)
        mkdir_command = "mkdir %s/%s" % (final_containment_removed_dir, filebase.replace('.fasta', ''))
        filebase_command = "%s/%s/%s_commands.sh" % (final_containment_removed_dir, filebase.replace('.fasta', ''), filebase)
        run_parallel_script = "%s/%s/%s_parallel.pbs" % (final_containment_removed_dir, filebase.replace('.fasta', ''), filebase)
        print filebase_command
        print run_parallel_script
        os.system(mkdir_command)
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
    os.system(i)
    done = "done"
    return done

def generate_parse_coord_aggregate_jobs(jobrun, filenames_array, temp_cmd):
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

def generate_parse_coord_db_aggregate_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd):
    f = open(temp_cmd, 'w+')
    command_array = []
    job_directory = args.out + "/" + "temp_jobs"
    make_sure_path_exists(job_directory)
    for folder in filenames_array:
        cmd = "~/anaconda/bin/python /nfs/esnitkin/bin_group/scripts/Scripts_v2.0/recombination_analysis/run_nucmer_coordinates_db.py -folder %s -out %s -dir %s -prokka_dir %s\n" % (folder, args.out, args.dir, args.prokka_dir)
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
    print len(command_array)
    return command_array

def generate_parse_containments_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd):
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
    print len(command_array)
    return command_array

def generate_parse_containments_final_jobs(jobrun, filenames_array, filenames_db_array, temp_cmd):
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
    print len(command_array)
    return command_array

#Main Steps
if __name__ == '__main__':
    """Start Timer"""
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now().strftime('%Y-%m-%d')

    logger = generate_logger(args.out, args.analysis_name, start_time_2)
    keep_logging('The Script started at: %s\n' % start_time, 'The Script started at: %s\n' % start_time, logger, 'info')


    print "\nThis pipeline will run the following steps:\n\n" \
          "1: Align all fasta input file against each other using Nucmer.\n" \
          "2: Parses the Nucmer generated coordinates files and generate annotated results of aligned fragments.\n" \
          "3: Generate a database of these extracted aligned regions by deduplicating and removing containments.\n" \
          "4: Performs nucmer alignment between input fasta file and extracted regions to generate a matrix of alignment score.\n" \

    # GENERATE FASTA FILENAMES ARRAY
    filenames_array = []

    # Copy fasta files to /tmp/ directory.
    temp_fasta_dir = "/tmp/fasta_file"
    make_sure_path_exists(temp_fasta_dir)
    copy_fasta_tmp_cmd = "cp %s/* %s/" % (args.dir, temp_fasta_dir)
    os.system(copy_fasta_tmp_cmd)

    with open(args.filename) as fp:
        for line in fp:
            line = line.strip()
            #line = args.dir + "/" + line
            line = temp_fasta_dir + "/" + line
            filenames_array.append(line)

    # CHECK FOR STEPS ARGUMENTS

    # Run pipeline steps
    if args.steps:
        # split values provided with -steps argument and decide the starting point of pipeline
        steps_list = args.steps.split(',')

        if "1" in steps_list:
            # Run nucmer on all the input fasta pseudomolecule assembly fasta file. compare All-vs-All.
            keep_logging('Running Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n', 'Running Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n', logger,
                         'info')
            make_sure_path_exists(args.out)
            temp_cmd = "%s/temp_commands" % args.out
            command_array_local, command_array_cluster = run_nucmer_parallel(temp_cmd)
            if args.jobrun:
                if args.jobrun == "parallel-local":
                    create_job(args.jobrun, command_array_local)
                elif args.jobrun == "cluster":
                    create_job(args.jobrun, command_array_cluster)

            keep_logging('End Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n',
                         'End Step 1: Run Nucmer on all input fasta files to compare All-vs-All\n', logger,
                         'info')

        if "2" in steps_list:
            # Step 2
            keep_logging('Running Step 2: Parsing Nucmer Coordinate results generated by Nucmer Step 1 and extracting all the aligned regions.',
                         'Running Step 2: Parsing Nucmer Coordinate results generated by Nucmer Step 1 and extracting all the aligned regions.', logger,
                         'info')
            temp_cmd_parse_coord = "%s/temp_commands_parse_coord" % args.out
            command_array = generate_parse_coord_aggregate_jobs(args.jobrun, filenames_array, temp_cmd_parse_coord)
            if args.jobrun:
                create_job_parse(args.jobrun, command_array)

            # Move Aggregate score results. No use of these results. Remove when done. Pending...
            aggregate_results = args.out + "/aggregate_results/"
            make_sure_path_exists(aggregate_results)
            os.system("mv %s/*.aligned %s/*.score %s" % (args.out, args.out, aggregate_results))

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
            print "Running: %s" % dedupe_cmd
            keep_logging('Running dedupe command: %s\n' % dedupe_cmd,
                         'Running dedupe command: %s\n' % dedupe_cmd, logger,
                         'info')
            os.system(dedupe_cmd)

            # Run dedupe to remove containment sequences from deduplicated preliminary extracted aligned region database results.
            dedupe_containment_cmd = "/nfs/esnitkin/bin_group/bbmap/dedupe.sh in=%s/Extracted_aligned_region_dedup.fasta out=%s/Extracted_aligned_region_dedup_cluster_99.fasta minidentity=99" % (
            database_directory, database_directory)
            os.system(dedupe_containment_cmd)
            keep_logging('Running dedupe containment command: %s\n' % dedupe_containment_cmd,
                         'Running dedupe containment command: %s\n' % dedupe_containment_cmd, logger,
                         'info')
            # Prepare or extract metadata annotations for sequences that remained after dedupe deduplicate and containment removal steps
            tab = "\\t"
            get_meta_annotations = "for i in `grep '>' %s/Extracted_aligned_region_dedup_cluster_99.fasta | sed 's/>//g'`; do echo \"grep '^$i' %s/Extracted_aligned_region.fasta_meta.tsv | awk -F'\\t' -v OFS='%s' '{print \$1,\$2,\$3,\$4,\$5,\$6}'\"; done > %s/grep_meta.sh" % (
                database_directory, database_directory, tab, database_directory)

            keep_logging('Extracted Metadata annotations for deduplicated and containment removed sequences: %s\n' % get_meta_annotations,
                         'Extracted Metadata annotations for deduplicated and containment removed sequences: %s\n' % get_meta_annotations, logger,
                         'info')
            keep_logging(
                'Extracted Metadata annotations for deduplicated and containment removed sequences: bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv\n' % (database_directory, database_directory),
                'Extracted Metadata annotations for deduplicated and containment removed sequences: bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv\n' % (database_directory, database_directory),
                logger,
                'info')

            os.system(get_meta_annotations)
            os.system("bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv" % (
            database_directory, database_directory))

            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
            make_sure_path_exists(dedup_database_dir)

            get_deduped_extracted_aligned_region_database_files = "for i in `awk -F'\\t' '{print $1}' %s/Extracted_aligned_region_dedup_cluster_99.tsv`; do  echo \"$i\" > %s/$i.txt; done" % (
            database_directory, dedup_database_dir)

            print "Running: %s" % get_deduped_extracted_aligned_region_database_files
            keep_logging(
                'Extract sequence names from deduplicated and containment removed sequences: %s\n' % get_deduped_extracted_aligned_region_database_files,
                'Extract sequence names from deduplicated and containment removed sequences:: %s\n' % get_deduped_extracted_aligned_region_database_files,
                logger,
                'info')

            os.system(get_deduped_extracted_aligned_region_database_files)

            extract_regions_from_file = "for i in `ls %s/*.txt`; do base=`echo $i | sed 's/.txt//g'`; ~/bin/seqtk/seqtk subseq %s/Extracted_aligned_region_dedup_cluster_99.fasta $i > $base.fasta; done" % (
            dedup_database_dir, database_directory)

            keep_logging(
                'Extract individual sequences from deduplicated and containment removed sequence mfa file: %s\n' % extract_regions_from_file,
                'Extract individual sequences from deduplicated and containment removed sequence mfa file: %s\n' % extract_regions_from_file,
                logger,
                'info')
            os.system(extract_regions_from_file)

        if "4" in steps_list:
            keep_logging('Running Step 4: Remove containment sequences from the preliminary database fragments. Part one: Run Nucmer all deduped extracted sequences against each other.',
                         'Running Step 4: Remove containment sequences from the preliminary database fragments. Part one: Run Nucmer all deduped extracted sequences against each other.', logger,
                         'info')

            # Containment sequence: a sequence which is exactly/99% similar to another sequence or a sequence which is contained inside a larger sequence
            # To remove these containment sequences, run nucmer on extracted fragments vs extracted fragments
            database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
            os.chdir(dedup_database_dir)

            get_extracted_aligned_region_dedup_cluster_99_fasta_files = "ls *.fasta > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % (
                dedup_database_dir)

            keep_logging(
                'Save deduped extracted aligned individual fasta filenames: ls *.fasta > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt\n' % dedup_database_dir,
                'Save deduped extracted aligned individual fasta filenames: ls *.fasta > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt\n' % dedup_database_dir,
                logger,
                'info')
            os.system(get_extracted_aligned_region_dedup_cluster_99_fasta_files)

            # Run Step 1 nucmer method again but with only Extracted_aligned_region_dedup_cluster_99_fasta_filenames
            output_folder = dedup_database_dir

            # GENERATE FASTA FILENAMES ARRAY
            # Copy fasta files to /tmp/ directory.
            temp_fasta_dir_containment = "/tmp/containment_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment)
            copy_fasta_tmp_cmd = "cp %s/*.fasta %s/" % (dedup_database_dir, temp_fasta_dir_containment)
            os.system(copy_fasta_tmp_cmd)


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
            command_array = run_nucmer_parallel_containment(temp_cmd)
            if args.jobrun:
                create_job_containments(args.jobrun, command_array)

        if "5" in steps_list:
            keep_logging(
                'Running Step 5: Remove containment sequences from the preliminary database fragments. Part two. Parse Nucmer coordinates file of deduped extracted aligned region vs each other file and find contained sequences.',
                'Running Step 5: Remove containment sequences from the preliminary database fragments. Part two. Parse Nucmer coordinates file of deduped extracted aligned region vs each other file and find contained sequences.',
                logger,
                'info')
            temp_cmd = "%s/temp_commands_parse_containments.sh" % args.out

            # Definition of Containment sequence: a sequence which is exactly/99% similar to another sequence or a smaller sequence which is contained inside a larger sequence.
            # To remove these containment sequences, run nucmer on extracted fragments vs  extracted fragments
            # Parse the nucmer generated coordinates file, generate a matrix of alignment ratios, find out...

            # Set up deduped database directories
            database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
            dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
            os.chdir(dedup_database_dir)


            # GENERATE FASTA FILENAMES ARRAY AGAIN IN CASE THIS STEP IS RUN ON ITS OWN.
            # Copy fasta files to /tmp/ directory.
            temp_fasta_dir_containment = "/tmp/containment_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment)
            copy_fasta_tmp_cmd = "cp %s/*.fasta %s/" % (dedup_database_dir, temp_fasta_dir_containment)
            os.system(copy_fasta_tmp_cmd)

            filenames_containments_array = []
            with open("%s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % dedup_database_dir) as fp:
                for line in fp:
                    line = line.strip()
                    # line = dedup_database_dir + "/" + line
                    line = temp_fasta_dir_containment + "/" + line
                    filenames_containments_array.append(line)
            fp.close()

            # Parse Nucmer Coordinates file for extracted aligned regions
            command_array = generate_parse_containments_jobs(args.jobrun, filenames_containments_array, filenames_containments_array, temp_cmd)
            if args.jobrun:
                create_job_parse(args.jobrun, command_array)

            # Read in the Matrix; find alignments with containments
            os.system("cd %s" % dedup_database_dir)
            os.system("ls *.fasta | sed 's/.fasta//g' > %s/rownames" % dedup_database_dir)
            newline = "\n"
            os.system(
                "ls *.fasta | sed 's/.fasta//g' | tr '\\n' '\\t' | sed 's/^/\\t/g' | sed 's/$/\\n/g' > %s/header" % (
                    dedup_database_dir))
            os.system("paste %s/rownames %s/*.score > %s/containment_matrix_temp.csv" % (
            dedup_database_dir, dedup_database_dir, dedup_database_dir))
            os.system("cat %s/header %s/containment_matrix_temp.csv > %s/containment_matrix.csv" % (
            dedup_database_dir, dedup_database_dir, dedup_database_dir))

            c_reader = csv.reader(open('%s/containment_matrix.csv' % dedup_database_dir, 'r'), delimiter='\t')
            columns = list(zip(*c_reader))
            counts = 1
            end = len(filenames_containments_array) + 1
            print end
            hit_dict = defaultdict(list)
            v = []
            for i in xrange(1, end, 1):
                cluster_rep = []
                for ind, val in enumerate(columns[i]):
                    if ":" not in val:
                        if float(val) >= 0.99:
                            if columns[i][0] != columns[0][ind]:
                                # print str(columns[i][0]) + "," + str(columns[0][ind])
                                if columns[0][ind] not in cluster_rep:
                                    cluster_rep.append(columns[0][ind])
                print cluster_rep
                hit_length = 0
                hit_name = ""
                for hits in cluster_rep:
                    hits_split = hits.split(':')
                    if hits_split[3] > hit_length:
                        hit_length = hits_split[3]
                        hit_name = hits
                for hits in cluster_rep:
                    if hit_name not in hit_dict.values():
                        hit_dict[hit_name].append(hits)
                        if hits != hit_name:
                            v.append(hits)
            f_rep = open("%s/containment_rep.txt" % dedup_database_dir, 'w+')
            for key in hit_dict.keys():
                if key not in v:
                    print_string = "%s\t" % key
                    for containments in hit_dict[key]:
                        print str(containments)
                        print_string = print_string + "," + str(containments)
                    print_string = print_string + "\n"
                    f_rep.write(print_string)
                    # print str(key) + "\t" + str(hit_dict[key]) + "\t" + str(len(hit_dict[key])) + "\n"
            f_rep.close()
            # print len(v)
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


            # GENERATE FASTA FILENAMES ARRAY
            temp_fasta_dir_containment_removed = "/tmp/containment_removed_fasta_files"
            make_sure_path_exists(temp_fasta_dir_containment_removed)


            containment_removed_db_filenames = []
            #print "here"
            with open("%s/containment_removed_db_filenames.txt" % final_containment_removed_dir) as fpp:
                #print "here"
                for line in fpp:
                    line = line.strip()
                    print line
                    copy_fasta_tmp_cmd = "cp %s/%s %s/" % (dedup_database_dir, line, temp_fasta_dir_containment_removed)
                    os.system(copy_fasta_tmp_cmd)
                    #line = dedup_database_dir + "/" + line
                    line = temp_fasta_dir_containment_removed + "/" + line
                    containment_removed_db_filenames.append(line)
            fpp.close()

            print len(containment_removed_db_filenames)
            print containment_removed_db_filenames

            # Run Nucmer original fasta vs final fragments database
            temp_cmd = "%s/temp_commands" % final_containment_removed_dir
            command_array = run_nucmer_db_final(temp_cmd)
            if args.jobrun:
                create_job_containments_final(args.jobrun, command_array)

            temp_cmd = "%s/temp_commands_parse_containments.sh" % final_containment_removed_dir

            command_array = generate_parse_containments_final_jobs(args.jobrun, filenames_array, containment_removed_db_filenames, temp_cmd)
            if args.jobrun:
                create_job_parse(args.jobrun, command_array)

        remove_temp_directories = "rm -rf /tmp/fasta_file /tmp/$PBS_JOBID"
        #os.system(remove_temp_directories)
        print "%s" % remove_temp_directories
    else:
        keep_logging('Running All Default Pipeline Steps...\n', 'Running All Default Pipeline Steps...\n', logger,
                     'info')
        keep_logging('Running Step 1: Nucmer All-vs-All', 'Running Step 1: Nucmer All-vs-All', logger,
                     'info')
        # Step 1
        make_sure_path_exists(args.out)
        temp_cmd = "%s/temp_commands" % args.out
        command_array = run_nucmer_parallel(temp_cmd)
        if args.jobrun:
            create_job(args.jobrun, command_array)

        # Step 2
        keep_logging('Running Step 2: Parsing Nucmer All-vs-All coordinates results',
                     'Running Step 2: Parsing Nucmer All-vs-All coordinates results', logger,
                     'info')
        temp_cmd_parse_coord = "%s/temp_commands_parse_coord" % args.out
        command_array = generate_parse_coord_aggregate_jobs(args.jobrun, filenames_array, temp_cmd_parse_coord)
        if args.jobrun:
            create_job_parse(args.jobrun, command_array)

        # Move Aggregate score results
        aggregate_results = args.out + "/aggregate_results/"
        make_sure_path_exists(aggregate_results)
        os.system("mv %s/*.aligned %s/*.score %s" % (args.out, args.out, aggregate_results))

        # Add Nursing Home steps here later

        # Step 3
        keep_logging('Running Step 3: Generate a preliminary database from the extracted aligned region.',
                     'Running Step 3: Generate a preliminary database from the extracted aligned region.', logger,
                     'info')
        # A lot of the code and path to executables in this step are hard coded here
        database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
        dedupe_cmd = "/nfs/esnitkin/bin_group/bbmap/dedupe.sh in=%s/Extracted_aligned_region.fasta out=%s/Extracted_aligned_region_dedup.fasta" % (database_directory, database_directory)
        print "Running: %s" % dedupe_cmd
        os.system(dedupe_cmd)
        dedupe_containment_cmd = "/nfs/esnitkin/bin_group/bbmap/dedupe.sh in=%s/Extracted_aligned_region_dedup.fasta out=%s/Extracted_aligned_region_dedup_cluster_99.fasta minidentity=99" % (database_directory, database_directory)
        print "Running: %s" % dedupe_containment_cmd
        os.system(dedupe_containment_cmd)
        tab = "\\t"
        get_meta_annotations = "for i in `grep '>' %s/Extracted_aligned_region_dedup_cluster_99.fasta | sed 's/>//g'`; do echo \"grep '^$i' %s/Extracted_aligned_region.fasta_meta.tsv | awk -F'\\t' -v OFS='%s' '{print \$1,\$2,\$3,\$4,\$5,\$6}'\"; done > %s/grep_meta.sh" % (
        database_directory, database_directory, tab, database_directory)

        print "Running: %s" % get_meta_annotations
        print "Running: bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv" % (database_directory, database_directory)
        os.system(get_meta_annotations)
        os.system("bash %s/grep_meta.sh > %s/Extracted_aligned_region_dedup_cluster_99.tsv" % (database_directory, database_directory))

        dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
        make_sure_path_exists(dedup_database_dir)

        get_deduped_extracted_aligned_region_database_files = "for i in `awk -F'\\t' '{print $1}' %s/Extracted_aligned_region_dedup_cluster_99.tsv`; do  echo \"$i\" > %s/$i.txt; done" % (database_directory, dedup_database_dir)

        print "Running: %s" % get_deduped_extracted_aligned_region_database_files
        os.system(get_deduped_extracted_aligned_region_database_files)

        extract_regions_from_file = "for i in `ls %s/*.txt`; do base=`echo $i | sed 's/.txt//g'`; ~/bin/seqtk/seqtk subseq %s/Extracted_aligned_region_dedup_cluster_99.fasta $i > $base.fasta; done" % (dedup_database_dir, database_directory)
        print "Running: %s" % extract_regions_from_file
        os.system(extract_regions_from_file)

        ####################################################
        # Step 4
        database_directory = args.out + "/2018_Recombination_analysis_Results/All_vs_All/database/"
        dedup_database_dir = database_directory + "/deduped_extracted_aligned_region_database/"
        os.chdir(dedup_database_dir)

        get_extracted_aligned_region_dedup_cluster_99_fasta_files = "ls *.fasta > %s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % (dedup_database_dir)
        print "Running: %s" % get_extracted_aligned_region_dedup_cluster_99_fasta_files
        os.system(get_extracted_aligned_region_dedup_cluster_99_fasta_files)

        # Run Step 1 again but with only Extracted_aligned_region_dedup_cluster_99_fasta_filenames
        output_folder = dedup_database_dir

        # GENERATE FASTA FILENAMES ARRAY
        filenames_containments_array = []
        with open("%s/Extracted_aligned_region_dedup_cluster_99_fasta_filenames.txt" % dedup_database_dir) as fp:
            for line in fp:
                line = line.strip()
                line = dedup_database_dir + "/" + line
                filenames_containments_array.append(line)
        fp.close()
        # temp_cmd = "%s/temp_commands" % output_folder
        # command_array = run_nucmer_parallel_containment(temp_cmd)
        # if args.jobrun:
        #     create_job_containments(args.jobrun, command_array)


        # Step 5
        temp_cmd = "%s/temp_commands_parse_containments.sh" % args.out

        # command_array = generate_parse_containments_jobs(args.jobrun, filenames_containments_array, filenames_containments_array, temp_cmd)
        # if args.jobrun:
        #     create_job_parse(args.jobrun, command_array)

        # Read in the Matrix
        os.system("cd %s" % dedup_database_dir)
        os.system("ls *.fasta | sed 's/.fasta//g' > %s/rownames" % dedup_database_dir)
        newline = "\n"
        os.system("ls *.fasta | sed 's/.fasta//g' | tr '\\n' '\\t' | sed 's/^/\\t/g' | sed 's/$/\\n/g' > %s/header" % (dedup_database_dir))
        os.system("paste %s/rownames %s/*.score > %s/containment_matrix_temp.csv" % (dedup_database_dir, dedup_database_dir, dedup_database_dir))
        os.system("cat %s/header %s/containment_matrix_temp.csv > %s/containment_matrix.csv" % (dedup_database_dir, dedup_database_dir, dedup_database_dir))

        c_reader = csv.reader(open('%s/containment_matrix.csv' % dedup_database_dir, 'r'), delimiter='\t')
        columns = list(zip(*c_reader))
        counts = 1
        end = len(filenames_containments_array) + 1
        print end
        hit_dict = defaultdict(list)
        v = []
        for i in xrange(1, end, 1):
            cluster_rep = []
            for ind, val in enumerate(columns[i]):
                if ":" not in val:
                    if float(val) >= 0.99:
                        if columns[i][0] != columns[0][ind]:
                            # print str(columns[i][0]) + "," + str(columns[0][ind])
                            if columns[0][ind] not in cluster_rep:
                                cluster_rep.append(columns[0][ind])
            print cluster_rep
            hit_length = 0
            hit_name = ""
            for hits in cluster_rep:
                hits_split = hits.split(':')
                if hits_split[3] > hit_length:
                    hit_length = hits_split[3]
                    hit_name = hits
            for hits in cluster_rep:
                if hit_name not in hit_dict.values():
                    hit_dict[hit_name].append(hits)
                    if hits != hit_name:
                        v.append(hits)
        f_rep = open("./containment_rep.txt", 'w+')
        for key in hit_dict.keys():
            if key not in v:
                print_string = "%s\t" % key
                for containments in hit_dict[key]:
                    print str(containments)
                    print_string = print_string + "," + str(containments)
                print_string = print_string + "\n"
                f_rep.write(print_string)
                # print str(key) + "\t" + str(hit_dict[key]) + "\t" + str(len(hit_dict[key])) + "\n"
        f_rep.close()
        #print len(v)
        f = open("./containment_removed_matrix.csv", 'w+')
        f_containments = open("./Only_containments_matrix.csv", 'w+')
        c_reader = csv.reader(open('containment_matrix.csv', 'r'), delimiter='\t')
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

        final_containment_removed_dir = "%s/Final_results/" % dedup_database_dir
        make_sure_path_exists(final_containment_removed_dir)

        f_containment_removed_db_filenames = open("%s/containment_removed_db_filenames.txt" % final_containment_removed_dir, 'w+')
        c_reader_2 = csv.reader(open('%s/containment_removed_matrix.csv' % dedup_database_dir, 'r'), delimiter=',')
        for row in c_reader_2:
            next(c_reader_2, None)
            if row[0] != "":
                f_containment_removed_db_filenames.write("%s.fasta\n" % str(row[0]))

        # GENERATE FASTA FILENAMES ARRAY
        #########Stuck here
        containment_removed_db_filenames = []
        print "here"
        with open("%s/containment_removed_db_filenames.txt" % final_containment_removed_dir) as fpp:
            print "here"
            for line in fpp:
                line = line.strip()
                print line
                print "here"
                line = dedup_database_dir + "/" + line
                containment_removed_db_filenames.append(line)
        fpp.close()
        print len(containment_removed_db_filenames)
        # temp_cmd = "%s/temp_commands" % final_containment_removed_dir
        # command_array = run_nucmer_db_final(temp_cmd)
        # if args.jobrun:
        #     create_job_containments_final(args.jobrun, command_array)

        # temp_cmd = "%s/temp_commands_parse_containments.sh" % args.out
        #
        # command_array = generate_parse_containments_final_jobs(args.jobrun, filenames_containments_array, filenames_containments_array, temp_cmd)
        # if args.jobrun:
        #     create_job_parse(args.jobrun, command_array)