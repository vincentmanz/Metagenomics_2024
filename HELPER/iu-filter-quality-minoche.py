#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2011, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.


# This script takes an .ini file and analyzes quality of reads based on the recipe given in
# this paper:
#
#         http://genomebiology.com/2011/12/11/R112
#
# This is how the expected .ini file looks like:
#
#-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<------------
#    [general]
#    project_name = (project name)
#    researcher_email = (your e-mail address)
#    input_directory = (directory from which the input files will be read from)
#    output_directory = (directory to store output)
#
#
#    [files]
#    pair_1 = (pair 1 files, comma separated)
#    pair_2 = (pair 2 files, comma separated. must be ordered based on pair 1 files)
#-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<-------8<------------
#

import os
import sys
import configparser

E = os.path.exists
J = os.path.join

import IlluminaUtils.lib.fastqlib as u
from IlluminaUtils.utils.runconfiguration import RunConfiguration
from IlluminaUtils.utils.runconfiguration import RunConfigError
from IlluminaUtils.utils.helperfunctions import ReadIDTracker
from IlluminaUtils.utils.helperfunctions import QualityScoresHandler
from IlluminaUtils.utils.helperfunctions import visualize_qual_stats_dict


def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def IsHighQuality(sequence, quality, p):
    # boolean list of 'base quality is higher than q':
    bq = [True if _q > 2 else False for _q in quality]
    trim_to = None

    # check for B-tail
    for trim_to in range(len(bq) - 1, 0, -1):
        if bq[trim_to]:
            trim_to += 1
            break

    trim_to = None if len(sequence) == trim_to else trim_to

    if trim_to is not None:
        if trim_to * 1.0 / len(sequence) < p:
            # Exclude Sequence
            return (False, trim_to, 'FAILED_REASON_P')

        sequence = sequence[0:trim_to]
        quality = quality[0:trim_to]

    if sequence.count('N'):
        # Exclude Sequence
        return (False, trim_to, 'FAILED_REASON_N')

    # C33
    half_length = len(quality) // 2
    if len([True for _q in quality[0:half_length] if _q >= 30]) < 0.66 * half_length:
        # Exclude Sequence
        return (False, trim_to, 'FAILED_REASON_C33')

    # Accept sequence
    return (True, trim_to, 'PASSED')


def main(config, args = None, visualize_quality_curves = False, ignore_deflines = False, store_read_fate=False, p = 0.75):

    def get_percent(x, y):
        if y == 0:
            return 0.0
        else:
            return x * 100.0 / y

    if args:
        p = args.p
        ignore_deflines = args.ignore_deflines
        visualize_quality_curves = args.visualize_quality_curves
        limit_num_pairs = args.limit_num_pairs
        print_qual_scores = args.print_qual_scores
        store_read_fate = args.store_read_fate


    #####################################################################################
    # dealing with output file pointers..
    #####################################################################################

    GetFilePath = lambda p: os.path.join(config.output_directory, config.project_name + '-' + p)

    errors_fp                 = open(GetFilePath('STATS.txt'), 'w')
    quality_passed_r1         = u.FastQOutput(GetFilePath('QUALITY_PASSED_R1.fastq'))
    quality_passed_r2         = u.FastQOutput(GetFilePath('QUALITY_PASSED_R2.fastq'))
    id_tracker_output_path    = GetFilePath('READ_IDs.cPickle.z')

    if visualize_quality_curves:
        qual_dict_output_path = GetFilePath('Q_DICT.cPickle.z')


    #####################################################################################
    # some useful variables before we begin..
    #####################################################################################

    number_of_pairs = 0
    total_pairs_failed = 0
    total_pairs_passed = 0
    both_pairs_failed_qual = 0
    only_pair_1_failed_qual = 0
    only_pair_2_failed_qual = 0
    total_pair_1s_trimmed = 0
    total_pair_2s_trimmed = 0

    if visualize_quality_curves:
        qual_dict = QualityScoresHandler()

    if store_read_fate:
        read_id_tracker = ReadIDTracker()

    #####################################################################################
    # main loop per file listed in config:
    #####################################################################################
    COMPRESSED = lambda x: os.path.basename(config.pair_1[index]).endswith('.gz')

    for index in range(0, len(config.pair_1)):
        try:
            pair_1 = u.FastQSource(config.pair_1[index], lazy_init=True)
            pair_2 = u.FastQSource(config.pair_2[index], lazy_init=True)

        except u.FastQLibError as e:
            print("FastQLib is not happy.\n\n\t", e, "\n")
            sys.exit()


        #####################################################################################
        # main loop per read:
        #####################################################################################

        while pair_1.next(raw = ignore_deflines) and pair_2.next(raw = ignore_deflines):
            number_of_pairs += 1

            if limit_num_pairs and number_of_pairs == limit_num_pairs:
                break

            s1, q1 = pair_1.entry.sequence, pair_1.entry.process_Q_list()
            s2, q2 = pair_2.entry.sequence, pair_2.entry.process_Q_list()

            if pair_1.p_available:
                pair_1.print_percentage()

            #####################################################################################
            # compare reads
            #####################################################################################

            p1_passed_qual, p1_trim_to, p1_fate = IsHighQuality(s1, q1, p)
            p2_passed_qual, p2_trim_to, p2_fate = IsHighQuality(s2, q2, p)

            if (not p1_passed_qual) or (not p2_passed_qual):
                # either of these pairs failed quality check.
                total_pairs_failed += 1
                if (not p1_passed_qual) and p2_passed_qual:
                    only_pair_1_failed_qual += 1
                    fate = p1_fate
                elif p1_passed_qual and (not p2_passed_qual):
                    only_pair_2_failed_qual += 1

                    if store_read_fate:
                        read_id_tracker.update(pair_1, pair_2, p2_fate)

                    fate = p2_fate
                else:
                    both_pairs_failed_qual += 1
                    fate = p1_fate

                if store_read_fate:
                    read_id_tracker.update(pair_1, pair_2, fate)

                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, fate)

            else:
                total_pairs_passed += 1

                if p1_trim_to or p2_trim_to:
                    if p1_trim_to:
                        pair_1.entry.trim(trim_to = p1_trim_to)
                        total_pair_1s_trimmed += 1
                    if p2_trim_to:
                        pair_2.entry.trim(trim_to = p2_trim_to)
                        total_pair_2s_trimmed += 1

                if store_read_fate:
                    read_id_tracker.update(pair_1, pair_2, p1_fate)

                if visualize_quality_curves:
                    qual_dict.update(pair_1, pair_2, p1_fate)

                quality_passed_r1.store_entry(pair_1.entry)
                quality_passed_r2.store_entry(pair_2.entry)

            if print_qual_scores:
                print('')
                print('%s' % pair_1.entry.header_line)
                print('READ 1 [%s]: %s' % (p1_fate, ', '.join([str(_) for _ in q1])))
                print('READ 2 [%s]: %s' % (p2_fate, ', '.join([str(_) for _ in q2])))
                print('')
                print(' --')


    errors_fp.write('number of pairs analyzed      : %d\n' % (number_of_pairs))
    errors_fp.write('total pairs passed            : %d (%%%.2f of all pairs)\n' % (total_pairs_passed, get_percent(total_pairs_passed, number_of_pairs)))
    errors_fp.write('  total pair_1 trimmed        : %d (%%%.2f of all passed pairs)\n' % (total_pair_1s_trimmed, get_percent(total_pair_1s_trimmed, total_pairs_passed)))
    errors_fp.write('  total pair_2 trimmed        : %d (%%%.2f of all passed pairs)\n' % (total_pair_2s_trimmed, get_percent(total_pair_2s_trimmed, total_pairs_passed)))
    errors_fp.write('total pairs failed            : %d (%%%.2f of all pairs)\n' % (total_pairs_failed, get_percent(total_pairs_failed, number_of_pairs)))
    errors_fp.write('  pairs failed due to pair_1  : %d (%%%.2f of all failed pairs)\n' % (only_pair_1_failed_qual, get_percent(only_pair_1_failed_qual, total_pairs_failed)))
    errors_fp.write('  pairs failed due to pair_2  : %d (%%%.2f of all failed pairs)\n' % (only_pair_2_failed_qual, get_percent(only_pair_2_failed_qual, total_pairs_failed)))
    errors_fp.write('  pairs failed due to both    : %d (%%%.2f of all failed pairs)\n' % (both_pairs_failed_qual, get_percent(both_pairs_failed_qual, total_pairs_failed)))
    if store_read_fate:
        for fate in [f for f in read_id_tracker.fates if f.startswith('FAILED')]:
            errors_fp.write('  %s%s: %d (%%%.2f of all failed pairs)\n' % (fate,
                                                                           (' ' * (28 - len(fate))),
                                                                           len(read_id_tracker.ids[fate]),
                                                                           get_percent(len(read_id_tracker.ids[fate]), total_pairs_failed)))
    print()

    if store_read_fate:
        sys.stderr.write('\rRead ID tracker dict is being stored ...')
        sys.stderr.flush()
        read_id_tracker.store(id_tracker_output_path)
        sys.stderr.write('\n')

    if visualize_quality_curves:
        sys.stderr.write('\rQuality scores dict is being stored ...')
        sys.stderr.flush()
        qual_dict.store_dict(qual_dict_output_path)

        for entry_type in qual_dict.entry_types:
            sys.stderr.write('\rQuality scores visualization in progress: %s%s' % (entry_type, ' ' * 20))
            sys.stderr.flush()
            visualize_qual_stats_dict(qual_dict.data[entry_type], GetFilePath(entry_type),\
                            title = 'Mean PHRED scores for pairs tagged as "%s"' % entry_type)
        sys.stderr.write('\n')

    errors_fp.close()
    quality_passed_r1.close()
    quality_passed_r2.close()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Implementation of "http://genomebiology.com/content/12/11/R112"')
    parser.add_argument('user_config', metavar = 'CONFIG_FILE',
                                        help = 'User configuration to run. See the source code to\
                                                see an example.')
    parser.add_argument('-p', default = 0.75, type=float, metavar = 'FLOAT',
                                        help = 'Minimum high-quality read length (default: %(default)s)')
    parser.add_argument('--ignore-deflines', action = 'store_true', default = False,
                                        help = 'If FASTQ files are not CASAVA outputs, parsing the header info\
                                                may go wrong. This flag tells the software to skip parsing\
                                                deflines.')
    parser.add_argument('--visualize-quality-curves', action = 'store_true', default = False,
                                        help = 'When set, mean quality score for individual bases will be\
                                                stored and visualized for each group of reads.')
    parser.add_argument('--limit-num-pairs', default = None, type=int, metavar = 'INTEGER',
                                        help = 'Put a limit to the number of pairs to analyze. For testing purposes.')
    parser.add_argument('--print-qual-scores', action = 'store_true', default = False,
                                        help = 'When set, the script will print out the Q-scores the way it sees it in\
                                                the FASTQ file. This flag will generate a lot of useless output to the\
                                                stdout, and you should not use it if you are not testing something.')
    parser.add_argument('--store-read-fate', action = 'store_true', default = False,
                                        help = 'As it goes through your raw reads, this program keeps track of the read fate\
                                                so you can learn what happened to a given read ID in your raw input data once\
                                                the analysis is done. This output can become extremely large, and often is\
                                                utterly useless to you unless you have a very specific benchmarking or debugging\
                                                interestes, hence, it is not stored by default. You can change that behavior by\
                                                using this flag, and ask illumina-utils to store this data on your disk.')

    args = parser.parse_args()
    user_config = configparser.ConfigParser()
    user_config.read(args.user_config)

    try:
        config = RunConfiguration(user_config)
    except RunConfigError as e:
        print(e)
        sys.exit()


    try:
        sys.exit(main(config, args))
    except RunConfigError as e:
        print(e)
        sys.exit(-1)
    except u.FastQLibError as e:
        print(e)
        sys.exit(-2)