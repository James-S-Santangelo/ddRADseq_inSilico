from Bio import SeqIO
from Bio import Restriction as rst
from matplotlib import pyplot as plt
import pandas as pd
from statistics import mean, median
import os
import time
import argparse


def create_directory(directory):
    """Creates directory if it does not exist

    Args:
        directory ('str'): Directory to be created

    Returns:
        None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)


def ddRAD_digest(sequence_records, RE_1, RE_2):
    """Performs a double digest of DNA sequences

    Args:
        sequence_records (:obj:'list' of :obj:'SeqRecord'): List with SeqRecords imported using BioPython's SeqIO module
        RE_1 (str): Name of first restriction enzyme
        RE_2 (str): Name of second restriction enzyme

    Returns:
        fragment_lengths (Dict of str: :obj:'dict' of str: :obj:'list' of obj:'int'): A complicated nested dictionary.
        Primary keys are sequence record ID's. Primary values are dictionaries where keys are strings representing the
        REs used in the analysis and the values are lists with values corresponding to the cut positions of the enzyme
        along the DNA sequences. Additional values are lists with integers representing the length of fragments bounded
        the same enzyme (2 keys: one for each RE) or bounded on either end by a different enzyme (1 key).
    """

    # Counters
    fragment_counter = 0
    bound_by_RE_1 = 0
    bound_by_RE_2 = 0
    dual_fragment_counter = 0
    singletons = 0

    # Retrieve RE's
    RE_1 = getattr(rst, RE_1)
    RE_2 = getattr(rst, RE_2)

    # Create restriction batch object
    RE_batch = RE_1 + RE_2

    # Initialize dictionary
    cut_positions = {}
    fragment_lengths = {}

    # Create strings from REs that will later become dictionary keys
    RE_1_both_sides_key = "{0}_2".format(RE_1)
    RE_2_both_sides_key = "{0}_2".format(RE_2)
    dual_key = "{0}_{1}".format(RE_1, RE_2)

    # Identify overhangs produced by digest and determine index in recognition site
    # where the cut occurs. Used to determine which RE cut on either end of the fragments
    # so they can be placed in the appropriate bin (above)
    RE_1_FivePrime_cut_index = RE_1.elucidate().find('^')
    RE_2_FivePrime_cut_index = RE_2.elucidate().find('^')
    RE_1_overhang_length = len(RE_1.site[RE_1_FivePrime_cut_index:])
    RE_2_overhang_length = len(RE_2.site[RE_2_FivePrime_cut_index:])

    # Iterate over sequence records and perform a restriction analysis
    # Return an analysis object where the cut positions of all REs
    # can be extracted. Add Analysis object as value to dictionary where the
    # key is the sequence record ID.
    for record in sequence_records:
        cut_positions[record.id] = rst.Analysis(RE_batch, record.seq, linear=True)

    # Iterate over dictionary and for each record ID, extract the cut positions for
    # both enzymes and create a list that merges these and sorts them in ascending order.
    # Add merged list as additional key:value pair to dictionary
    for ID, ana in cut_positions.items():
        ana.full()["merged_cuts"] = sorted(ana.full()[RE_1] + ana.full()[RE_2])
        cut_positions[ID] = ana.full()

    # Iterate over sequence records
    for record in sequence_records:

        sequence = str(record.seq) # Get sequence

        # Extract list with cut positions for both enzymes. # None added to cut_positions
        # to make sure last fragment is included
        merged_cut_positions = [0] + cut_positions[record.id]["merged_cuts"] + [None]

        # Create lists in dictionary to strore fragment lengths.
        fragment_lengths[record.id] = {RE_1_both_sides_key: [],
                                       RE_2_both_sides_key: [],
                                       dual_key: [],
                                       "Singles": []}

        fragments = []
        for i, j in zip(merged_cut_positions, merged_cut_positions[1:]):
            # print i, j
            if i == 0 and j is not None:
                fragments.append(sequence[i:j - 1])
            elif j is None:
                fragments.append(sequence[i - 1:j])
            else:
                fragments.append(sequence[i - 1:j - 1])

        for fragment in fragments:
            fragment_counter += 1
#             print fragment

            # If fragment is bounded on both sides by RE_1 sequence
            if fragment[0:RE_1_overhang_length] == RE_1.site[RE_1_FivePrime_cut_index:] and fragment[-1] == RE_1.site[:RE_1_FivePrime_cut_index]:
#                 print "RE_1 on both sides"
                bound_by_RE_1 += 1
                fragment_lengths[record.id][RE_1_both_sides_key].append(len(fragment))

            # If fragment is bounded on both sides by RE_2 sequence
            elif fragment[0:RE_2_overhang_length] == RE_2.site[RE_2_FivePrime_cut_index:] and fragment[-1] == RE_2.site[:RE_2_FivePrime_cut_index]:
#                 print "RE_2 on both sides"
                bound_by_RE_2 += 1
                fragment_lengths[record.id][RE_2_both_sides_key].append(len(fragment))

            # If fragment is bounded on left (i.e. 5 prime) by RE_1 sequence and right (i.e. 3 prime) by RE_2 sequence
            elif fragment[0:RE_1_overhang_length] == RE_1.site[RE_1_FivePrime_cut_index:] and fragment[-1] == RE_2.site[:RE_2_FivePrime_cut_index]:
#                 print "Dual"
                dual_fragment_counter += 1
                fragment_lengths[record.id][dual_key].append(len(fragment))

            # If fragment is bounded on left by RE_2 sequence and right by RE_1 sequence
            elif fragment[0:RE_2_overhang_length] == RE_2.site[RE_2_FivePrime_cut_index:] and fragment[-1] == RE_1.site[:RE_1_FivePrime_cut_index]:
#                 print "Dual"
                dual_fragment_counter += 1
                fragment_lengths[record.id][dual_key].append(len(fragment))

            # No cut site on one end (e.g. ends of sequences)
            else:
#                 print "End piece"
                singletons += 1
                fragment_lengths[record.id]["Singles"].append(len(fragment))

    print "The digest generated a total of {0} fragments".format(fragment_counter)
    time.sleep(1)
    print "The digest generated {0} fragments bound on both sides by {1}".format(bound_by_RE_1, RE_1)
    time.sleep(1)
    print "The digest generated {0} fragments bound on both sides by {1}".format(bound_by_RE_2, RE_2)
    time.sleep(1)
    print "The digest generated {0} fragments bound by different cut sites on either side".format(dual_fragment_counter)
    time.sleep(1)
    print "The digest generated {0} fragments with a cut side on only one side".format(singletons)
    time.sleep(1)

    return fragment_lengths


def ddRAD_fragments(ddRAD_dict, RE_1, RE_2):
    """Following double digest, extract length of fragments pooled across
    all digested sequence records

    Args:
        ddRAD_dict (Dict of str: :obj:'dict' of str: :obj:'list' of obj:'int'): A complicated nested dictionary.
        Primary keys are sequence record ID's. Primary values are dictionaries where keys are strings representing
        the REs used in the analysis and the values are lists with values corresponding to the cut positions of the
        enzyme along the DNA sequences. Additional values are lists with integers representing the length of fragments
        bounded the same enzyme (2 keys: one for each RE) or bounded on either end by a different enzyme (1 key).
        RE_1 (str): Name of first restriction enzyme
        RE_2 (str): Name of second restriction enzyme

    Returns:
        Frag_list: Sorted list of integers representing the length of each fragment that is bound on both ends
        by cuts from a different RE. Pooled across all digested sequence records.
    """

    # Initialize list
    Frag_list = []

    # Key for extraxting required values from nested dictionary
    dual_key = "{0}_{1}".format(RE_1, RE_2)

    # Iterate over nested dictionary
    for sequence, cut_dict in ddRAD_dict.items():
        # Interate fragment length integers in the list composed of fragment lengths for fragments bounded by two different REs
        for fragment_length in ddRAD_dict[sequence][dual_key]:
            Frag_list.append(fragment_length)

    return sorted(Frag_list)


def summary_statistics(fragment_lengths, min_bin_size, max_bin_size, step, RE_1, RE_2, outfile):
    """Generate summary statistics for fragment sizes generated
    from dual digest

    Args:
        fragment_lengths (:obj: of 'list'): List with length of fragments bound
        by different cut sites on either side.
        min_size (int): Minimum fragment size to include in output
        max_size (int): Maximum fragment size to include in output
        step (int): Size of bins (in bp)
        RE_1 (str): Name of first restriction enzyme
        RE_2 (str): Name of second restriction enzyme

    Returns:
        None: Writes output to disk
    """

    #Basic summary stats
    num_fragments = len(fragment_lengths)
    min_fragment_size = min(fragment_lengths)
    max_fragment_size = max(fragment_lengths)
    mean_size = mean(fragment_lengths)
    median_size = median(fragment_lengths)

    series = pd.Series(fragment_lengths)
    subset = series[series.between(min_bin_size, max_bin_size, inclusive=True)]
    ranges = [i for i in range(min_bin_size, max_bin_size + step, step)]
    range_count = subset.groupby(pd.cut(subset, ranges)).count()
    total_in_range = sum(range_count[i] for i in range(len(ranges) - 1))
    total_length_in_range = sum(fragment_lengths[i] for i in range(len(fragment_lengths) - 1) if fragment_lengths[i] >= min_bin_size and fragment_lengths[i] <= max_bin_size)

    with open(outfile, "w") as f:
        f.write("The digest with {0} annd {1} generated a total of {2} fragments\n".format(RE_1, RE_2, num_fragments))
        f.write("The smallest fragments generated by the digest was {0} bp\n".format(min_fragment_size))
        f.write("The largest fragments generated by the digest was {0} bp\n".format(max_fragment_size))
        f.write("The mean size of fragments generated from the digest was {0} bp\n".format(mean_size))
        f.write("The median size of fragments generated from the digest was {0} bp\n".format(median_size))
        f.write("The digest generated the following number of fragments within the specified ranges:\n")
        f.write("{0}\n".format(range_count))
        f.write("The digest generated a total of {0} fragments within {1} to {2} bp, totalling {3} bp in length".format(total_in_range, min_bin_size, max_bin_size, total_length_in_range))


if __name__ == "__main__":

        # Create directories
    print("Creating output directories")
    plot_dir = "./plots/"
    output_dir = "./output/"
    create_directory(plot_dir)
    create_directory(output_dir)
    time.sleep(1)

    # Define command line arguments
    parser = argparse.ArgumentParser(prog="From a multi-fasta file (e.g. genome sequence, perform and in-silico dual digest with specified restriction enzymes, returning the distribution of fragment lengths bound by cuts sites both enzymes.")
    parser.add_argument('fasta_dir', help="Path to directory containing multi-fasta file", type=str)
    parser.add_argument('RE_1', help="Name of first restriction enzyme used in digest", type=str)
    parser.add_argument('RE_2', help="Name of second restriction enzyme used in digest", type=str)
    args = vars(parser.parse_args())

    # Retrieve arguments passed from command line
    fasta_dir = args['fasta_dir']
    RE_1 = args['RE_1']
    RE_2 = args['RE_2']

    # Retrieve multi-fasta file
    multi_fasta = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta') or f.endswith('.fa')]
    if len(multi_fasta) == 1:
        multi_fasta = fasta_dir + multi_fasta[0]
    else:
        sys.exit("There is more than one fasta file in the directory. Exiting!")

    # Load multi-fasta file
    sequence_records = list(SeqIO.parse(multi_fasta, "fasta"))

    print "Performing dual digest of multi-fasta file now"
    fragment_lengths = ddRAD_digest(sequence_records, RE_1, RE_2)

    time.sleep(3)
    print "Retrieving lengths of fragments bound by different cut sites on either side"
    dual_frag_lengths = ddRAD_fragments(fragment_lengths, RE_1, RE_2)

    print "Creating histogram showing distribution of fragment lengths"
    time.sleep(1)

    plt.hist(dual_frag_lengths, bins=100)
    plt.title('{0} and {1} Digest'.format(RE_1, RE_2), fontsize=20)
    plt.xlabel("Fragment size", fontsize=18)
    plt.ylabel('Count', fontsize=18)

    plot_name = plot_dir + '{0}-{1}_digest.png'.format(RE_1, RE_2)
    plt.savefig(plot_name, dpi=100)

    print "The histogram is located in {0}".format(plot_name)
    time.sleep(1)

    print "Creating output with summary statistics now!"

    print ""
    print "PLEASE ANSWER TO THE QUESTIONS BELOW"

    min_bin_size = int(raw_input("What is the minimum fragment size you would like included in the summary statistics? "))
    max_bin_size = int(raw_input("What is the maximum fragment size you would like included in the summary statistics? "))
    step = int(raw_input("What size bins (in bp) would you like to group fragments into for the summary statistics? "))

    time.sleep(1)

    outfile = output_dir + "{0}-{1}_digest_summaryStats.txt"
    summary_statistics(dual_frag_lengths, min_bin_size, max_bin_size, step, RE_1, RE_2, outfile)


