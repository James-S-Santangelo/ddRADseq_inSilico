from Bio import SeqIO
from Bio import Restriction as rst
import os


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
        RE_1 (:obj: 'Restriction'); Restriction enzyme imported from BioPython's Restriction module
        RE_2 (:obj: 'Restriction'); Restriction enzyme imported from BioPython's Restriction module

    Returns:
        split_sequences (Dict of str: :obj:'dict' of str: :obj:'list' of obj:'int'): A complicated nested dictionary.
        Primary keys are sequence record ID's. Primary values are dictionaries where keys are strings representing the
        REs used in the analysis and the values are lists with values corresponding to the cut positions of the enzyme
        along the DNA sequences. Additional values are lists with integers representing the length of fragments bounded
        the same enzyme (2 keys: one for each RE) or bounded on either end by a different enzyme (1 key).
    """

    # Retrieve RE's
    RE_1 = getattr(rst, RE_1)
    RE_2 = getattr(rst, RE_2)

    # Create restriction batch object
    RE_batch = RE_1 + RE_2

    # Initialize dictionary
    split_sequences = {}

    # Create strings from REs that will later become dictionary keys
    RE_1_both_sides_key = "{0}_2".format(RE_1)
    RE_2_both_sides_key = "{0}_2".format(RE_2)
    RE_2_RE_1_key = "{0}_{1}".format(RE_2, RE_1)

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
        split_sequences[record.id] = rst.Analysis(RE_batch, record.seq, linear=True)

    # Iterate over dictionary and for each record ID, extract the cut positions for
    # both enzymes and create a list that merges these and sorts them in ascending order.
    # Add merged list as additional key:value pair to dictionary
    for ID, ana in split_sequences.items():
        ana.full()["merged"] = sorted(ana.full()[RE_1] + ana.full()[RE_2])
        split_sequences[ID] = ana.full()

    # Iterate over sequence records
    for record in sequence_records:

        sequence = record.seq # Get sequence

        # Extract list with cut positions for both enzymes. # None added to cut_positions
        # to make sure last fragment is included
        cut_positions = split_sequences[record.id]["merged"] + [None]

        # Create lists in dictionary to strore fragment lengths.
        split_sequences[record.id][RE_1_both_sides_key] = []
        split_sequences[record.id][RE_2_both_sides_key] = []
        split_sequences[record.id][RE_2_RE_1_key] = []
        split_sequences[record.id]["Singles"] = []

        # Create iterator based on number of cut positions
        for i in range(len(cut_positions) - 1):

            # If not the last fragments along the sequence
            if cut_positions[i + 1] != None:

                # Generate fragment as the sequence between two consecutive cutpositions
                fragment = str(sequence[cut_positions[i] - 1:cut_positions[i + 1] - 1])

                # If fragment is bounded on both sides by RE_1 sequence
                if fragment[0:RE_1_overhang_length] == RE_1.site[RE_1_FivePrime_cut_index:] and fragment[-1] == RE_1.site[:RE_1_FivePrime_cut_index]:
                    split_sequences[record.id][RE_1_both_sides_key].append(len(fragment))

                # If fragment is bounded on both sides by RE_2 sequence
                elif fragment[0:RE_2_overhang_length] == RE_2.site[RE_2_FivePrime_cut_index:] and fragment[-1] == RE_2.site[:RE_2_FivePrime_cut_index]:
                    split_sequences[record.id][RE_2_both_sides_key].append(len(fragment))

                # If fragment is bounded on left (i.e. 5 prime) by RE_1 sequence and right (i.e. 3 prime) by RE_2 sequence
                elif fragment[0:RE_1_overhang_length] == RE_1.site[RE_1_FivePrime_cut_index:] and fragment[-1] == RE_2.site[:RE_2_FivePrime_cut_index]:
                    split_sequences[record.id][RE_2_RE_1_key].append(len(fragment))

                # If fragment is bounded on left by RE_2 sequence and right by RE_1 sequence
                elif fragment[0:RE_2_overhang_length] == RE_2.site[RE_2_FivePrime_cut_index:] and fragment[-1] == RE_1.site[:RE_1_FivePrime_cut_index]:
                    split_sequences[record.id][RE_2_RE_1_key].append(len(fragment))

                # No cut site on one end (e.g. ends of sequences)
                else:
                    split_sequences[record.id]["Singles"].append(len(fragment))

            # If it is the last fragment.
            else:
                fragment = str(sequence[cut_positions[i] - 1:cut_positions[i + 1]])
                split_sequences[record.id]["Singles"].append(len(fragment))

    return split_sequences


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
    enz_1 = args['RE_1']
    enz_2 = args['RE_2']

    # Retrieve multi-fasta file
    multi_fasta = [f for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
    if len(multi_fasta) == 1:
        multi_fasta = fasta_dir + multi_fasta[0]
    else:
        sys.exit("There is more than one fasta file in the directory. Exiting!")

    # Load multi-fasta file
    sequence_records = list(SeqIO.parse(multi_fasta, "fasta"))

    print "Performing dual digest of multi-fasta file now"
    ddRAD_digest(sequence_records, RE_1, RE_2)
