import os
import gzip
from collections import Counter
from multiprocessing import Pool, cpu_count

from Bio import SeqIO
from tqdm import tqdm
import pandas as pd


def process_file(args):
    file, target_set = args
    sample = file.replace('.fastq.gz', '').replace('.fastq', '')
    total_reads = 0
    correct_sequences = 0
    counts = Counter()

    open_func = gzip.open if file.endswith('.gz') else open
    with open_func(file, 'rb') as handle:
        readline = handle.readline
        target_local = target_set
        counts_local = counts

        while True:
            header = readline()
            if not header:
                break
            seq = readline().rstrip(b'\n\r')
            readline()  # plus
            readline()  # qual

            total_reads += 1
            seq20 = seq[:20]
            if seq20 in target_local:
                counts_local[seq20] += 1
                correct_sequences += 1

    return sample, total_reads, correct_sequences, counts


if __name__ == '__main__':
    # We first open the library file as it was synthesized to obtain the sequence of all guides.
    with open('_library.txt', 'r') as library_handle:
        library = list(SeqIO.parse(library_handle, 'fasta'))

    os.chdir('./FASTQ')
    files = [f for f in os.listdir(os.getcwd()) if f.endswith(('.fastq', '.fastq.gz'))]

    guides = [str(l.seq[:20]).encode() for l in library]
    target_set = set(guides)
    sample_names = [f.replace('.fastq.gz', '').replace('.fastq', '') for f in files]

    # We then create a dictionary to store count data:
    # It contains 1 key per guide, then each value is a dictionary with 1 key and value per sample.
    targets = {}
    for guide in guides:
        guide_str = guide.decode()
        targets[guide_str] = {}
        for sample in sample_names:
            targets[guide_str][sample] = 0

    # We then go through each fastq sequencing file (1 per condition and replicate).
    # If the read is in the library, then we store it in the dictionary.
    nproc = min(len(files), cpu_count())
    results = []
    with Pool(processes=nproc) as pool:
        for result in tqdm(
            pool.imap_unordered(process_file, [(f, target_set) for f in files]),
            total=len(files)
        ):
            results.append(result)

    with open("./1.1_count_stats.txt", "w") as stats_out:
        for sample, total_reads, correct_sequences, counts in results:
            for guide, count in counts.items():
                targets[guide.decode()][sample] = count

            stats_out.write(sample + ":\n")
            stats_out.write(str(total_reads) + " reads\n")
            stats_out.write(str(round(100 * correct_sequences / total_reads, 2)) + "% correct reads\n")

    # The counts table is finally saved as a dataframe.
    counts_table = pd.DataFrame.from_dict(targets, orient='index')
    os.chdir('./..')
    counts_table.to_csv('./1.0_count_table.txt', sep='\t')