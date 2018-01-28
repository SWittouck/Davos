#!/usr/bin/python3

import argparse
import os

# class for reads
class Read:

    def __init__(self, sequence, index, count):
        self.sequence = sequence
        self.index = index
        self.count = count
        self.cost = 0
        self.read_prev = None

# class for dna files
class DnaFile:

    def __init__(self, name, reads_nested):
        self.name = name
        self.reads_nested = reads_nested

    # take subset of reads per position
    def subset_reads_per_index(self, n):
        for reads in self.reads_nested:
            reads.sort(key = lambda read: read.count, reverse = True)
            del reads[n:]

    # dynamic programming: update cost and read_prev for a given read
    def _update_read(self, read):
        reads_nested = self.reads_nested
        candidate_costs = []
        sequence = read.sequence
        for read_prev in reads_nested[read.index - 1]:
            sequence_prev = read_prev.sequence
            cost_prev = read_prev.cost
            cost_candidate = cost_prev + compare_seqs(sequence, sequence_prev, 25)
            try:
                read_prev_2 = read_prev.read_prev
                sequence_prev_2 = read_prev_2.sequence
                cost_candidate += compare_seqs(sequence, sequence_prev_2, 50)
                read_prev_3 = read_prev_2.read_prev
                sequence_prev_3 = read_prev_3.sequence
                cost_candidate += compare_seqs(sequence, sequence_prev_3, 75)
            except AttributeError:
                pass
            candidate_costs.append(cost_candidate)
        cost = min(candidate_costs)
        read_prev_ix = candidate_costs.index(cost)
        read.cost = cost
        read.read_prev = reads_nested[read.index - 1][read_prev_ix]

    # dynamic programming: update costs and read_prevs for all reads
    def _update_reads(self):
        reads_nested = self.reads_nested
        if len(reads_nested[0]) == 0:
            reads_nested[0].append(Read('C' * 100, 0, 1))
        for index in range(1, len(reads_nested)):
            if len(reads_nested[index]) == 0:
                reads_nested[index].append(Read('C' * 100, index, 1))
            for read in reads_nested[index]:
                self._update_read(read)

    # dynamic programming: perform backtracking
    def _perform_backtracking(self):
        reads_nested = self.reads_nested
        read_current = min(reads_nested[-1], key = lambda read: read.cost)
        for index in range(len(reads_nested) - 1, 0, -1):
            reads_nested[index] = [read_current]
            read_current = read_current.read_prev
        reads_nested[0] = [read_current]

    # perform dynamic programming
    def perform_dynamic_programming(self):
        self._update_reads()
        self._perform_backtracking()

    # assemble a nested list of reads using a consensus
    # per position (cpp) strategy
    def assemble_reads_cpp(self):

        reads_nested = self.reads_nested

        # make empty assembly (list of 'N' characters)
        assembly = ['N'] * (len(reads_nested) * 25 + 75)

        # fill in empty assembly with consensus characters 
        for position in range(0, len(assembly)):
            base_counts = {'A':0, 'T':0, 'C':0, 'G':0}
            index_ix_read_4 = position // 25
            position_in_read_4 = position % 25
            for k in range(0, 4):
                index_ix = index_ix_read_4 - k
                position_in_read = position_in_read_4 + k * 25
                try:
                    if index_ix < 0:
                        raise IndexError
                    for read in reads_nested[index_ix]:
                        base = read.sequence[position_in_read]
                        count = read.count
                        base_counts[base] += count
                except IndexError:
                    pass
            try:
                consensus_char = max(base_counts, key = base_counts.get)
            except IndexError:
                consensus_char = 'C'
            assembly[position] = consensus_char

        # return the assembly
        return(assembly)

# compare two sequences
def compare_seqs(sequence_1, sequence_2, offset):
    n = len(sequence_1) - offset
    cost = 0
    for i in range(n):
        char_1 = sequence_1[i] 
        char_2 = sequence_2[i + offset]
        if char_1 != char_2:
            cost += 1
    return(cost)

# parse tsv file into dictionary where the keys are the dna file numbers
# and the values are lists containing all reads of that dna file
def parse_tsv(fin):
    hin = open(fin, "r")
    next(hin) # skip header
    dna_files = {}
    for line in hin:
        line = line.strip().split("\t")
        dna_file_number = line[0]
        read = Read(line[3], int(line[1]), int(line[2]))
        dna_files.setdefault(dna_file_number, []).append(read)
    hin.close()
    return(dna_files)

# determine the number of indices in a dna file given a list of reads
# when no reads have been observed for 10 or more consecutive indices,
# higher indices will be dropped
def determine_n_indices(reads):
    indices = [read.index for read in reads]
    indices.sort()
    n_indices = 0
    index_prev = 0
    for index in indices:
        if (index - index_prev > 10):
            n_indices = index_prev
            break
        index_prev = index
    return(n_indices)
   
# convert list of reads to nested list (sublist per index)
def nest_reads(reads, n_indices):
    reads_nested = [list([]) for _ in range(0, n_indices)]
    for read in reads:
        try:
            reads_nested[read.index].append(read)
        except IndexError:
            pass
    return(reads_nested)

# replace homodinucleotides
def replace_homodinucleotides(assembly):

    # replace homodinucleotides in positions O - (l-3)
    bases = {'A', 'C', 'T', 'G'}
    for position in range(len(assembly) - 2):
        if assembly[position] == assembly[position + 1]:
            base_repl = bases - {assembly[position], assembly[position + 2]}
            base_repl = list(base_repl)[0]
            assembly[position + 1] = base_repl

    # replace possible final homodinucleotide
    if assembly[len(assembly) - 2] == assembly[len(assembly) - 1]:
        base_repl = bases - {assembly[len(assembly) - 2]}
        base_repl = list(base_repl)[0]
        assembly[len(assembly) - 1] = base_repl

    return(assembly)

# parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--max_indices_per_file",
        help = "set a maximum number of indices that one file can have",
        type = int, default = "1500"
    )
    parser.add_argument(
        "-n", "--max_sequences_per_index",
        help = "take the top n most frequent sequences per index",
        type = int
    )
    parser.add_argument(
        "-d", "--dynamic_programming",
        help = "identify the optimal sequence per index using dynamic programming",
        action = "store_true"
    )
    parser.add_argument(
        "-i", "--fin_sequences",
        help = "path to the tvs file with sequences"
    )
    parser.add_argument(
        "-o", "--dout_assemblies",
        help = "path to output directory to store assemblies"
    )
    args = parser.parse_args()
    return(args) 

if __name__ == "__main__":

    args = parse_arguments()

    max_indices = args.max_indices_per_file
    max_sequences = args.max_sequences_per_index
    dynamic_programming = args.dynamic_programming
    fin = args.fin_sequences
    dout = args.dout_assemblies

    if not os.path.exists(dout):
        print("creating output directory")
        os.makedirs(dout, exist_ok = True)

    print("parsing tsv file with sequences")
    dna_files = parse_tsv(fin)

    for dna_file_number, dna_file in dna_files.items():
        print("file number %s" % (dna_file_number))
        n_indices = determine_n_indices(dna_file)
        print("number of indices found: %i" % (n_indices))
        if n_indices > max_indices:
            print("reducing number of indices to max of %i" % (max_indices))
            n_indices = max_indices
        reads_nested = nest_reads(dna_file, n_indices)
        dna_file = DnaFile(dna_file_number, reads_nested)
        if not max_sequences is None:
            print("taking the %i most abundant sequences per index" % (max_sequences))
            dna_file.subset_reads_per_index(max_sequences)
        if dynamic_programming:
            print("performing dynamic programming")
            dna_file.perform_dynamic_programming()
        print("assembling reads")
        assembly = dna_file.assemble_reads_cpp()
        print("replacing homodinucleotides")
        assembly = replace_homodinucleotides(assembly)
        assembly = "".join(assembly)
        # with open(dout + "/assembly_cpp_file_" + str(dna_file_number) + ".dna", "w") as file:
        with open("%s/assembly_file_%s.dna" % (dout, dna_file_number), "w") as file:
            file.write(assembly)
