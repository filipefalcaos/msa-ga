import math
import random
import nwalign3 as nw


# Prints a chromosome's sequences
def print_chromosome(chromosome):
    for seq in chromosome:
        print(seq)


# Calc the initial maximum columns m
def calc_m(lines):
    m_aux = 0
    lengths = []

    for line in lines:
        length = len(line)
        lengths.append(length)

        if length >= m_aux:
            m_aux = length

    return m_aux


# Add the required gaps
def add_gaps(lines):

    # Complete with gaps
    for i in range(len(lines)):
        diff = m_1 - len(lines[i])
        for j in range(diff):
            lines[i] += "-"

    # Add the scale factor
    k = 1.2
    m = math.ceil(m_1 * k)

    # Add the scale factor columns
    diff = m - m_1
    for i in range(diff):
        for j in range(len(lines)):
            lines[j] += "-"

    # Add the offset columns
    offset = math.ceil(m_1 * 0.2)
    for i in range(offset):
        for j in range(len(lines)):
            lines[j] = "-" + lines[j]

    return lines


# Read input file string
with open("../data/data_1.txt") as f:
    input_str = f.read()

# Get lines
input_str = input_str.replace(" ", "")
lines_list = input_str.split("\n")

# Calc m1
m_1 = calc_m(lines_list)

# Add the gaps
lines_list = add_gaps(lines_list)

# Prints the original chromosome
print("Input matrix:")
print_chromosome(lines_list)

# Create the initial population
C = 100
pop = []
for c in range(C):
    # New chromosome
    lines_list_aux = []

    # Use nwalign to compute the pairwise alignments
    # by the Needleman-Wunsch algorithm
    for i in range(len(lines_list)):
        alignments = []

        # Compute the pairwise alignments
        for j in range(len(lines_list)):
            if i != j:
                curr_alignment = nw.global_align(lines_list[i], lines_list[j])
                alignments.append(curr_alignment)

        # Randomly select an alignment
        alignment = random.choice(alignments)
        alignment = alignment[0]
        lines_list_aux.append(alignment)

    # Add the chromosome generated and prints it
    pop.append(lines_list_aux)
    print("\nChromosome " + str(c + 1) + ":")
    print_chromosome(lines_list_aux)
