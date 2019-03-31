import math
import random
import nwalign3 as nw

# Global vars
k = 1.2
C = 100
e = 1.5
N = 100


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


# Add spaces only to complete
def add_spaces(lines):
    m = calc_m(lines)
    for i in range(len(lines)):
        diff = m - len(lines[i])
        for j in range(diff):
            lines[i] += " "

    return lines


# Add the required gaps
def add_gaps(lines):
    # Complete with gaps
    for i in range(len(lines)):
        diff = m_1 - len(lines[i])
        for j in range(diff):
            lines[i] += "-"

    # Add the scale factor
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


# The g function defined in the paper
def g_func(x, y):
    if (x != y) and (x != "-") and (y != "-"):
        return 0
    elif (x != y) and (x == "-" or y == "-"):
        return 1
    elif x == y:
        return 2

    return 0


# Get the number of columns with only gaps
def get_gap_cols(lines):
    n_cols = 0
    only_gaps = True
    n = len(lines)

    for i in range(calc_m(lines)):
        for j in range(n):
            if (lines[j][i] != " ") and (lines[j][i] != "-"):
                only_gaps = False

        if only_gaps:
            n_cols += 1
        else:
            only_gaps = True

    return n_cols


# Apply the f function defined in the paper
def f_func(lines):
    value = 0
    lines = add_spaces(lines)
    x = get_gap_cols(lines)
    m = calc_m(lines)

    for i in range(m):
        for j in range(len(lines) - 1):
            value += g_func(lines[j][i], lines[j + 1][i])

    value = value - 2 * x * (len(lines) - 1) - e * x
    return value / m


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

count = 0
while count < 1:

    # Evaluate all chromosomes
    evals = []
    for p in pop:
        evals.append(f_func(p))

    # New population
    new_pop = []

    # Normalize the fitness
    pop_sum = sum(evals)
    evals_aux = []
    for ev in evals:
        evals_aux.append(ev / pop_sum)

    # Build the mating pool
    pool = []
    for i in range(len(pop)):
        prob = int(evals_aux[i] * 100)

        for j in range(prob):
            pool.append(pop[i])

    # Randomly select two parents
    p1 = random.choice(pool)
    p2 = random.choice(pool)
    print("\nParent 1: ")
    print_chromosome(p1)
    print("\nParent 2: ")
    print_chromosome(p2)

    # Select crossover method
    method = random.randint(0, 1)
    child = p1
    m = calc_m(pop[0])
    n = len(pop[0])

    # Apply vertical crossover
    if method == 0:
        rand_v = random.randint(1, m - 1)
        print("\nSlice vertically by: " + str(rand_v))

        for i in range(len(p1)):
            print(p1[i][0:rand_v] + "|", end="")
            print(p2[i][rand_v:])
            child[i] = p1[i][0:rand_v] + p2[i][rand_v:]

    # Apply horizontal crossover
    else:
        rand_h = random.randint(1, n - 1)
        print("\nSlice horizontally by: " + str(rand_h))
        print_chromosome(p1[:rand_h])
        [print("*", end="") for i in range(m)]
        print()
        print_chromosome(p2[rand_h:])
        child = p1[:rand_h] + p2[rand_h:]

    count = count + 1
