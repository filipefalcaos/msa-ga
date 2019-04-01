import math
import random
import copy
import nwalign3 as nw

# Global vars
k = 1.2
C = 100
e = 1.5
N = 300
mutation_rate = 0.03


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


# Remove spaces to standardize
def remove_spaces(pop):
    for chromosome in pop:
        for i in range(len(chromosome)):
            chromosome[i] = chromosome[i].replace(" ", "")

    return pop


# Remove useless gaps from an alignment
def remove_useless_gaps(lines):
    to_rm = []
    only_gaps = True
    m = calc_m(lines)

    # Find columns with gaps only
    for i in range(m):
        for j in range(n):
            if (lines[j][i] != " ") and (lines[j][i] != "-"):
                only_gaps = False

        if only_gaps:
            to_rm.append(i)
        else:
            only_gaps = True

    # Remove gap-only columns
    to_rm.reverse()
    for i in range(len(lines)):
        line = []
        for j in range(m):
            line.append(lines[i][j])

        for j in to_rm:
            line.pop(j)

        lines[i] = ''.join(line)

    # Return the new matrix
    return lines


# Add the required gaps
def add_gaps(lines):

    # Complete with gaps
    for i in range(len(lines)):
        diff = calc_m(lines) - len(lines[i])
        for j in range(diff):
            lines[i] += "-"

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
    return value


# Read input file string
with open("../data/data_1.txt") as f:
    input_str = f.read()

# Get lines
input_str = input_str.replace(" ", "")
lines_list = input_str.split("\n")

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
    lines_list_aux = add_gaps(lines_list_aux)
    pop.append(lines_list_aux)
    print("\nChromosome " + str(c + 1) + ":")
    print_chromosome(lines_list_aux)

# Repeat N times or until a good solution
# appears
best_chromosome = None
count = 0
new_pop = []
print()

while count < N:

    # Evaluate all chromosomes
    evals = []
    for p in pop:
        evals.append(f_func(p))

    # Repeat C times
    for w in range(C):

        # Normalize the fitness
        pop_sum = sum(evals)
        evals_aux = []
        for ev in evals:
            evals_aux.append(ev / pop_sum)

        # Build the mating pool
        pool = []
        for i in range(len(pop)):
            prob = math.ceil(evals_aux[i] * 100)

            for j in range(prob):
                pool.append(pop[i])

        # Randomly select two parents
        p1 = random.choice(pool)
        p2 = random.choice(pool)
        # print("\nParent 1: ")
        # print_chromosome(p1)
        # print("\nParent 2: ")
        # print_chromosome(p2)

        # Select crossover method
        method = random.randint(0, 1)
        child = p1
        m = calc_m(pop[0])
        n = len(pop[0])

        # Apply vertical crossover
        if method == 0:
            rand_v = random.randint(1, m - 1)
            # print("\nSlice vertically by: " + str(rand_v))

            for i in range(len(p1)):
                child[i] = p1[i][:rand_v] + p2[i][rand_v:]

        # Apply horizontal crossover
        else:
            rand_h = random.randint(1, n - 1)
            # print("\nSlice horizontally by: " + str(rand_h))
            child = p1[:rand_h] + p2[rand_h:]

        # Select mutation method and apply
        if round(random.uniform(0, 1), 2) < mutation_rate:
            method = random.randint(0, 1)

            # Apply removal
            if method == 0:
                cell_i = random.randint(1, n - 1)
                cell_j = random.randint(1, len(child[cell_i]) - 1)

                while child[cell_i][cell_j] != "-":
                    cell_i = random.randint(0, n - 1)
                    cell_j = random.randint(0, len(child[cell_i]) - 1)

                child[cell_i] = child[cell_i][:cell_j] + child[cell_i][cell_j + 1:]
                # print("\nMutate: removed a gap from (" + str(cell_i + 1) + ", " + str(cell_j + 1) + ")")

            # Apply addition
            else:
                cell_i = random.randint(1, n - 1)
                cell_j = random.randint(1, len(child[cell_i]) - 1)

                child[cell_i] = child[cell_i][:cell_j] + "-" + child[cell_i][cell_j:]
                # print("\nMutate: added a gap in (" + str(cell_i + 1) + ", " + str(cell_j + 1) + ")")

        # Add the child to the new population
        child = add_gaps(child)
        child = remove_useless_gaps(child)
        new_pop.append(child)
        # print_chromosome(child)

    # Get the best chromosome
    best_val = 0
    best_chromosome = None
    for chromosome in new_pop:
        curr_val = f_func(chromosome)
        if curr_val >= best_val:
            best_val = curr_val
            best_chromosome = chromosome

    # Print stats
    print("Generation " + str(count + 1) + ": " + str(best_val))

    # Update population
    pop = copy.deepcopy(new_pop)
    pop = remove_spaces(pop)
    new_pop = []
    count = count + 1

# Best solution
print("\nBest solution:")
print_chromosome(best_chromosome)
