# Class that implements the proposed genetic algorithm
# The initial population, selection, substitution, adaptation function,
# crossover and mutation operators are implemented in this class.

import nwalign3 as nw
import numpy as np
import random
import copy
import math
import utils


class GA:

    def __init__(self, chromosomes, generations, min_generations, mutation_rate, e):
        self.chromosomes = chromosomes
        self.generations = generations
        self.min_generations = min_generations
        self.mutation_rate = mutation_rate
        self.e = e

    # The g function defined in the paper
    @staticmethod
    def g_func(x, y):
        if (x != y) and (x != "-") and (y != "-"):
            return 0
        elif (x != y) and (x == "-" or y == "-"):
            return 1
        elif x == y:
            return 2

        return 0

    # Apply the f function defined in the paper
    def f_func(self, lines):
        value = 0
        lines = utils.Utils.add_spaces(lines)
        x = utils.Utils.get_gap_cols(lines)
        m = utils.Utils.calc_m(lines)

        for i in range(m):
            for j in range(len(lines) - 1):
                value += self.g_func(lines[j][i], lines[j + 1][i])

        value = value - 2 * x * (len(lines) - 1) - self.e * x
        return value

    # Returns if there are no relevant changes
    # after 100 generations
    def no_change(self, best):
        if len(best) < self.min_generations:
            return False
        else:
            percent = int(0.2 * len(best))
            last = best[-percent:]

            if np.var(last) < 1.05:
                return True

        return False

    # Create the initial population
    def init_pop(self, lines_list):
        pop = []
        for c in range(self.chromosomes):

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
            lines_list_aux = utils.Utils.add_gaps(lines_list_aux)
            pop.append(lines_list_aux)
            print("\nChromosome " + str(c + 1) + ":")
            utils.Utils.print_chromosome(lines_list_aux)

        # Initial population
        return pop

    # Evaluate all chromosomes
    def eval_all(self, pop):
        evaluations = []
        for p in pop:
            evaluations.append(self.f_func(p))

        return evaluations

    # Select two parents for a crossover operation
    @staticmethod
    def select_parents(pop, evaluations):

        # Normalize the fitness
        pop_sum = sum(evaluations)
        evaluations_aux = []
        for ev in evaluations:
            evaluations_aux.append(ev / pop_sum)

        # Build the mating pool
        pool = []
        for i in range(len(pop)):
            prob = math.ceil(evaluations_aux[i] * 100)

            for j in range(prob):
                pool.append(pop[i])

        # Randomly select two parents
        p1 = random.choice(pool)
        p2 = random.choice(pool)
        return [p1, p2]

    # Apply a crossover operation on p1 and p2
    @staticmethod
    def apply_crossover(pop, p1, p2):

        # Select crossover method
        method = random.randint(0, 1)
        child = p1
        m = utils.Utils.calc_m(pop[0])
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

        # Return the child
        return child

    # Apply a mutation on child
    def apply_mutation(self, pop, child):

        n = len(pop[0])

        # Select mutation method and apply
        if round(random.uniform(0, 1), 2) < self.mutation_rate:
            rand = round(random.uniform(0, 1), 2)

            # Apply removal
            if rand < 0.3:
                cell_i = random.randint(1, n - 1)
                cell_j = random.randint(1, len(child[cell_i]) - 1)

                while child[cell_i][cell_j] != "-":
                    cell_i = random.randint(0, n - 1)
                    cell_j = random.randint(0, len(child[cell_i]) - 1)

                child[cell_i] = child[cell_i][:cell_j] + child[cell_i][cell_j + 1:]
                # print("\nMutate: removed a gap from (" + str(cell_i + 1) + ", " + str(cell_j + 1) + ")")

            # Apply k addition
            else:
                cell_i = random.randint(1, n - 1)
                cell_j = random.randint(1, len(child[cell_i]) - 1)
                k = random.randint(2, math.ceil(0.1 * utils.Utils.calc_m(child)))
                to_add = ""

                for i in range(k):
                    to_add += "-"

                child[cell_i] = child[cell_i][:cell_j] + to_add + child[cell_i][cell_j:]
                # print("\nMutate: added a gap in (" + str(cell_i + 1) + ", " + str(cell_j + 1) + ")")

        # Return the child
        return child

    def run_ga(self, input_path):

        # Read input file string
        lines_list = utils.Utils.prepare_input(input_path)

        # Prints the original chromosome
        print("Input matrix:")
        utils.Utils.print_chromosome(lines_list)

        # Create the initial population
        pop = self.init_pop(lines_list)

        # Repeat for all generations or until a good solution
        # appears
        best_chromosome = None
        best_chromosomes = []
        count = 0
        new_pop = []
        print()

        while count < self.generations:

            # Evaluate all chromosomes
            evaluations = self.eval_all(pop)

            # Repeat for all chromosomes
            for w in range(self.chromosomes):

                # Select two parents
                p1, p2 = self.select_parents(pop, evaluations)

                # Apply a crossover operation on p1 and p2
                child = self.apply_crossover(pop, p1, p2)

                # Apply a mutation on child
                child = self.apply_mutation(pop, child)

                # Add the child to the new population
                child = utils.Utils.add_gaps(child)
                child = utils.Utils.remove_useless_gaps(child)
                new_pop.append(child)
                # utils.Utils.print_chromosome(child)

            # Get the best chromosome
            best_val = 0
            best_chromosome = None
            for chromosome in new_pop:
                curr_val = self.f_func(chromosome)
                if curr_val >= best_val:
                    best_val = curr_val
                    best_chromosome = chromosome

            # Print stats
            print("Generation " + str(count + 1) + ": " + str(best_val))

            # Add best to the list
            best_chromosomes.append(best_val)

            # Break the execution if there are no
            # relevant changes
            if self.no_change(best_chromosomes):
                print("\nAbort: no variation!", end="")
                break

            # Update population
            pop = copy.deepcopy(new_pop)
            pop = utils.Utils.remove_spaces(pop)
            new_pop = []
            count = count + 1

        # Best solution
        print("\nBest solution:")
        utils.Utils.print_chromosome(best_chromosome)
