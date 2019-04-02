# Class that handles utility methods
#
# More in detail, it handles mostly the addition, counting and
# removal of gaps and spaces from a matrix. It also reads and prepares
# the input matrix for the GA.

import math


class Utils:

    # Prints a chromosome's sequences
    @staticmethod
    def print_chromosome(chromosome):
        for seq in chromosome:
            print(seq)

    # Calc the maximum number of columns m
    @staticmethod
    def calc_m(lines):
        m_aux = 0
        lengths = []

        for line in lines:
            length = len(line)
            lengths.append(length)

            if length >= m_aux:
                m_aux = length

        return m_aux

    # Add spaces only to complete the columns of
    # a matrix and make calculations easier
    @staticmethod
    def add_spaces(lines):
        m = Utils.calc_m(lines)
        for i in range(len(lines)):
            diff = m - len(lines[i])
            for j in range(diff):
                lines[i] += " "

        return lines

    # Remove spaces to standardize
    @staticmethod
    def remove_spaces(pop):
        for chromosome in pop:
            for i in range(len(chromosome)):
                chromosome[i] = chromosome[i].replace(" ", "")

        return pop

    # Add the required gaps to produce an
    # initial alignment
    @staticmethod
    def add_gaps(lines):

        # Complete with gaps
        for i in range(len(lines)):
            diff = Utils.calc_m(lines) - len(lines[i])
            for j in range(diff):
                lines[i] += "-"

        return lines

    # Add the offset columns
    @staticmethod
    def add_offset(lines):
        offset = math.ceil(Utils.calc_m(lines) * 0.2)

        for i in range(offset):
            for j in range(len(lines)):
                lines[j] = "-" + lines[j]

        return lines

    # Remove useless gaps from an alignment
    @staticmethod
    def remove_useless_gaps(lines):
        to_rm = []
        only_gaps = True
        m = Utils.calc_m(lines)

        # Find columns with gaps only
        for i in range(m):
            for j in range(len(lines)):
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

    # Get the number of columns with only gaps
    @staticmethod
    def get_gap_cols(lines):
        n_cols = 0
        only_gaps = True
        n = len(lines)

        for i in range(Utils.calc_m(lines)):
            for j in range(n):
                if (lines[j][i] != " ") and (lines[j][i] != "-"):
                    only_gaps = False

            if only_gaps:
                n_cols += 1
            else:
                only_gaps = True

        return n_cols

    @staticmethod
    def prepare_input(input_path):

        # Read input file string
        with open(input_path) as f:
            input_str = f.read()

        # Get lines
        input_str = input_str.replace(" ", "")
        lines_list = input_str.split("\n")

        # Add the gaps
        lines_list = Utils.add_gaps(lines_list)
        return lines_list
