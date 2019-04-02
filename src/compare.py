from utils import Utils
from ga import GA


# ClustalW output
clustal_1 = Utils.prepare_input("data/clustal_2.txt")
genetic_algorithm = GA(None, None, None, None, e=1.5)
print(genetic_algorithm.f_func(clustal_1))
