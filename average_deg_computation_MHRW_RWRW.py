'''
@EMDC 2011 - 2013 

Homework for ID2220	Advanced areas of distributed systems

Estimation of the average node degree with MHRW and RWRW

@authors: Vaidas Brundza & Aras Tarhan
'''
import random
import numpy
import sys
import networkx as netX

# Min/Max incomes predefined
min_income = 1000
max_income = 8000
step_income = 100
# ---------------
 
# a function for loading the file into memory as a graph
def load_graph(filename):
    G = netX.read_edgelist(filename,  delimiter='\t',  nodetype=int)
    return G

# a function to generate incomes randomly for a given graph
# the range is {min_income, max_income} with a step of 100
def generateNodesIncome(G):
    Income = {}
    for x in G.__iter__():
        Income[x] = random.randrange(min_income, max_income, step_income)
    return Income

# a function to perform Re-Weighted Random Walk
def reWeightedRW(G, income, sampleSize):
    node = random.choice(G.nodes())

    sampling = list()
    node_degrees = list()
    node_incomes = list()
    prob = list()

    for i in range(sampleSize):
        sampling.append(node)
        node_degrees.append(len(G[node]))
        node_incomes.append(income[node])

        # Select a random neighbor of node
        node = random.choice(G.neighbors(node))

    # The normal random walk. biased, without correction.
    biased_average_degrees = numpy.average(node_degrees)
    biased_average_incomes = numpy.average(node_incomes)

    # Correcting the random walk sampling with inverted-node-degree probability
    summation_of_denominator = 0.0
    for x in node_degrees:
        summation_of_denominator += (1.0 / x)

	numerator=0.0
	
	for i in range(len(node_degrees)):
		numerator += (node_incomes[i] * 1.0) / node_degrees[i]	
	reweighted_average_degrees = len(node_degrees)/summation_of_denominator
    reweighted_average_incomes = numerator/summation_of_denominator
	
    return [biased_average_degrees, reweighted_average_degrees, biased_average_incomes, reweighted_average_incomes]

# a function to perform Metropolis-Hasting Random Walk
def metropolisHastingsRW(G, incomes, sampleSize):
    # Sets the starting node for random walk
    node = random.choice(G.nodes())

    sampling = list()
    node_degrees = list()
    node_incomes = list()

    # Actual random walk 
    for i in range(sampleSize):
        sampling.append(node)
        node_degrees.append(len(G[node]))
        node_incomes.append(incomes[node])

        # Select a random neighbor of node for RW 
        neighbor = random.choice(G.neighbors(node))

        # Perform Metropolis-Hastings algorithm with additional coefficient: if probability is less than 20%, we set it to 20%
        if (len(G[node]) > len(G[neighbor])):
            node = neighbor
        else:
            rand = random.random()
            prob = float(len(G[node])) / len(G[neighbor])
            if (rand < prob):
                node = neighbor
            else:
                node = node
                
    # Print "average with duplicates"
    avg_degrees_w_duplicate = numpy.average(node_degrees)
    avg_incomes_w_duplicate = numpy.average(node_incomes)
    
    # Remove duplicates from sampled list
    degrees_wo_duplicate = dict(zip(sampling, node_degrees))
    incomes_wo_duplicate = dict(zip(sampling, node_incomes))

    # Print "average without duplicate"
    avg_degrees_wo_duplicate = numpy.average(degrees_wo_duplicate.values())
    avg_incomes_wo_duplicate = numpy.average(incomes_wo_duplicate.values())

    return [avg_degrees_w_duplicate, avg_degrees_wo_duplicate, avg_incomes_w_duplicate, avg_incomes_wo_duplicate]


# a function for gathering and writing RWRW and MHRW results to file
def resultPrint(G, income, number_of_tests, sample_sizes):
    source = open('results.txt', 'w')
    real_avg_degree = 2.0 * G.number_of_edges()/G.number_of_nodes()
    real_avg_income = numpy.average(income.values())
    
    # Data structures for results handling
    biased_degree = list()
    re_degree = list()
    biased_incomes = list()
    re_incomes = list()
    
    degree = list()
    degree_wo_duplicates = list()
    incomes = list()
    incomes_wo_duplicates = list()
    
    # Writes results of Re-Weighted Random Walk to file
    source.write("------------Results of Re-Weighted Random Walk------------\n\n")
    source.write('#sample mean_biased_deg mean_rw_deg     std_biased_deg  std_rw_deg      real_avg_deg    mean_biased_inc mean_rw_inc     std_biased_inc  std_rw_inc      real_avg_inc\n')
   
    for sample in sample_sizes:
        for i in range(number_of_tests):
            re_result = reWeightedRW(G, income, sample)
            biased_degree.append(re_result[0])
            re_degree.append(re_result[1])
            biased_incomes.append(re_result[2])
            re_incomes.append(re_result[3])
            
        # Mean of the estimated averages RWRW
        mean_biased_degree = numpy.average(biased_degree)
        mean_re_degree = numpy.average(re_degree)
        mean_biased_incomes = numpy.average(biased_incomes)
        mean_re_incomes = numpy.average(re_incomes)
        
        # Standard deviation of the estimated averages RWRW
        std_biased_degrees = numpy.std(biased_degree)
        std_re_degrees = numpy.std(re_degree)
        std_biased_incomes = numpy.std(biased_incomes)
        std_re_incomes = numpy.std(re_incomes)

        source.write('{0:0.0f}\t{1:2.7f}\t{2:2.7f}\t{3:2.7f}\t{4:2.7f}\t{5:2.7f}\t{6:2.7f}\t{7:2.7f}\t{8:2.7f}\t{9:2.7f}\t{10:2.7f}\t\n'.format(sample, mean_biased_degree, mean_re_degree, std_biased_degrees, std_re_degrees, real_avg_degree, mean_biased_incomes, mean_re_incomes, std_biased_incomes, std_re_incomes, real_avg_income))
         
    # Writes results of Metropolis-Hasting Random Walk to file
    source.write("\n------------Results of Metropolis-Hasting Random Walk------------\n\n")
    source.write('#sample mean_deg_wd     mean_deg_wod    std_deg_wd      std_deg_wod     real_avg_deg    mean_inc_wd     mean_inc_wod    std_inc_wd      std_inc_wod     real_avg_inc\n')
   
    for sample in sample_sizes:
        for i in range(number_of_tests):
            mh_result = metropolisHastingsRW(G, income, sample)
            degree.append(mh_result[0])
            degree_wo_duplicates.append(mh_result[1])
            incomes.append(mh_result[2])
            incomes_wo_duplicates.append(mh_result[3])

        # Mean of the estimated averages MHRW
        mean_degree_wd = numpy.average(degree)
        mean_degree_wod = numpy.average(degree_wo_duplicates)
        mean_incomes_wd = numpy.average(incomes)
        mean_incomes_wod = numpy.average(incomes_wo_duplicates)

        # Standard deviation of the estimated averages MHRW
        std_degree_wd = numpy.std(degree)
        std_degree_wod = numpy.std(degree_wo_duplicates)
        std_incomes_wd = numpy.std(incomes)
        std_incomes_wod = numpy.std(incomes_wo_duplicates)
        
        source.write('{0:0.0f}\t{1:2.7f}\t{2:2.7f}\t{3:2.7f}\t{4:2.7f}\t{5:2.7f}\t{6:2.7f}\t{7:2.7f}\t{8:2.7f}\t{9:2.7f}\t{10:2.7f}\t\n'.format(sample, mean_degree_wd, mean_degree_wod, std_degree_wd, std_degree_wod, real_avg_degree, mean_incomes_wd, mean_incomes_wod, std_incomes_wd, std_incomes_wod, real_avg_income))

def main(argv):
	# Input arguments: #noOfTests #noOfSamples & sampleSizes;
	number_of_tests = int(sys.argv[1])
	
	sample_sizes=list()

	for i in range(int(sys.argv[2])):
		sample_sizes.append(int(sys.argv[i+3]))
			
	G = load_graph("p2p-Gnutella31.txt")
	print("---> Generating income")
	income = generateNodesIncome(G)
	print("---> Printing results")
	resultPrint(G, income, number_of_tests, sample_sizes)
	print("---> Results have been published !")

if __name__ == '__main__':
    main(sys.argv)