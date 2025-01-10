import random
import matplotlib.pyplot as plt

initial_pop_list = []    #empty list to put 25 AA, 50 Aa, and 25 aa into to create our population
popsize = 100    #the size of the population being sampled

for i in range(0, int(popsize/4)):   #add AAs to the list, the first 25
    initial_pop_list.append("AA")
for i in range(int(popsize/4), int(popsize-(popsize/4))):   #add Aas to the list, the middle 50
    initial_pop_list.append("Aa")
for i in range(int(popsize-(popsize/4)), int(popsize)): #add aas to the list, the last 25
    initial_pop_list.append("aa")

def sim_pop_gen_drift_task2(starting_population): #define a function to simulate the drift of A and a alleles
    freq_AA = [] # [gen1, gen2, .. ]
    freq_Aa = []
    freq_aa = []
    generation_nums = 500
    
    prior_population = starting_population # So we can input the starting population into new_pop_list in for loop (initializing)
    
    count_AA = 0 #count of AA, Aa, and aa genotypes in our initial population
    count_Aa = 0
    count_aa = 0
    for allele in starting_population: #count up the number of genotypes in the population
            if allele == "AA":
                count_AA += 1 # count the AAs
            if allele == "Aa":
                count_Aa += 1 # count the Aas
            if allele == "aa": 
                count_aa += 1 # count the aas
    freq_AA.append(count_AA/100) #add the frequencies to the freq_AA, freq_Aa, and freq_aa lists as the starting frequency
    freq_Aa.append(count_Aa/100)
    freq_aa.append(count_aa/100)
    
    for generation in range(generation_nums):
        new_population = []
        for individual in prior_population:
            mate_pair = random.sample(prior_population, 2) # random sample to sample w/out replacement as individual can't mate with self
            random_number = random.randint(1,5) #choose random number from 1-5, each number has a 20% chance of being selected
            for mate in mate_pair:
                #doing this for-loop because only 80% of the aa individuals will survive to sexual maturity and can be in a mating pair
                if (mate == "aa") and random_number == 1:
                    mate_pair = random.sample(prior_population, 2) #if the number is 1 (20% weight), repeat the sampling as this mate can't reproduce
                elif (mate == "aa") and random_number != 1: #if the number is 2-5 (80% weight), continue making offspring
                    continue 
                
            offspring = mate_pair[0][0] + mate_pair[1][0] #make offspring by choosing 1 allele from each parent
            new_population.append(offspring) #add the offspring to the new population
        
        count_AA = 0
        count_Aa = 0
        count_aa = 0
        for allele in new_population: #count up the number of alleles in each population
            if allele == "AA":
                count_AA += 1
            if (allele == "Aa") or (allele == "aA"): #there could be some aA and Aa combinations, but they are essentially the same for our purposes so put them together
                count_Aa += 1
            if allele == "aa":
                count_aa += 1
        freq_AA.append(count_AA/100) #add the frequencies to the lists so they can be plotted
        freq_Aa.append(count_Aa/100) 
        freq_aa.append(count_aa/100)
        
        if (count_Aa == 0) and (count_aa == 0): 
            break #stop sim when a disappears from the population
        
        prior_population = new_population #put the current generation back into the loop to run for n generations
           
    return freq_AA, freq_Aa, freq_aa # return list of frequencies of each genotype to be plotted

freq_AA, freq_Aa, freq_aa = sim_pop_gen_drift_task2(initial_pop_list) #so we know that we are returning n AA frequencies,
                                                                        #n Aa frequencies, and n aa frequencies

x_data = range(0, len(freq_AA))  # x-axis is the length of freq_AA, which should be the number of generations
y_data = freq_AA  # y data for freq_AA
y_data2 = freq_Aa  # y data for freq_Aa
y_data3 = freq_aa # y data for freq_aa

plt.plot(x_data, y_data, color = 'r', label="Frequency of allele AA")   # Plot freq_AA
plt.plot(x_data, y_data2, 'b', label="Frequency of allele Aa") # Plot freq_Aa
plt.plot(x_data, y_data3, 'g', label="Frequency of allele aa") # Plot freq_aa

plt.xlabel("Generations")
plt.ylabel("Frequency")   # Add axis labels
plt.title("Frequency of AA, Aa, and aa genotypes over " + str(len(freq_AA)) + " generations\n in a breeding population ")
plt.legend(loc="best")
plt.show()


