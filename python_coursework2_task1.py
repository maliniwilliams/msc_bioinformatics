import random
import matplotlib.pyplot as plt

initial_pop_list = []    #empty list to put 50 A and 50 B into to create our population
popsize = 100    #the size of the population being sampled

for i in range(0, int(popsize/2)):   #add As to the list
    initial_pop_list.append("A")
for i in range(int(popsize/2), popsize):   #add Bs to the list
    initial_pop_list.append("B")
    
def sim_pop_gen_drift_task1(starting_population): #define a function to simulate the drift of A and B alleles
    freq_A = [] # [gen1, gen2, .. ]
    freq_B = []
    generation_nums = 1000
    
    prior_population = starting_population # So we can input the starting population into new_pop_list in for loop (initializing)
    
    count_A = 0 #count of A and B alleles in our initial population
    count_B = 0 
    for allele in starting_population: #count up the number of alleles in the starting population
            if allele == "A":
                count_A += 1 # count the As
            if allele == "B":
                count_B += 1 # count the Bs
    freq_A.append(count_A/100) #add the frequencies to the freq_A and freq_B lists as the starting frequency
    freq_B.append(count_B/100)
                      
    for generation in range(generation_nums): # this gives us 1000 generations of sampling
        new_pop_list = random.choices(prior_population, k=100) #sim 100 random samples with replacement from population as new_population
        count_A = 0 #count of A alleles in the new population
        count_B = 0 #count of B alleles in the new population
        for allele in new_pop_list: #count up the number of alleles in each population
            if allele == "A":
                count_A += 1
            if allele == "B":
                count_B += 1
        freq_A.append(count_A/100) #add the frequencies to the lists so they can be plotted
        freq_B.append(count_B/100)   
        prior_population = new_pop_list # Prior population (population of generation n) is now "prior_population" for generation n+1

        if (count_A == 0) or (count_B == 0):
            break
            
    return freq_A, freq_B

freq_A, freq_B = sim_pop_gen_drift_task1(initial_pop_list) #so we know that we are returning n A frequencies and n B frequencies

x_data = range(0, len(freq_A))  # x-axis is the length of freq_A, which should be the number of generations (same as len(freq_b))
y_data = freq_A  # y data for freq_A
y_data2 = freq_B  # y data for freq_B

plt.plot(x_data, y_data, 'r', label="Frequency of allele A")   # Plot freq_A
plt.plot(x_data, y_data2, 'b', label="Frequency of allele B")   # Plot freq_B

plt.xlabel("Generations")
plt.ylabel("Frequency")   # Add axis labels
plt.title("Frequency of A and B alleles over 1000 generations \n with random sampling")      # Add a title
plt.legend(loc='best')       # Add the key in the best position in the plot

plt.show()  # Display the plot