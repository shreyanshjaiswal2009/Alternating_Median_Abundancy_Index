import numpy as np
import heapq
from decimal import Decimal, getcontext
import time


#Introduction to Terminology: Limit. The limit is the value of N, till where one wishes to compute the alternating Median abundancy index.
#The higher the limit, the more better the value.
#We set a precision of 100 decimals for accurate calculations.
getcontext().prec = 100

# Precompute divisors for all numbers up to the limit. We precompute divisors and store it in the RAM, so that, computing indices becomes faster.
def precompute_divisors(limit):
    divisors = [[] for _ in range(limit + 1)]
    for i in range(1, limit + 1):
        for j in range(i, limit + 1, i):
            divisors[j].append(i)
    return divisors
#The divisors are arranged in an ascending order. Due to this, we essentially compute the 'reverse' written $\chi_+(n)$ (as mentioned in the manuscript), in the following function. We then take abs(\chi_+(n)) to get = $\chi(n)$

#Compute $\chi(n), which is the alternating divisor sum$.
def alternating_divisor_sum(n,divisors):
    total = 0
    for i, d in enumerate(divisors[n], start=1):
        if i % 2 == 0: #Eveny 2nd term in the divisor sum is a -ve.
            total -= d
        else:
            total += d
    return abs(total) #To get $\chi(n)$.

# Compute the alternating abundancy index.
def alternating_abundancy_index(n,divisors):
    divisor_sum = alternating_divisor_sum(n,divisors)
    return Decimal(divisor_sum) / Decimal(n) if n != 0 else Decimal(0)

# Main function to compute the median for the alternating abundancy indices.
def compute_median_index(limit):
    divisors = precompute_divisors(limit)
    min_heap = []  #The heaps are used for finding the median.
    max_heap = []  
    
    # Calculate specific powers of 10 to check. We print the progress of the alternating median abundancy index at 10^k, for k = 0,1,2..
    log_checkpoints = [10**i for i in range(int(np.log10(limit)) + 1)]
    
    median_candidate = Decimal(0) #Start with a median of 0
    for n in range(1, limit + 1):
        I_k = alternating_abundancy_index(n,divisors) #Calculate the alternating abundancy indices (for the given n)
        
        # Insert into heaps for dynamic median calculation
        if not max_heap or I_k <= -max_heap[0]:
            heapq.heappush(max_heap, -I_k)
        else:
            heapq.heappush(min_heap, I_k)
        
        # Balance the heaps
        if len(max_heap) > len(min_heap) + 1:
            heapq.heappush(min_heap, -heapq.heappop(max_heap))
        elif len(min_heap) > len(max_heap):
            heapq.heappush(max_heap, -heapq.heappop(min_heap))
        
        # Calculate median
        if len(max_heap) == len(min_heap):
            median_candidate = (-max_heap[0] + min_heap[0]) / 2
        else:
            median_candidate = -max_heap[0]
        
        # Print only at the logarithmic checkpoints
        if n in log_checkpoints:
            proportion_greater = Decimal(len(min_heap)) / Decimal(n)
            if proportion_greater == 0.5: #The output value of the proportion of indices greater than the median index should exactly be 0.5, as per the definition of the median index.
                print(f"Results verified.","N =",n, "Median candidate = ",median_candidate) #We print the median candidate at regular intervals.
            elif n == 1 and proportion_greater == 0: #One expection case for N = 1 is manually filtered out.
                print(f"Results verified.","N =",n, "Median candidate = ",median_candidate) 
            else:
                print("Unexpected error.")
            #If this proportion is not 0.5, indicates an issue in the computation. This feature acts as a simple embedded error detection mechanism.
    
    return median_candidate

#Final changes to be done here.
if __name__ == "__main__":
    limit = 10**6  # Adjust the limit as needed. Please be cautious, the program Will consume a majority of the system's RAM. Example limit set here = 10^5.
    #Limit should always be 10**N for a natural number N. We place this requirement, as the program only prints for N=10^k, for k=0,1,2,3... etc.
    start_time = time.time() #Start the time to obtain further data.
    median_result = compute_median_index(limit)
    end_time = time.time()
    print(f"\n\nFinal median:{median_result:.50f} after taking limit N =",limit) #Print the final median index with a precision of 50 decimal places.
    print(f"Execution time:{end_time - start_time:.2f} seconds\n") #Prints the total time used by the computation.
