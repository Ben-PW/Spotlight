# Running this code yourself

This code should be trivial to run on your own machine. 
The only change you will need to make is to uncomment a section of code in 'Data_process.R' if you are running it for the first time.
The code in this repo is a demo version, it runs a scaled down design 2 x 2 x 4 x 3 = 48 levels, based on 2 empirical datasets with 100 simulated datasets, per condition,
per empirical dataset. 
The empirical datasets consist of a 16 node network and another network of approximately 200 nodes.
Once you have commented out the relevant section of code in 'Data_process.R', you can just run 'Spotlight_main.R'. The code should take no more than
2 - 3 minutes, and use no more than 1GB. Network level output is stored in the dataframe 'testResults'. Node level output is stored in 'nodeResults'.

# Planned improvements

I plan to update the loop so that, instead of 'nodeResults' being a dataframe of every centrality statistic for every node in every network (which, for this reduced version of the run, resulted in a dataset of over 1 million rows), the loop instead calculates the node level outcome measures and only stores them per network. This should substantially reduce the memory required to run the loop.

I have considered parallelisation for the main loop, however this might be difficult as it currently relies in incrementing the 'kg' and 'kn' counters. Multiple workers would likely conflict in incrementing these counters, meaning a re-write would be necessary. Currently I have saved on compute by computing spotlight assignment as early in the loop as possible, which means later combinations of missingness and b can be calculated on the same networks. I plan to add timing functions to the loop to give estimated completion times, as this will let me know if parallelisation is even necessary.

If memory becomes an issue, I plan to batch the calculations, or have the loop save results to disk instead of storing them in RAM. Saving to disk will likely be easier.

# Key design choices

Parameters
- spotlight_pcts (percentage of nodes selected to be spotlit)
- alphas (exponent selected for weighting spotlit assignment by node degree)
- ml (missingness level selected)
- b (weight assigned to non-spotlit ties, directly affects sampling probability for tie deletion)

For assigning spotlight by attribute, I have a slightly hacky plan. Attribute will be proxied by a dummy variable labelled 'Degree', as the spotlight assignment code expects a variable of that name to weight assignment probabilties by. 

The main design choice currently is how edge deletion probability is assigned. As is, the 'b' parameter controls the sampling weight assigned to each non-spotlit tie. Spotlit ties are assigned a weight of 1. The probability of a tie being sampled for deletion is, therefore, the product of the proportion of non spotlit ties, the missingness level selected, and the b parameter. Any feedback on this mechanism is welcome, as it is relatively easy to change at this stage of the code.
