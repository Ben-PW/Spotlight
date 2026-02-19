# Running this code yourself

This code should be trivial to run on your own machine. 
The only change you will need to make is to uncomment a section of code in 'Data_process.R' if you are running it for the first time.
The code in this repo is a demo version, it runs a scaled down design 2 x 2 x 4 x 3 = 48 levels, based on 2 empirical datasets with 100 simulated datasets, per condition,
per empirical dataset. 
The empirical datasets consist of a 16 node network and another network of approximately 200 nodes.
Once you have commented out the relevant section of code in 'Data_process.R', you can just run 'Spotlight_main.R'. The code should take no more than
2 - 3 minutes, and use no more than 1GB. Network level output is stored in the dataframe 'testResults'. Node level output is stored in 'nodeResults'.
