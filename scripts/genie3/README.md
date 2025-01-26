Steps for running genie3 for a specific tissue (here, the tissue is liver)

1. Run './create_genie3_batches.sh liver'
2. Set tissue name to liver through 'tissue="liver"' line in genie3.sbatch
3. Run genie3.sbatch
4. Set tissue name to liver through 'tissue="liver"' line in genie3_post_processing.sh
5. Run genie3_post_processing.sh to generate matrix
6. Set tissue name to liver through 'tissue <- 'liver'' line in histogram_stuff.R
7. Run histogram_stuff.R (PLEASE RUN - IT CONTAINS A NECESSARY CHECK TO PREVENT AUTOREGULATION)
