# Distributed Computation Research

## Description
folders:
regularCodes - holds all paper-based codes:
  allAtOne - m workers each compute m functions together, take the time from the fastest worker
    uncoded - m workers each compute one function, take the time for all functions to be computed
    FR - there are m/d groups and in each group, d workers each compute d functions and take the time from the fastest worker of each group
    LT - m workers each compute me/m encoded symbols, take the time for all functions to be recoverable (when rank is equal to m)
  ratelessCodes - holds all improved codes:
    FR_RR - same as FR but upon finishing, each worker chooses the next function to compute using round robin scheme
    FR_rand - same as FR but upon finishing, each worker chooses the next function to compute by randomly picking
    uncoded_RR - special case of FR_RR where d = 1
    uncoded_rand - special case of FR_rand where d = 1
    BCC - there are m/d batches and each worker randomly picks a batch to compute, each batch is computed by the fastest worker, take the time for all batches to be computed
    LT_gaus - m workers each continuously compute encoded symbols until all functions are recoverable (when rank is equal to m)
    LT_oneVec - m workers each continuously compute encoded symbols until the vector of all 1s can be derived from the matrix
    SR_gaus - m workers each compute their assigned systematic symbol and then continuously compute encoded symbols until all functions are recoverable (when rank is equal to m)
    SR_oneVec - m workers each compute their assigned systematic symbol and then continuously compute encoded symbols until the vector of all 1s can be derived from the matrix
