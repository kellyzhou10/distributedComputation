# Distributed Computation Research

## Description

## Contents
### Regular Codes
#### allAtOne: m workers each compute m functions together, take the time from the fastest worker
#### uncoded: m workers each compute one function, take the time for all functions to be computed
#### FR: there are m/d groups and in each group, d workers each compute d functions and take the time from the fastest worker of each group
#### LT: m workers each compute me/m encoded symbols, take the time for all functions to be recoverable (when rank is equal to m)
### Rateless Codes
#### FR_RR: same as FR but upon finishing, each worker chooses the next function to compute using round robin scheme
#### FR_rand: same as FR but upon finishing, each worker chooses the next function to compute by randomly picking
#### BCC: there are m/d batches and each worker randomly picks a batch to compute, each batch is computed by the fastest worker, take the time for all batches to be computed
#### LT_largem: m workers each continuously compute encoded symbols until the vector of all 1s can be derived from the matrix; uses m = 1000
#### LT_smallm: same as LT_largem but uses m = 30
#### SR_largem: m workers each compute their assigned systematic symbol and then continuously compute encoded symbols until the vector of all 1s can be derived from the matrix; uses m = 1000
#### SR_smallm: same as SR_largem but uses m = 30
