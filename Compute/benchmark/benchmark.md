# Benchmarking SEISMIC-RNA


## Datasets

### Simulated Datasets

#### Length of reference
- 300
- 600
- 1200

#### Length of reads
- 150 (fragmented)
- 150x150 (amplicon)
- 150x150 (fragmented)

#### Clusters
- 1 cluster
- 2 clusters (50/50)
- 2 clusters (75/25)
- 2 clusters (90/10)
- 2 clusters (98/2)
- 3 clusters (33/33/33)
- 3 clusters (60/30/10)
- 3 clusters (50/45/5)
- 4 clusters (25/25/25/25)
- 4 clusters (40/30/20/10)
- 4 clusters (90/5/3/2)

#### Depth
- 10,000 reads
- 100,000 reads
- 1,000,000 reads

#### Mutation rate of unpaired bases
- 1%
- 3%
- 6%


## Benchmarks

### Features

- Record of settings
- Summary of results
- Replicate pooling
- Positional masking
- Read masking
- Clustering
- Structure modeling
- Graphs
- Data export


### Accuracy

- Population average reactivities
- Population average
- Clusters, only reads with no two mutations closer than 3 nt
- Correct number of clusters


### Speed

#### SEISMIC-RNA

- No masking
- Single masking
- Multiple maskings


### File sizes

#### RNA Framework

- No MM files
- With MM files