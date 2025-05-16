step 1 we extract our dataset from the MetaboLigths database and load it as an `MsExperiment` object. 

![img1](%5Cimages)

step 2 We next configure the parallel processing setup. Raw `.mzML` files contain millions of individual data points (mass-to-charge ratios, intensities, and retention times). processing them one by one on a single CPU core, can take hours. This code splits the workload so that **2 files are processed at the exact same time** using 2 separate cores of the processor.
img2

step 3
