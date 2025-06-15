

# **End-to-End Untargeted Metabolomics Workflow**

I am learning by using this end-to-end workflow for LC-MS-based untargeted metabolomics experiments, conducted entirely within R using the Bioconductor ecosystem and base R functionality.

## ***Project Overview***

The samples were analyzed using ultra-high-performance liquid chromatography (UHPLC) coupled to a Q-TOF mass spectrometer (TripleTOF 5600+), utilizing hydrophilic interaction liquid chromatography (HILIC) for separation.

*Datasets used:*

 - LC-MS (MS1): Quantification of small polar metabolites in human
   plasma.
   
 - LC-MS/MS: Selected samples for feature identification and annotation.

*Experimental Design:*

 - Study Groups: Cardiovascular disease (CVD) patients vs. healthy
   controls (CTR).
 -  Samples: 3 CVD, 3 CTR, and 4 Quality Control (QC) samples.

 - Database Reference: MetaboLights MTBLS8735.

## **Step 1: Data Import**

We extract the dataset from MetaboLights and load it as an MsExperiment object. This object acts as a container for both metadata and spectral data.
![Data Import](images/img1.png)

## **Step 2: Data Organization & Preprocessing**

Raw .mzML files contain millions of data points. To improve efficiency, we configure parallel processing to split the workload across CPU cores.
![Parallel Processing](images/img2.png)
Metadata Sanitization
We extract the raw metadata, rename columns to human-readable formats (Sample, Type, Phenotype, Age), and integrate these into the lcms1 object.
![Metadata Table](images/img4.png)
Visual Encoding
To ensure consistency across all downstream plots, we apply a global color mapping using the RColorBrewer package, assigning specific colors to QC, CVD, and CTR groups.
![Color Palette](images/img5.png)
Spectral Inspection
The MS data is stored as a Spectra object. We inspect the metadata variables, including msLevel, rtime (retention time), precursorMz, and intensity.

![Spectra Object](images/img6.png)

![BPC Plot](images/img7.png)

Refining the chromatogram 

![Refined BPC](images/img8.png)

## **Step 3: Visualization of data (QA)**

BPC (Base Peak Chromatogram): Evaluates LC consistency.
TIC (Total Ion Chromatogram): Evaluates overall instrument sensitivity.

Effective visualization is paramount for assessing data quality. We used two internal standards:

First importing them both

![Importing IS](images/img9.png)

Then plotting

![EIC Comparison](images/img10.png)
we immediately notice clear concentration difference between QCs and study samples for the isotope labeled cystine ion. 
Also, the labeled methionine internal standard exhibits a signal amidst some noise and a noticeable retention time shift between samples.

## **Step 4 Data preprocessing**

3 major stages:

*1- Chromatographic Peak Detection*

To achieve this, we employ the findChromPeaks() function within xcms.
The preferred algorithm in this case, CentWave, utilizes continuous wavelet transformation (CWT)-based peak detection .
It is an effective way in handling non-Gaussian shaped chromatographic peaks or peaks with varying retention time widths, which are commonly encountered in HILIC separations.

Below we apply the CentWave algorithm  for cysteine and methionine ions and evaluate the results

![centwave](images/img11.png)

*peakwidth*: the minimal and maximal expected width of the peaks in the retention time dimension. dependent on the chromatographic settings 

*ppm*: The maximal allowed difference of mass peaks’ m/z values (in parts-per-million) in consecutive scans to consider them representing signal from the same ion.

*integrate*: the algorithm used to determine the peak boundaries. 

1 = assumes a symmetric, near gaussian peak shape while 
2 = determines the peak boundaries based on the observed intensities. 
We use integrate = 2 as it more accurately models non-symmetric chromatographic peaks with a longer right tail which are frequently observed in LC-MS data.

For this dataset the peak widths appear to be around 2 to 10 seconds. We choose a range that is not too wide or too narrow.

To determine the ppm, a deeper analysis of the dataset is needed. 
ppm depends on the instrument, 

Accurate determination of ppm : The following steps involve generating a highly restricted MS area with a single mass peak per spectrum, representing the cystine ion. The m/z of these peaks is then extracted, their absolute difference calculated and finally expressed in ppm

The following steps involve generating a highly restricted MS area with a single mass peak per spectrum, representing the cystine ion. The m/z of these peaks is then extracted, their absolute difference calculated and finally expressed in ppm.


![ppm](images/img11b.png)

We therefore, choose a value close to the maximum within this range for parameter `ppm`, i.e., 15 ppm.

Now that we have calibrated our parameters, we can perform the chromatographic peak detection on our EICs. When working with these restricted data windows, there are three critical factors to keep in mind:

 1. Retention Time Window Width**

Peak detection algorithms require sufficient "baseline" data points—the signal outside the peak—to accurately estimate background noise. If `findChromPeaks` fails to detect a feature in an EIC, the first troubleshooting step is to widen the retention time range of the EIC to capture more of the baseline.

2. Adjusting the Signal-to-Noise Threshold (`snthresh`)

When processing a whole LC-MS dataset, the algorithm has thousands of data points to calculate a robust noise level. Within a single, narrow EIC, that pool of data is much smaller. Consequently, you should often **reduce the `snthresh` value** for EICs to prevent the algorithm from prematurely rejecting real signals.

3. Handling Sparse Data Points

If peak appears "jagged" or consists of too few MS1 data points, the algorithm may struggle to model the peak shape. In such cases, setting `extendLengthMSW = TRUE` in the `CentWaveParam` object can help the algorithm extend the search area for wavelet scales, improving detection for sparser signals.

![eval](images/img13.png)

 - The columns `into` (integrated area) and `rt` (retention time) are  the primary quantitative outputs.  
   
 - L There is a consistent `sn`(signal-to-noise) ratio across all 10 samples indicating that `CentWaveParam` configuration is stable and capturing the signal  reliably.

*2- Retention Time (RT) Alignment*
3- Correspondence (Grouping=Match features across different samples)
