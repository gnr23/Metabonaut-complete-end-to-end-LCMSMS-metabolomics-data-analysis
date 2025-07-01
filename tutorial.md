

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

## **Step 1: Data Import**, organization and visualization

We extract the dataset from MetaboLights and load it as an MsExperiment object. This object acts as a container for both metadata and spectral data.

![Data Import](images/img1.png)

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


BPC (Base Peak Chromatogram): Evaluates LC consistency.
TIC (Total Ion Chromatogram): Evaluates overall instrument sensitivity.

Effective visualization is paramount for assessing data quality. We used two internal standards:

First importing them both

![Importing IS](images/img9.png)

Then plotting

![EIC Comparison](images/img10.png)
we immediately notice clear concentration difference between QCs and study samples for the isotope labeled cystine ion. 
Also, the labeled methionine internal standard exhibits a signal amidst some noise and a noticeable retention time shift between samples.

## **Step 2 Data preprocessing**

3 major stages:

## *1- Chromatographic Peak Detection*

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

Now we have calibrated our parameters and we can perform the chromatographic peak detection on our EICs. 

there are three critical factors :

 *1. Retention Time Window Width***

Peak detection algorithms require sufficient "baseline" data points—the signal outside the peak—to accurately estimate background noise. If `findChromPeaks` fails to detect a feature in an EIC, the first troubleshooting step is to widen the retention time range of the EIC to capture more of the baseline.

*2. Adjusting the Signal-to-Noise Threshold (`snthresh`)*

When processing a whole LC-MS dataset, the algorithm has thousands of data points to calculate a robust noise level. Within a single, narrow EIC, that pool of data is much smaller. Consequently, you should often **reduce the `snthresh` value** for EICs to prevent the algorithm from prematurely rejecting real signals.

*3. Handling Sparse Data Points*

If peak appears "jagged" or consists of too few MS1 data points, the algorithm may struggle to model the peak shape. In such cases, setting `extendLengthMSW = TRUE` in the `CentWaveParam` object can help the algorithm extend the search area for wavelet scales, improving detection for sparser signals.


![eval](images/img13.png)

 - The columns `into` (integrated area) and `rt` (retention time) are  the primary quantitative outputs.  
   
 - There is a consistent `sn`(signal-to-noise) ratio across all 10 samples indicating that `CentWaveParam` configuration is stable and capturing the signal reliably.
 
moving from **targeted EIC-based tuning** to **global processing** of the entire LC-MS dataset. Reverting `snthresh` to the default is a professional choice, as the algorithm now has the full, statistically robust context of the entire run to distinguish peaks from noise.

**Global Chromatographic Peak Detection**

 - By applying `findChromPeaks()` to the `lcms1` object (the full dataset), `xcms` will now identify features across the entire experiment. 

 - The use of `chunkSize` is optimization, ensuring that memory usage remains stable even as your dataset scales.

Chromatographic peak detection on EICs after processing.

![chunksize](images/img14.png)

*conclusion*: Peaks were detected  in all samples for both ions indicates that the peak detection process over the entire dataset was effective.

***Refining Peak Detection Artifacts***

In complex LC-MS data, "peak splitting" is a common occurrence. This usually happens when the algorithm identifies a single chromatographic elution event as two distinct peaks due to a momentary dip in signal intensity (often caused for example by electronic noise or poor ionization efficiency).

Refining peak detection is a vital step in metabolomics, 
`refineChromPeaks` with `MergeNeighboringPeaksParam`, aims at merging such split peaks.moving beyond simple detection into structural data validation.

peak artefacts detection 

![arrtefacts](images/img15.png)

![arrtefacts2](images/img16.png)

-   **`expandMz`:** we keep this very, very small. If too big, might accidentally grab a different isotope that belongs to a different molecule.
    
-   **`expandRt`:** 2.5 seconds before /after peak. That is about half the length of a normal peak,

![mergeneigh](images/img17.png)

Examples of CentWave peak detection artifacts after merging

![mergeneigh2](images/img18.png)

Therefore, we next apply these settings to the entire dataset. 
After peak merging, column `"merged"` in the result object’s `

![wholedataset](images/img19.png)

[chromPeakData((https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)` data frame can be used to evaluate which peaks in the result represent signal from merged, or from the originally identified chromatographic peaks.

![chropmpeakdata](images/img20.png)

Before alignment we verify that your Internal Standards are consistently detected
 
By using `chromatogram()` to extract peaks based on`intern_standard` table, we are effectively "tagging" EICs with the refined peak data. = easier inspect of specific ions later.

![chromatogram1](images/img21.png)

Evaluation of problematic samples 

Below we can see presentation of samples and number of identified chromatographic peaks.

![eval2](images/img22.png)

***Summary of chromatographic peak detection results per sample.***

![summary](images/img23.png)

Summary of chromatographic peak detection results per sample.

***



## *2- Peak Retention Time (RT) Alignment*

Even with ideal instrument conditions, experimental random variations will cause the peaks RT to shift across the chromatogram.

![plotallign](images/img24.png)

![plotallign2](images/img25.png)

BPC of QC samples.

Theoretically, we proceed in two steps: first we select only the QC samples of our dataset and do a first alignment on these, using the so-called _anchor peaks_. In this way we can assume a linear shift in time, 

Despite having external QCs in our data set, we still use the subset-based alignment assuming retention time shifts to be independent of the different sample matrix (human serum or plasma) and instead are mostly instrument-dependent. _[xcms](https://bioconductor.org/packages/3.22/xcms)_ package.

For the initial correspondence, we use the _PeakDensity_ approach ([Louail et al. 2025](https://rformassspectrometry.github.io/Metabonaut/articles/a-end-to-end-untargeted-metabolomics.html#ref-louail_xcms_2025)) that groups chromatographic peaks with similar _m/z_ and retention time into LC-MS features. The parameters for this algorithm, that can be configured using the `PeakDensityParam` object, are `sampleGroups`, `minFraction`, `binSize`, `ppm` and `bw`.

-   `binSize`  and  `ppm`  define the required similarity of  _m/z_  values. 
    
-   High density areas are identified using the base R  `[density()](https://rdrr.io/r/stats/density.html)`  function, for which  `bw`  is a parameter: higher values define wider retention time areas, lower values require chromatographic peaks to have more similar retention times. 

This code effectively visualizes the **High Density Areas** (the black line) vs. the actual **Chromatographic Peaks** detected in the samples.

![plotallign3](images/img26.png)

To optimize the `PeakDensityParam`, we performed a targeted visual inspection of _Cystine_. By setting `bw = 2` and `minFraction = 0.9`, we ensured that the peak grouping algorithm correctly identified the high-density retention time regions across all QC samples, creating a robust set of anchor peaks for subsequent alignment.

By basing thealignment on the **QC subset** using `PeakGroupsParam`, we will map how the peaks retention time drifted over time.
parameters used will be:

The parameters for this algorithm are:

-   `subsetAdjust`  and  `subset`: Allows for subset alignment. Here we base the retention time alignment on the QC samples, i.e., retention time shifts will be estimated based on these repeatedly measured samples.  `

-  To align the samples, we need "anchor peaks"—peaks that are unambiguously present and identifiable across the samples.`minFraction` sets the stringency for these anchors. If we set `minFraction = 0.9`, a feature must be detected in at least 90% of the QC subset to be considered a "trustworthy" anchor. 
    
-   `span`: Not all RT shifts are linear. Sometimes the drift happens more aggressively at the beginning of the chromatogram than at the end. To correct this, we use **LOESS (Locally Estimated Scatterplot Smoothing)** regression.

![plotallign4](images/img27.png)

Once the alignment has been performed, the user should evaluate the results using the `[plotAdjustedRtime()](https://rdrr.io/pkg/xcms/man/plotAdjustedRtime.html)` function.

![plotallign5](images/img28.png)

Retention time alignment results. There is only a small retention time shift since all samples were measured within the same measurement run.

We can also compare the BPC before and after alignment. To get the original data, i.e. the raw retention times, we can use the `[dropAdjustedRtime()](https://rdrr.io/pkg/xcms/man/XCMSnExp-class.html)` function:

![Beforeafter](images/img29.png)

BPC before and after alignment.

We notice that the largest shift can be observed in the retention time range from 120 to 130s.  Apart from that retention time range, only little changes can be observed.

Evaluating the impact of alignment on the **Internal Standards (IS)** .

![alignIS](images/img30.png)

EICs of cystine and methionine before and after alignment.

We conducted a targeted inspection of internal standards (e.g., Cystine and Methionine). While the alignment impact on Cystine was subtle , the improvement in Methionine was obvious
    
we performed a statistical audit using the **MetaboAnnotation** package.
    
    -   **Matching:** we used `matchValues()` with `MzRtParam` to objectively isolate peaks within a strict 50 ppm (m/z) and 10-second (RT) window.
        
    -   **Filtering:** we applied `SingleMatchParam` to ensure only the highest-intensity peaks (`target_maxo`) were selected, effectively filtering out background noise.
        
-   the  `intern_standard_peaks` object acts as an index of high-confidence reference peaks across your entire sample set.
With parameters `mzColname` and `rtColname` we specify the column names in the query (our IS) and target (chromatographic peaks) that contain the _m/z_ and retention time values on which to match entities. We perform this matching separately for each sample. For each internal standard in every sample, we use the `[filterMatches()](https://rdrr.io/pkg/MetaboAnnotation/man/Matched.html)` function and the `[SingleMatchParam()](https://rdrr.io/pkg/MetaboAnnotation/man/Matched.html)` parameter to select the chromatographic peak with the highest intensity.

![metabotta](images/img31.png)

 We can now extract the retention times for these chromatographic peaks before and after alignment.
 
![plotallign5](images/img32.png)

We can now evaluate the impact of the alignment on the retention times of internal standards across the full data set:

![plotallign5](images/img33.png)

Retention time variation of internal standards before and after alignment.

Retention time alignment conclusion: Very slight reduction of the variation in retention times of I.S. by the alignment performed.

## *3- Correspondence (Grouping=Match features across different samples)*

