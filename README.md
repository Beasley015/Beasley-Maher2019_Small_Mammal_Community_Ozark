# Small mammal community composition varies among Ozark glades
Beasley, E.M. & Maher, S.P. In review. Small mammal community composition varies among Ozark glades. Journal of Mammalogy.

### For questions about data or code, please contact the corresponding author: Emily Beasley (ebeasley{at}uvm{dot}edu).

## Abstract
Island Biogeography Theory (IBT) explains and estimates large-scale ecological patterns among islands and isolated habitat patches. Specifically, IBT predicts that the number of species per habitat patch differs as a function of area and isolation as a result of local colonization and extinction. Accurate richness estimates are essential for testing predictions of IBT, but differences in detectability of species can lead to bias in empirical data sets. Hierarchical community models correct for imperfect detection by leveraging information from across the community to estimate species-specific occupancy and detection probabilities. Using the fragmented Ozark glades as our model system, we constructed a hierarchical community model to 1) estimate site-level and regional species richness of small mammals while correcting for detection error, and 2) determine environmental covariates driving occupancy. We sampled 16 glades in southwest Missouri in summer 2016–2017 and quantified mammal community structure within the glade network. The detected species pool included 8 species, and the model yielded a regional species estimate of 8.6 species, with a mean of 3.47 species per glade. Species richness increased with patch area but not isolation, and effects of patch shape varied between species in the community. The results reflect that the small mammal communities of the Ozark glades demonstrate patterns consistent with richness predictions of IBT but fails with respect to isolation. Discrepancies between the expectations and our data likely are the result of a permeable matrix reducing apparent isolation.

## Data
There are two zipped folders with data corresponding to the analyses. Raw data is contained in the folder MammIBTRawData. Standardized covariates and results of the Multi-Species Occupancy Model (MSOM) are contained in the folders ModelOutputsSampled (isolation measured between sampled glades only) and ModelOutputsAll (isolation measured between all glades in the sampling units). Text files with metadata are included in each folder.

## Code
MSOMFull.R - Full code for constructing and running a MSOM using R2OpenBugs. Code includes cleaning and preparing the raw data, writing
and running model script, and a few summary statistics. The model estimates species-specific occupancy and detection, evaluates effects 
of covariates on the above, and estimates total species richness N using data augmentation. Data is in the folder MammIBTRawData.

MSOMFigs.R - Code for generating figures and interpreting results of MSOM. Data is in the folders ModelOutputsSampled and ModelOutputsAll.
