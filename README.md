# BayesianVAR

# Origin of the real data set
We use data from ABIDE (Di Martino et al.,
2014;http://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html) preprocessed (Craddock et al., 2013;http://preprocessed-connectomes-project.org/abide/index.html) consisting of resting state fMRI data from 539 individuals
diagnosed with ASD and 573 healthy controls; we use randomly selected subsets
of 20 controls and 20 ASD patients from the data collected at New York University. The fMRI
data were collected using a 3 T Siemens Allegra scanner using a TR of 2 seconds. Each fMRI
dataset contains 180 time points. No motion scrubbing has been performed, but the first
four volumes were dropped in the processing to obtain 176 time points. The ABIDE Preprocessed
data have been processed with four different pipelines, and we use the data from the
CCS (connectome computation system) pipeline here. We use the data preprocessed without
global signal regression and without bandpass filtering, as bandpass filtering will substantially
change the autoregressive structure and we prefer to model it. Interested readers are
referred to ABIDE preprocessed for preprocessing details. As all the preprocessed data are
freely available, other researchers can reproduce our findings.

# References for the real data set
Craddock, C., Benhajali, Y., Chu, C., Chouinard, F., Evans, A., Jakab, A., Khundrakpam, B. S.,
Lewis, J. D., Li, Q., Milham, M., et al. (2013). The neuro bureau preprocessing initiative:
open sharing of preprocessed neuroimaging data and derivatives. Frontiers in Neuroinformatics,
7.

Di Martino, A., Yan, C.-G., Li, Q., Denio, E., Castellanos, F. X., Alaerts, K., Anderson, J. S.,
Assaf, M., Bookheimer, S. Y., Dapretto, M., et al. (2014). The autism brain imaging data
exchange: towards a large-scale evaluation of the intrinsic brain architecture in autism.
Molecular psychiatry, 19(6):659–667.
