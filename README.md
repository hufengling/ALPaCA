# Automated Segmentation of Multiple Sclerosis Lesions, Paramagnetic Rims, and Central Vein Sign on MRI Provides Reliable Diagnostic Biomarkers

--------
**Maintainer**: Fengling Hu, fengling.hu@pennmedicine.upenn.edu

## Table of content
- [1. Installation](#id-section1)
- [2. Background](#id-section2)

<div id='id-section1'/>

## 1. Installation
The R package uploaded here does not contain the pretrained neural network parameters due to GitHub space constraints. The full package can be installed from the following DOI link to **Zenodo**:

10.5281/zenodo.17215591

Once downloaded, please unzip the "ALPaCA.zip" package. Then, the package can be installed in R using the following command.

```
# install.packages("devtools")
devtools::install_local("/path/to/file/ALPaCA")
```

Then, you can load this package via:

```
library(ALPaCA)
```

## 2. Background

This GitHub provides the package for "Automatic Segmentation and Classification of Multiple Sclerosis Lesions and Subtypes on Multi-Modal Magnetic Resonance Imaging Using a Convolutional Neural Network." Multiple sclerosis (MS) is a demyelinating, inflammatory disorder characterized by central nervous system lesions detectable via magnetic resonance imaging (MRI). Presence of classical MS lesions and two MS lesion subtypes – paramagnetic rim lesions (PRLs) and central vein sign lesions (CVSs) – are important for MS diagnosis and prognosis. However, manual segmentation of MS lesions, PRLs and CVS lesions is time-consuming and rater-dependent. We propose a fully automated method for segmenting MS lesions and subtypes, called Automated Lesion, PRL, and CVS Analysis network (ALPaCA).

Dependencies can be found in the NAMESPACE file – most of the packages required are standard neuroimaging packages. "torch" is the R version of PyTorch -- both packages use the libtorch backend. There are two main functions for the package. The first is “preprocess_images()” which will preprocess images from raw niftis. The function performs n4 bias correction, brain extraction, lesion candidate segmentation, and lesion splitting. The next function is “make_predictions()” which will take in the preprocessed images, create patches for each lesion candidate, and pass those through the ALPaCA network.

Additional documentation for the “preprocess_images()” and “make_predictions()” functions, as well as other non-user-facing helper functions can be accessed via the package documentation (i.e. "help()" function).
<div id='id-section3'/>
