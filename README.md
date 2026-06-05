# Automated Segmentation of Multiple Sclerosis Lesions, Paramagnetic Rims, and Central Vein Sign on MRI Provides Reliable Diagnostic Biomarkers
## Install from the most updated Release
--------
**Maintainer**: Fengling Hu, fengling.hu@pennmedicine.upenn.edu

## Table of content
- [1. Installation](#id-section1)
- [2. Background](#id-section2)

<div id='id-section1'/>

## 1. Installation
The R package uploaded here does NOT directly contain the pretrained neural network parameters due to GitHub space constraints. In order to download and install the whole package, including neural network parameters, please follow below instructions.

1) Download and unzip the ALPaCA.zip file
2) Download the models (models_1_5.zip and models_6_10.zip) from the GitHub Release. Because of GitHub file size requirements, these models could not be version-controlled and shipped with the rest of the package.
3) Unzip both of the model zip files. Each zip file should include 10 models (5 autoencoder_\*.pt and 5 predictor_\*.pt files). Directly copy all 20 .pt files to the inst/extdata/ folder. The final format of this folder should be: inst/extdata/autoencoder_1.pt, inst/extdata/autoencoder_2.pt, etc.
4) Install the package locally (either using devtools::install_local("/path/to/ALPaCA/folder") within R or `R CMD INSTALL /path/to/ALPaCA/folder` from bash.)

Then, you can load this package via:

```
library(ALPaCA)
```

## 2. Background

This GitHub provides the package for "Automatic Segmentation and Classification of Multiple Sclerosis Lesions and Subtypes on Multi-Modal Magnetic Resonance Imaging Using a Convolutional Neural Network." Multiple sclerosis (MS) is a demyelinating, inflammatory disorder characterized by central nervous system lesions detectable via magnetic resonance imaging (MRI). Presence of classical MS lesions and two MS lesion subtypes – paramagnetic rim lesions (PRLs) and central vein sign lesions (CVSs) – are important for MS diagnosis and prognosis. However, manual segmentation of MS lesions, PRLs and CVS lesions is time-consuming and rater-dependent. We propose a fully automated method for segmenting MS lesions and subtypes, called Automated Lesion, PRL, and CVS Analysis network (ALPaCA).

Dependencies can be found in the NAMESPACE file – most of the packages required are standard neuroimaging packages. "torch" is the R version of PyTorch -- both packages use the libtorch backend. There are two main functions for the package. The first is “preprocess_images()” which will preprocess images from raw niftis. The function performs n4 bias correction, brain extraction, lesion candidate segmentation, and lesion splitting. The next function is “make_predictions()” which will take in the preprocessed images, create patches for each lesion candidate, and pass those through the ALPaCA network.

Additional documentation for the “preprocess_images()” and “make_predictions()” functions, as well as other non-user-facing helper functions can be accessed via the package documentation (i.e. "help()" function).
<div id='id-section3'/>
