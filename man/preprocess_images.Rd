% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess_images.R
\name{preprocess_images}
\alias{preprocess_images}
\title{Preprocess T1, FLAIR, EPI, and PHASE inputs}
\usage{
preprocess_images(
  t1_path,
  flair_path,
  epi_path,
  phase_path,
  output_dir,
  reorient = T,
  cores = 1,
  verbose = FALSE,
  return_images = T
)
}
\arguments{
\item{t1_path}{Path to the T1-weighted MRI image (in .nii.gz format).}

\item{flair_path}{Path to the FLAIR MRI image (in .nii.gz format).}

\item{epi_path}{Path to the EPI MRI image (in .nii.gz format).}

\item{phase_path}{Path to the phase MRI image (in .nii.gz format).}

\item{output_dir}{Directory where preprocessed images and results will be saved.}

\item{reorient}{Logical, indicating whether to reorient the images (default is TRUE).}

\item{cores}{Number of CPU cores to use for processing (default is 1).}

\item{verbose}{Logical, indicating whether to display verbose output (default is FALSE).}

\item{return_images}{Logical, indicating whether a named list of output images should be returned in addition to writing to disk (default is TRUE)}
}
\value{
Saves the following images to disk: t1_final.nii.gz, flair_final.nii.gz, epi_final.nii.gz, phase_final.nii.gz, prob.nii.gz, labeled_candidates.nii.gz.
If return_images = TRUE, also returns named list containing the images with names: t1, flair, epi, phase, prob_map, labeled_candidates.nii.gz. Named list can be used as input to \code{make_predictions}. If return_images = FALSE, returns NULL.
}
\description{
This function preprocesses inputs to feed to the neural network. Standard preprocessing includes bias correction, registration, skull-stripping, and intensity normalization to N(0, 1).
Additionally, a lesion probability map is calculated using MIMoSA (Valcarcel et al., 2018). The lesion probability map is thresholded very liberally at 0.05. Lesions are then labeled by lesion number, with confluent lesions being split using the Hessian (Dworkin et al., 2018).
}
\examples{
\dontrun{
preprocess_images("t1_image.nii.gz", "flair_image.nii.gz",
"epi_image.nii.gz", "phase_image.nii.gz",
"output_directory")
}
}
