#' Preprocess T1, FLAIR, EPI, and PHASE inputs
#'
#' This function preprocesses inputs to feed to the neural network. Standard preprocessing includes bias correction, registration, skull-stripping, and intensity normalization to N(0, 1).
#' Additionally, a lesion probability map is calculated using MIMoSA (Valcarcel et al., 2018). The lesion probability map is thresholded very liberally at 0.05. Lesions are then labeled by lesion number, with confluent lesions being split using the Hessian (Dworkin et al., 2018).
#'
#' @param t1_path Path to the T1-weighted MRI image (in .nii.gz format).
#' @param flair_path Path to the FLAIR MRI image (in .nii.gz format).
#' @param epi_path Path to the EPI MRI image (in .nii.gz format).
#' @param phase_path Path to the phase MRI image (in .nii.gz format).
#' @param output_dir Directory where preprocessed images and results will be saved.
#' @param brainmask_path Path to a brain mask (in .nii.gz format, default is NULL).
#' @param reorient Logical, indicating whether to reorient the images (default is TRUE).
#' @param cores Number of CPU cores to use for processing (default is 1).
#' @param verbose Logical, indicating whether to display verbose output (default is FALSE).
#' @param return_images Logical, indicating whether a named list of output images should be returned in addition to writing to disk (default is TRUE)
#'
#' @return Saves the following images to disk: t1_final.nii.gz, flair_final.nii.gz, epi_final.nii.gz, phase_final.nii.gz, prob.nii.gz, labeled_candidates.nii.gz.
#' If return_images = TRUE, also returns named list containing the images with names: t1, flair, epi, phase, prob_map, labeled_candidates.nii.gz. Named list can be used as input to \code{make_predictions}. If return_images = FALSE, returns NULL.
#'
#' @importFrom stats predict
#' @import neurobase
#' @import mimosa
#' @import oro.nifti
#' @import fslr
#' @import ANTsR
#' @import ANTsRCore
#' @import extrantsr
#' @import WhiteStripe
#'
#' @export
#'
#' @examples \dontrun{
#' preprocess_images("t1_image.nii.gz", "flair_image.nii.gz",
#' "epi_image.nii.gz", "phase_image.nii.gz",
#' "output_directory")
#' }

preprocess_images <- function(t1_path, flair_path, epi_path, phase_path,
                              output_dir, brainmask_path = NULL,
                              reorient = T, cores = 1, verbose = FALSE,
                              return_images = T) {
  if (class(c(t1_path, flair_path, epi_path, phase_path)) != "character") {
    stop("Must provide paths to .nii.gz files.")
  }
  if (!all(file.exists(c(t1_path, flair_path, epi_path, phase_path)))) {
    stop("Files are missing or paths or wrong")
  }
  if (!file.exists(output_dir)) {
    warning("Output directory does not exist. Making output directory")
    dir.create(output_dir, showWarnings = F, recursive = T)
  }

  if (reorient) {
    t1 <- oro2ants(read_rpi(t1_path, verbose = verbose))
    flair <- oro2ants(read_rpi(flair_path, verbose = verbose))
    epi <- oro2ants(read_rpi(epi_path, verbose = verbose))
    phase <- oro2ants(read_rpi(phase_path, verbose = verbose))
  } else {
    t1 <- check_ants(t1_path)
    flair <- check_ants(flair_path)
    epi <- check_ants(epi_path)
    phase <- check_ants(phase_path)
  }

  # N4 bias correction
  t1 <- n4BiasFieldCorrection(t1, verbose = verbose)
  flair <- n4BiasFieldCorrection(flair, verbose = verbose)
  epi <- n4BiasFieldCorrection(epi, verbose = verbose)
  phase <- n4BiasFieldCorrection(phase, verbose = verbose)

  # Register T1 and FLAIR to EPI space. Change phase metadata to EPI (since it can be a tiny bit off)
  t1_reg <- antsRegistration(epi, t1, typeofTransform = "Rigid")
  t1_reg <- antsApplyTransforms(fixed = epi, moving = t1,
                                transformlist = c(t1_reg$fwdtransforms),
                                interpolator = "lanczosWindowedSinc")
  flair_reg <- antsRegistration(epi, flair, typeofTransform = "Rigid")
  flair_reg <- antsApplyTransforms(fixed = epi, moving = flair,
                                   transformlist = c(flair_reg$fwdtransforms),
                                   interpolator = "lanczosWindowedSinc")
  phase <- antsCopyImageInfo(epi, phase)

  # Brain extraction
  ## 6/26/25 - EAH
  ## fslbet isn't working well. trying fslbet_robust
  #mask <- fslbet(flair_reg) != 0
  if (is.null(brainmask_path)) {
     mask <- fslbet_robust(t1_reg) > 0
  } else {
     mask <- check_ants(brainmask_path)
  }
  t1_reg <- t1_reg * mask
  flair_reg <- flair_reg * mask
  epi <- epi * mask
  phase <- phase * mask

  # Calculate normalized images for DL network and write to storage
  t1_dist <- c(mean(t1_reg[mask]), sd(t1_reg[mask]))
  t1_final <- ((t1_reg - t1_dist[1]) / t1_dist[2]) * mask
  antsImageWrite(t1_final, file.path(output_dir, "t1_final.nii.gz"))

  flair_dist <- c(mean(flair_reg[mask]), sd(flair_reg[mask]))
  flair_final <- ((flair_reg - flair_dist[1]) / flair_dist[2]) * mask
  antsImageWrite(flair_final, file.path(output_dir, "flair_final.nii.gz"))

  epi_dist <- c(mean(epi[mask]), sd(epi[mask]))
  epi_final <- ((epi - epi_dist[1]) / epi_dist[2]) * mask
  antsImageWrite(epi_final, file.path(output_dir, "epi_final.nii.gz"))

  phase_dist <- c(mean(phase[mask]), sd(phase[mask]))
  phase_final <- ((phase - phase_dist[1]) / phase_dist[2]) * mask
  antsImageWrite(phase_final, file.path(output_dir, "phase_final.nii.gz"))

  # WhiteStripe T1 and FLAIR for MIMoSA
  t1_reg_oro <- ants2oro(t1_reg)
  t1_ws <- whitestripe_norm(t1_reg_oro,
                            whitestripe(t1_reg_oro,
                                        "T1", stripped = TRUE,
                                        verbose = verbose)$whitestripe.ind)
  flair_reg_oro <- ants2oro(flair_reg)
  flair_ws <- whitestripe_norm(flair_reg_oro,
                               whitestripe(flair_reg_oro,
                                           "T2", stripped = TRUE,
                                           verbose = verbose)$whitestripe.ind)

  # Run MIMoSA
  mimosa <- mimosa_data(brain_mask = mask,
                        FLAIR = flair_ws, T1 = t1_ws,
                        gold_standard = NULL, normalize = "no",
                        cores = cores, verbose = verbose)
  predictions_WS <- predict(mimosa_model,
                            mimosa$mimosa_dataframe,
                            type = "response")
  predictions_nifti_WS <- niftiarr(mimosa$top_voxels, 0)
  predictions_nifti_WS[mimosa$top_voxels == 1] <- predictions_WS
  prob <- oro2ants(
    fslsmooth(predictions_nifti_WS, sigma = 1.25,
              mask = mimosa$tissue_mask,
              retimg = TRUE, smooth_mask = TRUE, verbose = verbose)
  )
  antsImageWrite(prob, file.path(output_dir, "prob.nii.gz"))

  # Threshold MIMoSA mask and identify/split confluent lesions
  prob_05 <- antsImageClone(prob > 0.05)
  if (sum(prob_05) == 0) {
    prob_05_labeled <- antsImageClone(prob_05)
    prob_05_erode <- antsImageClone(prob_05)
  } else {
    prob_05_labeled <- oro2ants(label_lesion(prob, prob_05, mincluster = 30))
    prob_05_erode <- iMath(prob_05_labeled, "GE", 1)
  }
  antsImageWrite(prob_05_labeled, file.path(output_dir, "labeled_candidates.nii.gz"))
  antsImageWrite(prob_05_erode, file.path(output_dir, "eroded_candidates.nii.gz"))

  return(list(t1 = t1_final,
              flair = flair_final,
              epi = epi_final,
              phase = phase_final,
              prob_map = prob,
              labeled_candidates = prob_05_labeled,
              eroded_candidates = prob_05_erode)
  )
}
