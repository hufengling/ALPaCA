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
#' @param brain_mask_path Path to a brain mask (in .nii.gz format, default is NULL).
#' @param custom_lesion_map Path to custom lesion probability_map to use instead of MIMoSA output (in .nii.gz format, default is NULL). Lesion probability map must be between 0 and 1. If a custom lesion map is supplied, MIMoSA will be skipped and subsequent lesion mask processing will be performed on the custom map. A desired threshold must be provided to use with the custom lesion map. ALPaCA classification is not guaranteed to generalize to custom lesion maps.
#' @param custom_map_threshold Value between 0 to 1 at which to threshold custom_lesion_map. For best chances that ALPaCA will generalize well, a relatively lower threshold (high sensitivity, low specificity) should be chosen. (Default is NULL).
#' @param reorient Logical, indicating whether to reorient the images (default is TRUE).
#' @param cores Number of CPU cores to use for processing (default is 1).
#' @param verbose Logical, indicating whether to display verbose output (default is FALSE).
#' @param return_images Logical, indicating whether a named list of output images should be returned in addition to writing to disk (default is TRUE)
#'
#' @return Saves the following images to disk: t1_final.nii.gz, flair_final.nii.gz, epi_final.nii.gz, phase_final.nii.gz, prob.nii.gz, labeled_candidates.nii.gz, eroded_candidates.nii.gz.
#' If return_images = TRUE, also returns named list containing the images with names: t1, flair, epi, phase, prob_map, labeled_candidates.nii.gz, and eroded_candidates.nii.gz. Named list can be used as input to \code{make_predictions}. If return_images = FALSE, returns NULL.
#'
#' @import ANTsR
#' @importFrom stats predict
#' @importFrom mimosa mimosa_data
#' @importFrom fslr fslsmooth
#' @importFrom neurobase niftiarr read_rpi
#' @importFrom extrantsr ants2oro oro2ants fslbet_robust
#' @importFrom WhiteStripe whitestripe whitestripe_norm
#'
#' @export
#'
#' @examples \dontrun{
#' preprocess_images(
#'   "t1_image.nii.gz", "flair_image.nii.gz",
#'   "epi_image.nii.gz", "phase_image.nii.gz",
#'   "output_directory"
#' )
#' }
preprocess_images <- function(t1_path, flair_path, epi_path, phase_path,
                              output_dir, brain_mask_path = NULL,
                              custom_lesion_map = NULL, custom_map_threshold = NULL,
                              reorient = TRUE, cores = 1, verbose = FALSE,
                              return_images = TRUE) {
  if (!inherits(c(t1_path, flair_path, epi_path, phase_path), "character")) {
    stop("Must provide paths to .nii.gz files.")
  }
  if (!all(file.exists(c(t1_path, flair_path, epi_path, phase_path)))) {
    stop("Files are missing or paths or wrong")
  }
  if (!file.exists(output_dir)) {
    warning("Output directory does not exist. Making output directory")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  if (!is.null(custom_lesion_map)) {
    if (!inherits(custom_lesion_map, "character")) {
      stop("Must provide paths to .nii.gz files.")
    }
    if (!file.exists(custom_lesion_map)) {
      stop("Files are missing or paths or wrong")
    }
    if (is.null(custom_map_threshold)) {
      stop("If custom_lesion_map is provided, custom_map_threshold cannot be NULL. Please provide a desired threshold.")
    }
    if (custom_map_threshold <= 0 || custom_map_threshold >= 1) {
      stop("custom_map_threshold must be between 0 and 1")
    }
  }
  if (!file.exists(output_dir)) {
    warning("Output directory does not exist. Making output directory")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # Read files
  if (reorient) {
    t1 <- oro2ants(read_rpi(t1_path, verbose = verbose))
    flair <- oro2ants(read_rpi(flair_path, verbose = verbose))
    epi <- oro2ants(read_rpi(epi_path, verbose = verbose))
    phase <- oro2ants(read_rpi(phase_path, verbose = verbose))
  } else {
    t1 <- antsImageRead(t1_path)
    flair <- antsImageRead(flair_path)
    epi <- antsImageRead(epi_path)
    phase <- antsImageRead(phase_path)
  }

  # N4 bias correction
  t1 <- n4BiasFieldCorrection(t1, verbose = verbose)
  flair <- n4BiasFieldCorrection(flair, verbose = verbose)
  epi <- n4BiasFieldCorrection(epi, verbose = verbose)
  phase <- n4BiasFieldCorrection(phase, verbose = verbose)

  # Register T1 and FLAIR to EPI space. Change phase metadata to EPI (since it can be a tiny bit off)
  t1_reg <- antsRegistration(epi, t1, typeofTransform = "Rigid")
  t1_reg <- antsApplyTransforms(
    fixed = epi, moving = t1,
    transformlist = c(t1_reg$fwdtransforms),
    interpolator = "lanczosWindowedSinc"
  )
  flair_reg <- antsRegistration(epi, flair, typeofTransform = "Rigid")
  flair_reg <- antsApplyTransforms(
    fixed = epi, moving = flair,
    transformlist = c(flair_reg$fwdtransforms),
    interpolator = "lanczosWindowedSinc"
  )
  phase <- antsCopyImageInfo(epi, phase)

  # Brain extraction
  if (is.null(brain_mask_path)) {
    mask <- fslbet_robust(t1_reg) > 0
  } else {
    mask <- antsImageRead(brain_mask_path)
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

  # Create/split lesion maps
  if (is.null(custom_lesion_map)) {
    # WhiteStripe T1 and FLAIR for MIMoSA
    t1_reg_oro <- ants2oro(t1_reg)
    t1_ws <- whitestripe_norm(
      t1_reg_oro,
      whitestripe(t1_reg_oro,
                  "T1",
                  stripped = TRUE,
                  verbose = verbose
      )$whitestripe.ind
    )
    flair_reg_oro <- ants2oro(flair_reg)
    flair_ws <- whitestripe_norm(
      flair_reg_oro,
      whitestripe(flair_reg_oro,
                  "T2",
                  stripped = TRUE,
                  verbose = verbose
      )$whitestripe.ind
    )

    # Run MIMoSA
    mimosa_output <- mimosa_data(
      brain_mask = mask,
      FLAIR = flair_ws, T1 = t1_ws,
      gold_standard = NULL, normalize = "no",
      cores = cores, verbose = verbose
    )
    predictions_WS <- predict(mimosa_model,
                              mimosa_output$mimosa_dataframe,
                              type = "response"
    )
    predictions_nifti_WS <- niftiarr(mimosa_output$top_voxels, 0)
    predictions_nifti_WS[mimosa_output$top_voxels == 1] <- predictions_WS
    prob_map <- oro2ants(
      fslsmooth(predictions_nifti_WS,
                sigma = 1.25,
                mask = mimosa_output$tissue_mask,
                retimg = TRUE, smooth_mask = TRUE, verbose = verbose
      )
    )
    antsImageWrite(prob_map, file.path(output_dir, "prob.nii.gz"))

    # Threshold MIMoSA mask and identify/split confluent lesions
    prob_thresh <- antsImageClone(prob_map > 0.05)
    if (sum(prob_thresh) == 0) {
      prob_labeled <- antsImageClone(prob_thresh)
      prob_eroded <- antsImageClone(prob_thresh)
    } else {
      prob_labeled <- oro2ants(label_lesion(prob_map, prob_thresh, mincluster = 30))
      prob_eroded <- iMath(prob_labeled, "GE", 1)
    }
    antsImageWrite(prob_labeled,
                   file.path(output_dir, "labeled_candidates.nii.gz"))
    antsImageWrite(prob_eroded,
                   file.path(output_dir, "eroded_candidates.nii.gz"))
  } else {
    custom_split_output <- split_custom_lesion_mask(custom_lesion_map,
                                                    custom_map_threshold,
                                                    output_dir,
                                                    verbose = verbose)
    prob_labeled <- custom_split_output$labeled_candidates
    prob_eroded <- custom_split_output$eroded_candidates
  }

  if (return_images) {
    list(
      t1 = t1_final,
      flair = flair_final,
      epi = epi_final,
      phase = phase_final,
      prob_map = prob_map,
      labeled_candidates = prob_labeled,
      eroded_candidates = prob_eroded
    )
  } else {
    NULL
  }
}
