#' Perform lesion splitting on custom lesion probability map
#'
#' This function allows for lesion splitting of a custom lesion probability
#' map. The probability map must be continuous from 0 to 1 and is originally
#' designed to work with MIMoSA probability maps. If deep learning
#' probability maps are used instead, range-constrained non-sigmoid
#' transformed outputs may work better, as sigmoid transforms may place
#' probability values too close to 0 and 1 for the Hessian-based approach
#' to appropriately split. The lesion splitting algorithm is from "An
#' Automated Statistical Technique for Counting Distinct Multiple Sclerosis
#' Lesions" (2018) by Jordan Dworkin et al.
#' 
#' @param custom_map_path Path to the custom probability map path (in .nii.gz format).
#' @param custom_map_threshold Value from 0 to 1 indicating threshold for voxel segmentation.
#' @param output_dir Directory where preprocessed images and results will be saved.
#' @param verbose Logical, indicating whether to display verbose output (default is FALSE).
#' @param return_images Logical, indicating whether a named list of output images should be returned in addition to writing to disk (default is TRUE)
#'
#' @return Saves the following images to disk: labeled_candidates.nii.gz, eroded_candidates.nii.gz
#' If return_images = TRUE, also returns named list containing the images with names: labeled_candidates.nii.gz, eroded_candidates.nii.gz.
#'
#' @import ANTsR
#' @importFrom stats predict
#' @import mimosa
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
split_custom_lesion_mask <- function(custom_map_path, custom_map_threshold,
                                     output_dir, verbose = FALSE,
                                     return_images = TRUE) {
  if (!inherits(custom_map_path, "character")) {
    stop("Must provide paths to .nii.gz files.")
  }
  if (!file.exists(custom_map_path)) {
    stop("Files are missing or paths or wrong")
  }
  if (custom_map_threshold <= 0 || custom_map_threshold >= 1) {
    stop("custom_map_threshold must be between 0 and 1")
  }
  if (!file.exists(output_dir)) {
    warning("Output directory does not exist. Making output directory")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  custom_map <- antsImageRead(custom_map_path)

  if (max(custom_map) > 1 || min(custom_map) < 0) {
    stop("custom_map values must be between 0 and 1")
  }

  # Threshold MIMoSA mask and identify/split confluent lesions
  prob_thresh <- antsImageClone(custom_map > custom_map_threshold)
  if (sum(prob_thresh) == 0) {
    prob_labeled <- antsImageClone(prob_thresh)
    prob_erode <- antsImageClone(prob_thresh)
  } else {
    prob_labeled <- oro2ants(label_lesion(custom_map,
                                          prob_thresh,
                                          mincluster = 30))
    prob_erode <- iMath(prob_labeled, "GE", 1)
  }
  antsImageWrite(prob_labeled,
                 file.path(output_dir, "labeled_candidates.nii.gz"))
  antsImageWrite(prob_erode,
                 file.path(output_dir, "eroded_candidates.nii.gz"))
  if (return_images) {
    list(
      labeled_candidates = prob_labeled,
      eroded_candidates = prob_erode
    )
  } else {
    NULL
  }
}
