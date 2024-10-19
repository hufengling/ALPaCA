#' #' Make Predictions using ALPaCA
#'
#' This function runs the pre-trained Automated Lesion, PRL (Paramagnetic Rim Lesion), and CVS (Central Vein Sign) Analysis network. This model takes in multi-modal data in the form of pre-processed T1, FLAIR, EPI, and phase images to generate predictions of whether lesion candidates (identified via MIMoSA) are true lesions, PRLs, or CVS.
#'
#' @param ants_list An optional named list containing input images or file paths to images. Recommended to be output from the preprocess_images() function. Names must be: t1, flair, epi, phase, and labeled_candidates. Either ants_list must be provided, or t1, flair, epi, phase, and labeled_candidates must be provided.
#' @param t1 antsImage or file path to .nii.gz representing T1-weighted MRI image.
#' @param flair antsImage or file path to .nii.gz representing the FLAIR MRI image.
#' @param epi antsImage or file path to .nii.gz representing the EPI MRI image.
#' @param phase antsImage or file path to .nii.gz representing the phase MRI image.
#' @param labeled_candidates antsImage or file path to .nii.gz representing labeled candidates for lesion regions.
#' @param eroded_candidates antsImage or file path to .nii.gz representing eroded candidates for lesion regions.
#' @param lesion_priority A character vector specifying priority for lesion prediction thresholds -- Youden's J, Specificity, Sensitivity. Thresholds are based on training set ROC curves from CV models. Default priority is Youden's J, with sensitivity \eqn{\approx} 0.83 and specificity \eqn{\approx} 0.86. 'Specificity' prioritizes specificity 3 times more than sensitivity, with sensitivity \eqn{\approx} 0.69 and specificity \eqn{\approx} 0.94. 'Sensitivity' prioritizes sensitivity 3 times more than specificity, with sensitivity \eqn{\approx} 0.92 and specificity \eqn{\approx} 0.70.
#' @param prl_priority A character vector specifying priority for PRL prediction thresholds. Same options and default as lesion_priority. For Youden's J, sensitivity \eqn{\approx} 0.76 and specificity \eqn{\approx} 0.83. For Specificity, sensitivity \eqn{\approx} 0.63 and specificity \eqn{\approx} 0.90. For Sensitivity, sensitivity \eqn{\approx} 0.86 and specificity \eqn{\approx} 0.64.
#' @param cvs_priority A character vector specifying priority for CVS prediction thresholds.Same options and default as lesion_priority. For Youden's J, sensitivity \eqn{\approx} 0.81 and specificity \eqn{\approx} 0.65. For Specificity, sensitivity \eqn{\approx} 0.27 and specificity \eqn{\approx} 0.91. For Sensitivity, sensitivity \eqn{\approx} 0.91 and specificity \eqn{\approx} 0.47.
#' @param return_raw_probabilities A logical flag indicating whether to return raw probability antsImages for each region. A raw probability lesion-wise dataframe is always returned. (Default is FALSE)
#' @param clear_discordant_predictions A logical flag indicating whether to clear discordant predictions (ie candidates where the model predicts "CVS"/"PRL, but not "lesion".) In training, lesion prediction was almost always more reliable under Youden's J threshold choice. (Default is TRUE)
#' @param n_patches An integer specifying the number of patches to sample for predictions. Coordinates are sampled from within each lesion and a patch is built around that center coordinate. (Default is 20)
#' @param n_models An integer specifying the number of cross-validation models to use for predictions (between 1 and 10). Default is to use all CV models to predict and average the final prediction across all 10 models.
#' @param rotate_patches A logical flag indicating whether to rotate the patches used for predictions. Useful for decreasing dependence of predictions for each sample. (Default is TRUE)
#' @param verbose A logical flag indicating whether to display verbose progress messages. (Default is FALSE)
#'
#' @return A list containing the ALPaCA mask (segmentation of lesions), predictions, and prediction uncertainties.
#'
#' @import ANTsR
#' @import ANTsRCore
#' @import extrantsr
#' @import torch
#'
#' @export
#'
#' @examples \dontrun{
#' # Make predictions using input images and default parameters.
#' predictions <- make_predictions(t1 = t1_image,
#'                                  flair = flair_image,
#'                                  epi = epi_image,
#'                                  phase = phase_image,
#'                                  labeled_candidates = labeled_candidates_image)
#'
#' # Make predictions using input images and return raw probabilities.
#' predictions <- make_predictions(t1 = t1_image,
#'                                  flair = flair_image,
#'                                  epi = epi_image,
#'                                  phase = phase_image,
#'                                  labeled_candidates = labeled_candidates_image,
#'                                  return_raw_probabilities = TRUE)
#' }

make_predictions <- function(ants_list = NULL,
                             t1 = NULL, flair = NULL, epi = NULL, phase = NULL,
                             labeled_candidates = NULL, eroded_candidates = NULL,
                             lesion_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             prl_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             cvs_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             return_raw_probabilities = FALSE, clear_discordant_predictions = TRUE,
                             n_patches = 20, n_models = 10, rotate_patches = TRUE,
                             verbose = FALSE) {
  # Error checking
  if (is.null(ants_list)) { # Make sure all images are provided
    if (any(is.null(t1),
            is.null(flair),
            is.null(epi),
            is.null(phase),
            is.null(labeled_candidates),
            is.null(eroded_candidates))) {
      stop("Images must either be provided via ants_list, or images must be provided for each of t1, flair, epi, phase, labeled_candidates, and eroded_candidates")
    }

    t1 <- check_ants(t1)
    flair <- check_ants(flair)
    epi <- check_ants(epi)
    phase <- check_ants(phase)
    labeled_candidates <- check_ants(labeled_candidates)
    eroded_candidates <- check_ants(eroded_candidates)
  }

  if (!is.null(ants_list)) { # Make sure all images are provided
    if (!all(c("t1", "flair", "epi", "phase", "labeled_candidates", "eroded_candidates") %in% names(ants_list))) {
      stop("If images are provided via ants_list, ants_list must be a named list with items: t1, flair, epi, phase, labeled_candidates, and eroded_candidates. Output from preprocess_images() function can be directly used with return_image = TRUE.")
    }
    t1 <- check_ants(ants_list$t1)
    flair <- check_ants(ants_list$flair)
    epi <- check_ants(ants_list$epi)
    phase <- check_ants(ants_list$phase)
    labeled_candidates <- check_ants(ants_list$labeled_candidates)
    eroded_candidates <- check_ants(ants_list$eroded_candidates)
  }

  if (!all(antsSameMetadata(t1, flair), # Make sure all images are registered
           antsSameMetadata(t1, epi),
           antsSameMetadata(t1, phase),
           antsSameMetadata(t1, labeled_candidates),
           antsSameMetadata(t1, eroded_candidates))) {
    stop("All images provided must be registered to the same space. Please run preprocess_images()")
  }

  # Make sure priorities are understood
  lesion_priority <- match.arg(lesion_priority, c("Youden's J", "Specificity", "Sensitivity"))
  prl_priority <- match.arg(prl_priority, c("Youden's J", "Specificity", "Sensitivity"))
  cvs_priority <- match.arg(cvs_priority, c("Youden's J", "Specificity", "Sensitivity"))

  if (n_patches < 1) {
    stop("n_patches must be a positive integer.")
  }

  if (n_models < 1 | n_models > 10) {
    stop("n_models must be a positive integer between 1 and 10, inclusive.")
  }

  # If there are no lesions, don't have to return anything
  if (sum(labeled_candidates) == 0) {
    warning("No lesion candidates detected.")
    return(NULL)
  }

  # Load CV models
  models_list <- lapply(1:n_models, function(i) {
    return(list(
      jit_load(file.path(system.file(package = "ALPaCA"),
                         "extdata", paste0("autoencoder_", i, ".pt"))),
      jit_load(file.path(system.file(package = "ALPaCA"),
                         "extdata", paste0("predictor_", i, ".pt"))))
    )
  })

  n_lesions <- max(labeled_candidates)
  # Pre-allocate memory
  prediction_tensor <- torch_zeros(c(n_lesions, 3))
  std_tensor <- torch_zeros_like(prediction_tensor)

  if (verbose) {
    print("Running patches through ALPaCA")
  }
  for (candidate_id in 1:n_lesions) {
    if (verbose) {
      print(paste0("Making predictions for lesion ", candidate_id,
                   " of ", n_lesions))
    }
    # Get indexes within lesion indexed by candidate_id
    candidate_coords <- which(labeled_candidates == candidate_id, arr.ind = TRUE)
    under_zero <- apply(candidate_coords - 12, 1, function(i) { # Check if patch bleeds into "nothing"
      any(i < 0)
    })
    over_dim <- apply(candidate_coords + 11, 1, function(i) { # Check if patch bleeds into "nothing" on other side
      any(i[1] > dim(t1)[1],
          i[2] > dim(t1)[2],
          i[3] > dim(t1)[3])
    })
    candidate_coords <- candidate_coords[!under_zero & !over_dim, ] # Check that all patches fit inside the images

    max_coords <- min(n_patches, nrow(candidate_coords)) # Sample some of the candidate_coords
    # If there are no coords with full patch, just guess 0 for everything
    if (max_coords == 0) {
      warning(paste0("No full patches could be extracted for lesion ", candidate_id, ". Default prediction of 0."))
      prediction_tensor[candidate_id, ] <- torch_zeros(c(1, 3))
      std_tensor[candidate_id, ] <- torch_zeros(c(1, 3))
      next
    }

    if (max_coords < n_patches & rotate_patches) { # If we are rotating patches, there is less dependence and resampling same point is still useful
      max_coords <- n_patches
      random_inds <- candidate_coords[sample(1:nrow(candidate_coords),
                                             n_patches, replace = TRUE), ]
    } else {
      random_inds <- candidate_coords[sample(1:nrow(candidate_coords), max_coords), ]
    }
    starts <- random_inds - 12
    ends <- random_inds + 11

    all_patch <- torch_zeros(c(max_coords, 4, 24, 24, 24)) # Pre-allocate memory
    for (patch_id in 1:max_coords) {
      all_patch[patch_id, , , ,] <- extract_patch(candidate_id,  # Extract patches centered at the candidate_coords above. Rotate and flip patches if desired to decrease dependency
                                                  starts[patch_id, ],
                                                  ends[patch_id, ],
                                                  t1, flair, epi, phase,
                                                  labeled_candidates,
                                                  eroded_candidates,
                                                  rotate_patches = rotate_patches)
    }

    all_output <- torch_zeros(c(max_coords * n_models, 3))
    for (model_id in 1:n_models) {
      encoder <- models_list[[model_id]][[1]]$encoder # Extract encoder
      predictor <- models_list[[model_id]][[2]] # Extract predictor

      with_no_grad({ # Run patches through model
        output <- predictor(encoder(all_patch))
      })

      all_output[((model_id - 1) * max_coords + 1):(model_id * max_coords), ] <- output
    }

    prediction_tensor[candidate_id, ] <- torch_mean(all_output, dim = 1, # Get the mean prediction for all coordinates and models
                                                    keepdim = TRUE)
    std_tensor[candidate_id, ] <- torch_std(all_output, dim = 1, # Get the standard deviation of all predictions for a sense of uncertainty
                                            keepdim = TRUE)
  }

  # Convert torch tensor to dataframe
  predictions <- as.data.frame(as.matrix(prediction_tensor))
  std <- as.data.frame(as.matrix(std_tensor))
  names(predictions) <- c("Lesion", "PRL", "CVS")
  names(std) <- c("Lesion", "PRL", "CVS")

  # Convert probability predictions to binary predictions based on thresholds learned from training data
  binary_predictions <- matrix(nrow = nrow(predictions),
                               ncol = ncol(predictions))
  # Lesion thresholding
  if (lesion_priority == "Youden's J") {
    binary_predictions[, 1] <- predictions[, 1] > 0.5517
  }

  if (lesion_priority == "Specificity") {
    binary_predictions[, 1] <- predictions[, 1] > 0.7243
  }

  if (lesion_priority == "Sensitivity") {
    binary_predictions[, 1] <- predictions[, 1] > 0.3787
  }
  #  PRL thresholding
  if (prl_priority == "Youden's J") {
    binary_predictions[, 2] <- predictions[, 2] > 0.0744
  }

  if (prl_priority == "Specificity") {
    binary_predictions[, 2] <- predictions[, 2] > 0.1135
  }

  if (prl_priority == "Sensitivity") {
    binary_predictions[, 2] <- predictions[, 2] > 0.0350
  }
  # CVS thresholding
  if (cvs_priority == "Youden's J") {
    binary_predictions[, 3] <- predictions[, 3] > 0.2094
  }

  if (cvs_priority == "Specificity") {
    binary_predictions[, 3] <- predictions[, 3] > 0.3500
  }

  if (cvs_priority == "Sensitivity") {
    binary_predictions[, 3] <- predictions[, 3] > 0.1102
  }

  binary_predictions[, 2] <- binary_predictions[, 2] * 2 # Allows for unique identification of the row-wise sum
  binary_predictions[, 3] <- binary_predictions[, 3] * 4 # for easier if statements in following section

  lesion_sums <- rowSums(binary_predictions) # 0 is (0, 0, 0), 1 is (1, 0, 0), 3 is (1, 1, 0), 5 is (1, 0, 1), 7 is (1, 1, 1)
  # If sum is 2 or 4, that means prediction is (0, 1, 0) or (0, 0, 1). Since lesion prediction is more reliable than PRL or CVS, convert these to (0, 0, 0)
  if (clear_discordant_predictions) {
    lesion_sums[lesion_sums %% 2 == 0] <- 0
  }

  alpaca_mask <- antsImageClone(labeled_candidates) * 0
  if (return_raw_probabilities) {
    lesion_prob_image <- antsImageClone(alpaca_mask)
    prl_prob_image <- antsImageClone(alpaca_mask)
    cvs_prob_image <- antsImageClone(alpaca_mask)
  }

  n_lesions <- max(labeled_candidates)
  for (j in 1:n_lesions) {
    tmp_lesion_mask <- labeled_candidates == j

    alpaca_mask <- alpaca_mask + tmp_lesion_mask * lesion_sums[j]
    if (return_raw_probabilities) {
      lesion_prob_image <- lesion_prob_image + tmp_lesion_mask * predictions[j, 1]
      prl_prob_image <- prl_prob_image + tmp_lesion_mask * predictions[j, 2]
      cvs_prob_image <- cvs_prob_image + tmp_lesion_mask * predictions[j, 3]
    }
  }

  if (return_raw_probabilities) {
    return(list(alpaca_mask = alpaca_mask,
                raw_probabilities = list(lesion_probs = lesion_prob_image,
                                         prl_probs = prl_prob_image,
                                         cvs_probs = cvs_prob_image),
                predictions = (binary_predictions == 0) * 1, # Convert back to 0s and 1s
                probabilities = predictions,
                prediction_uncertainties = std))
  }

  return(list(alpaca_mask = alpaca_mask,
              predictions = (binary_predictions == 0) * 1, # Convert back to 0s and 1s
              probabilities = predictions,
              prediction_uncertainties = std))
}
