# ── Internal helpers ──────────────────────────────────────────────────────────

#' Load a NIfTI file as a float torch tensor with a leading channel dimension
#'
#' @param path Character. Path to a `.nii` or `.nii.gz` file.
#' @return A `[1, X, Y, Z]` float32 torch tensor.
#' @noRd
load_nifti_tensor <- function(path) {
  arr <- antsImageRead(path)
  torch_tensor(as.array(arr), dtype = torch_float())$unsqueeze(1L)
}


#' Sample voxel coordinates for a single lesion candidate
#'
#' @param candidate_id Integer. Label value of the candidate in `lesion_mask`.
#' @param num_coords Integer. Maximum number of coordinates to sample.
#' @param lesion_mask `[1, X, Y, Z]` torch tensor. Labeled candidate mask.
#' @return A `[k, 4]` torch tensor of (channel, x, y, z) indices, where
#'   `k <= num_coords`.
#' @noRd
get_coords <- function(candidate_id, num_coords, lesion_mask) {
  if (as.numeric(sum(lesion_mask == candidate_id)) == 0) {
    return(NULL)
  }
  candidate_coords <- torch_nonzero(lesion_mask == candidate_id)
  n <- candidate_coords$size(1L)
  max_coords <- min(num_coords, n)
  random_inds <- sample.int(n, max_coords)
  candidate_coords[random_inds, ]
}


#' Build isolation and EPI masks for one candidate
#'
#' Scales down overlapping neighbours to 10 % of their intensity so the target
#' lesion dominates the patch, and constructs a binary EPI mask from the eroded
#' candidate mask.
#'
#' @param lesion_mask_patch  `[1, 24, 24, 24]` torch tensor.
#' @param lesion_erode_patch `[1, 24, 24, 24]` torch tensor.
#' @param candidate_id Integer. Target candidate label.
#' @return Named list with elements `isolation_mask` and `epi_tensor`, both
#'   `[1, 24, 24, 24]` torch tensors.
#' @noRd
isolate_lesion <- function(lesion_mask_patch,
                           lesion_erode_patch,
                           candidate_id) {
  tmp_lesion_mask <- lesion_mask_patch$clone()
  tmp_erode <- lesion_erode_patch$clone()

  lesion_ids <- unique(as.vector(as.array(tmp_lesion_mask)))
  lesion_ids <- lesion_ids[as.numeric(lesion_ids) != 0]

  tmp_mask <- (tmp_lesion_mask == candidate_id)
  epi_tensor <- tmp_mask + (tmp_erode == candidate_id)

  isolation_mask <- torch_ones_like(tmp_lesion_mask)
  if (length(lesion_ids) > 1L) {
    isolation_mask <- 0.1 + 0.9 * ((tmp_lesion_mask == 0) + tmp_mask)
  }
  list(isolation_mask = isolation_mask, epi_tensor = epi_tensor)
}

#' Pad a patch that clips the image boundary back to 24 × 24 × 24
#'
#' @param patch              `[C, x, y, z]` torch tensor (possibly < 24 on each side).
#' @param lesion_mask_patch  `[1, x, y, z]` torch tensor.
#' @param lesion_erode_patch `[1, x, y, z]` torch tensor.
#' @param start_ends         Integer vector of length 6:
#'   `c(x_start, x_end, y_start, y_end, z_start, z_end)` before boundary clamping.
#' @param t1                 `[1, X, Y, Z]` reference tensor (full volume).
#' @param n_channels         Integer. Number of image channels (default 4).
#' @return Named list: `patch`, `lesion_mask_patch`, `lesion_erode_patch`.
#' @noRd
pad_patches <- function(patch, lesion_mask_patch, lesion_erode_patch,
                        start_ends, t1, n_channels = 4L) {
  patch_pad_tensor <- torch_zeros(n_channels, 24L, 24L, 24L)
  mask_pad_tensor <- torch_zeros(1L, 24L, 24L, 24L)
  erode_pad_tensor <- torch_zeros(1L, 24L, 24L, 24L)

  starts <- start_ends[c(1L, 3L, 5L)]
  start_patch <- ifelse(starts < 0L, -starts, 0L)

  ends <- start_ends[c(2L, 4L, 6L)]
  vol_sizes <- c(t1$size(2L), t1$size(3L), t1$size(4L))
  end_patch <- ifelse(ends > vol_sizes, 24L - (ends - vol_sizes), 24L)

  patch_pad_tensor[
    ,
    (start_patch[1] + 1L):end_patch[1],
    (start_patch[2] + 1L):end_patch[2],
    (start_patch[3] + 1L):end_patch[3]
  ] <- patch

  mask_pad_tensor[
    ,
    (start_patch[1] + 1L):end_patch[1],
    (start_patch[2] + 1L):end_patch[2],
    (start_patch[3] + 1L):end_patch[3]
  ] <- lesion_mask_patch

  erode_pad_tensor[
    ,
    (start_patch[1] + 1L):end_patch[1],
    (start_patch[2] + 1L):end_patch[2],
    (start_patch[3] + 1L):end_patch[3]
  ] <- lesion_erode_patch

  list(
    patch = patch_pad_tensor,
    lesion_mask_patch = mask_pad_tensor,
    lesion_erode_patch = erode_pad_tensor
  )
}

# ── Exported function ─────────────────────────────────────────────────────────

#' #' #' Make Predictions using ALPaCA
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
#' @param output_dir Directory where results will be saved.
#' @param lesion_priority A character vector specifying priority for lesion prediction thresholds -- Youden's J, Specificity, Sensitivity. Thresholds are based on training set ROC curves from CV models. Default priority is Youden's J, with sensitivity \eqn{\approx} 0.83 and specificity \eqn{\approx} 0.86. 'Specificity' prioritizes specificity 3 times more than sensitivity, with sensitivity \eqn{\approx} 0.69 and specificity \eqn{\approx} 0.94. 'Sensitivity' prioritizes sensitivity 3 times more than specificity, with sensitivity \eqn{\approx} 0.92 and specificity \eqn{\approx} 0.70.
#' @param prl_priority A character vector specifying priority for PRL prediction thresholds. Same options and default as lesion_priority. For Youden's J, sensitivity \eqn{\approx} 0.76 and specificity \eqn{\approx} 0.83. For Specificity, sensitivity \eqn{\approx} 0.63 and specificity \eqn{\approx} 0.90. For Sensitivity, sensitivity \eqn{\approx} 0.86 and specificity \eqn{\approx} 0.64.
#' @param cvs_priority A character vector specifying priority for CVS prediction thresholds.Same options and default as lesion_priority. For Youden's J, sensitivity \eqn{\approx} 0.81 and specificity \eqn{\approx} 0.65. For Specificity, sensitivity \eqn{\approx} 0.27 and specificity \eqn{\approx} 0.91. For Sensitivity, sensitivity \eqn{\approx} 0.91 and specificity \eqn{\approx} 0.47.
#' @param return_raw_probabilities A logical flag indicating whether to return raw probability antsImages for each region. A raw probability lesion-wise dataframe is always returned. (Default is FALSE)
#' @param clear_discordant_predictions A logical flag indicating whether to clear discordant predictions (ie candidates where the model predicts "CVS"/"PRL, but not "lesion".) In training, lesion prediction was almost always more reliable under Youden's J threshold choice. (Default is TRUE)
#' @param n_patches An integer specifying the number of patches to sample for predictions. Coordinates are sampled from within each lesion and a patch is built around that center coordinate. (Default is 20)
#' @param model_ids Vector of CV models to use with default to use all CV models to predict and average the final prediction across all 10 models (Default 1:10 is all 10 models)
#' @param rotate_patches A logical flag indicating whether to rotate the patches used for predictions. Useful for decreasing dependence of predictions for each sample. (Default is TRUE)
#' @param verbose A logical flag indicating whether to display verbose progress messages. (Default is FALSE)
#'
#' @return A list containing the ALPaCA mask (segmentation of lesions), predictions, and prediction uncertainties.
#'
#' @importFrom ANTsR antsImageRead antsImageWrite antsImageClone
#' @importFrom utils write.csv
#' @import torch
#'
#' @export
#'
#' @examples \dontrun{
#' # Make predictions using input images and default parameters.
#' predictions <- make_predictions(
#'   t1 = t1_image,
#'   flair = flair_image,
#'   epi = epi_image,
#'   phase = phase_image,
#'   labeled_candidates = labeled_candidates_image
#' )
#'
#' # Make predictions using input images and return raw probabilities.
#' predictions <- make_predictions(
#'   t1 = t1_image,
#'   flair = flair_image,
#'   epi = epi_image,
#'   phase = phase_image,
#'   labeled_candidates = labeled_candidates_image,
#'   return_raw_probabilities = TRUE
#' )
#' }
make_predictions <- function(ants_list = NULL,
                             t1 = NULL, flair = NULL, epi = NULL, phase = NULL,
                             labeled_candidates = NULL,
                             eroded_candidates = NULL,
                             output_dir = NULL,
                             lesion_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             prl_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             cvs_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             return_raw_probabilities = FALSE,
                             clear_discordant_predictions = TRUE,
                             n_patches = 25,
                             model_ids = 1:10,
                             rotate_patches = TRUE,
                             verbose = FALSE) {
  if (n_patches < 1) {
    stop("n_patches must be a positive integer.")
  }

  if (any(model_ids) < 1 || any(model_ids) > 10) {
    stop("model_ids must be either a positive integer between 1 and 10, inclusive, or a vector of integers between 1 and 10.")
  }

  if (is.null(output_dir)) {
    stop("Must provide output_dir to write niftis and csvs")
  }

  # Make sure priorities are understood
  lesion_priority <- match.arg(lesion_priority, c("Youden's J", "Specificity", "Sensitivity"))
  prl_priority <- match.arg(prl_priority, c("Youden's J", "Specificity", "Sensitivity"))
  cvs_priority <- match.arg(cvs_priority, c("Youden's J", "Specificity", "Sensitivity"))

  device <- if (cuda_is_available()) "cuda" else "cpu"
  device <- torch_device(device)
  message("Using device: ", device)

  # ── Load images ─────────────────────────────────────────────────────────────
  message("Loading images")

  # Error checking
  if (!is.null(ants_list)) { # Make sure all images are provided
    if (!all(c("t1", "flair", "epi", "phase", "labeled_candidates", "eroded_candidates") %in% names(ants_list))) {
      stop("If images are provided via ants_list, ants_list must be a named list with items: t1, flair, epi, phase, labeled_candidates, and eroded_candidates. Output from preprocess_images() function can be directly used with return_image = TRUE.")
    }
    t1 <- load_nifti_tensor(ants_list$t1)
    flair <- load_nifti_tensor(ants_list$flair)
    epi <- load_nifti_tensor(ants_list$epi)
    phase <- load_nifti_tensor(ants_list$phase)

    labeled_candidates_ants <- antsImageRead(ants_list$labeled_candidates)
    labeled_candidates <- load_nifti_tensor(ants_list$labeled_candidates)
    eroded_candidates <- load_nifti_tensor(ants_list$eroded_candidates)
  }
  if (is.null(ants_list)) {
    if (any(
      is.null(t1),
      is.null(flair),
      is.null(epi),
      is.null(phase),
      is.null(labeled_candidates),
      is.null(eroded_candidates)
    )) {
      stop("Images must either be provided via ants_list, or images must be provided for each of t1, flair, epi, phase, labeled_candidates, and eroded_candidates")
    }
    t1 <- load_nifti_tensor(t1)
    flair <- load_nifti_tensor(flair)
    epi <- load_nifti_tensor(epi)
    phase <- load_nifti_tensor(phase)

    labeled_candidates_ants <- antsImageRead(labeled_candidates)
    labeled_candidates <- load_nifti_tensor(labeled_candidates)
    eroded_candidates <- load_nifti_tensor(eroded_candidates)
  }

  # If there are no lesions, don't have to return anything
  if (as.numeric(sum(labeled_candidates)) == 0) {
    warning("No lesion candidates detected.")
    return(NULL)
  }

  # ── Load models ─────────────────────────────────────────────────────────────
  model_dir <- system.file("extdata", package = "ALPaCA")

  if (verbose)
    message("Loading models from: ", model_dir)

  autoencoder_paths <- lapply(model_ids, \(i) {
    file.path(model_dir, paste0("autoencoder_", i, ".pt"))
  })

  predictor_paths <- lapply(model_ids, \(i) {
    file.path(model_dir, paste0("predictor_", i, ".pt"))
  })

  # ── Per-candidate inference ─────────────────────────────────────────────────
  n_candidates <- as.integer(labeled_candidates$max()$item())
  n_classes <- 3L

  if (verbose)
    message("Making predictions for ", n_candidates, " lesion candidates")

  output_tensor <- torch_zeros(n_candidates, n_classes)
  variance_tensor <- torch_zeros(n_candidates, n_classes)

  for (candidate_id in seq_len(n_candidates)) {
    if (verbose)
      message("  Candidate ", candidate_id, " / ", n_candidates)
    coords <- get_coords(candidate_id, n_patches, labeled_candidates)

    if (is.null(coords)) {
      next
    }

    if (coords$size(1L) > 0L) {
      all_patch <- torch_zeros(
        coords$size(1L), 4L,
        24, 24, 24
      )

      for (i in seq_len(coords$size(1L))) {
        coord <- coords[i, ]
        patch <- extract_patch(
          coord, candidate_id,
          t1, flair, phase, epi,
          labeled_candidates, eroded_candidates,
          rotate_patches = rotate_patches
        )
        all_patch[i, , , , ] <- patch
      }

      all_patch <- all_patch$to(device = device)

      prediction_sum <- torch_zeros(nrow(coords), 3)
      for (model_ind in seq_along(model_ids)) {
        prl_autoencoder <- jit_load(
          autoencoder_paths[[model_ind]][1],
          device = "cpu"
        )
        prl_autoencoder$to(device = device)

        prl_predictor <- jit_load(
          predictor_paths[[model_ind]][1],
          device = "cpu"
        )
        prl_predictor$to(device = device)

        with_no_grad({
          encoded <- prl_autoencoder$encoder(all_patch)
          raw_output <- prl_predictor(encoded)
        })
        prediction <- raw_output$to(device = "cpu")
        prediction_sum <- prediction_sum + prediction

        rm(prl_autoencoder)
        rm(prl_predictor)
        gc()
      }
      prediction_mean <- prediction_sum / length(model_ids)

      if (n_patches > 1L) {
        output_tensor[candidate_id, ] <- torch_mean(prediction_mean, dim = 1L)
        variance_tensor[candidate_id, ] <- torch_var(prediction_mean, dim = 1L)
      } else {
        output_tensor[candidate_id, ] <- prediction_mean
        variance_tensor[candidate_id, ] <- torch_zeros(1L, n_classes)
      }
    }
    # Candidates with no sampled coordinates retain their initialised zeros
  }

  # ── Build output data.frame ─────────────────────────────────────────────────
  combined <- torch_cat(
    list(output_tensor, variance_tensor),
    dim = 2L
  )
  output_df <- as.data.frame(as.matrix(combined$detach()))

  names(output_df) <- c("lesion_prob", "PRL_prob", "CVS_prob",
                        "lesion_prob_variance", "PRL_prob_variance",
                        "CVS_prob_variance")

  if (verbose)
    message("Assembling predictions")
  # Convert probability predictions to binary predictions based on thresholds learned from training data
  binary_predictions <- matrix(
    nrow = nrow(output_df),
    ncol = 3
  )
  # Lesion thresholding
  if (lesion_priority == "Youden's J") {
    binary_predictions[, 1] <- output_df[, 1] > 0.5517
  }

  if (lesion_priority == "Specificity") {
    binary_predictions[, 1] <- output_df[, 1] > 0.7243
  }

  if (lesion_priority == "Sensitivity") {
    binary_predictions[, 1] <- output_df[, 1] > 0.3787
  }
  #  PRL thresholding
  if (prl_priority == "Youden's J") {
    binary_predictions[, 2] <- output_df[, 2] > 0.0744
  }

  if (prl_priority == "Specificity") {
    binary_predictions[, 2] <- output_df[, 2] > 0.1135
  }

  if (prl_priority == "Sensitivity") {
    binary_predictions[, 2] <- output_df[, 2] > 0.0350
  }
  # CVS thresholding
  if (cvs_priority == "Youden's J") {
    binary_predictions[, 3] <- output_df[, 3] > 0.2094
  }

  if (cvs_priority == "Specificity") {
    binary_predictions[, 3] <- output_df[, 3] > 0.3500
  }

  if (cvs_priority == "Sensitivity") {
    binary_predictions[, 3] <- output_df[, 3] > 0.1102
  }

  binary_predictions[, 2] <- binary_predictions[, 2] * 2 # Allows for unique identification of the row-wise sum
  binary_predictions[, 3] <- binary_predictions[, 3] * 4 # for easier if statements in following section

  lesion_sums <- rowSums(binary_predictions) # 0 is (0, 0, 0), 1 is (1, 0, 0), 3 is (1, 1, 0), 5 is (1, 0, 1), 7 is (1, 1, 1)
  # If sum is 2, 4, or 6, that means prediction is (0, 1, 0), (0, 0, 1), or (0, 1, 1). Since lesion prediction is more reliable than PRL or CVS, convert these to (0, 0, 0)
  if (clear_discordant_predictions) {
    lesion_sums[lesion_sums %% 2 == 0] <- 0
    binary_predictions[lesion_sums == 0, 1] <- 0
    binary_predictions[lesion_sums == 0, 2] <- 0
    binary_predictions[lesion_sums == 0, 3] <- 0
  }

  alpaca_mask <- antsImageClone(labeled_candidates_ants) * 0
  if (return_raw_probabilities) {
    lesion_prob_image <- antsImageClone(alpaca_mask)
    prl_prob_image <- antsImageClone(alpaca_mask)
    cvs_prob_image <- antsImageClone(alpaca_mask)
  }

  n_lesions <- max(labeled_candidates_ants)
  for (j in 1:n_lesions) {
    tmp_lesion_mask <- labeled_candidates_ants == j

    alpaca_mask <- alpaca_mask + tmp_lesion_mask * lesion_sums[j]
    if (return_raw_probabilities) {
      lesion_prob_image <- lesion_prob_image + tmp_lesion_mask * output_df[j, 1]
      prl_prob_image <- prl_prob_image + tmp_lesion_mask * output_df[j, 2]
      cvs_prob_image <- cvs_prob_image + tmp_lesion_mask * output_df[j, 3]
    }
  }
  # saving alpaca_mask, raw probability images, and prediction outputs
  antsImageWrite(alpaca_mask, file.path(output_dir, "alpaca_mask.nii.gz"))

  if (return_raw_probabilities) {
    antsImageWrite(lesion_prob_image,
                   file.path(output_dir, "lesion_prob.nii.gz"))
    antsImageWrite(prl_prob_image, file.path(output_dir, "prl_prob.nii.gz"))
    antsImageWrite(cvs_prob_image, file.path(output_dir, "cvs_prob.nii.gz"))
  }

  write.csv((binary_predictions != 0) * 1,
            file.path(output_dir, "predictions.csv"))
  write.csv(output_df,
            file.path(output_dir, "probabilities.csv")) # Assuming probabilities are the same here

  if (return_raw_probabilities) {
    return(list(
      alpaca_mask = alpaca_mask,
      raw_probabilities = list(
        lesion_probs = lesion_prob_image,
        prl_probs = prl_prob_image,
        cvs_probs = cvs_prob_image
      ),
      predictions = (binary_predictions != 0) * 1, # Convert back to 0s and 1s
      probabilities = output_df
    ))
  }

  list(
    alpaca_mask = alpaca_mask,
    predictions = (binary_predictions != 0) * 1, # Convert back to 0s and 1s
    probabilities = output_df
  )
}
