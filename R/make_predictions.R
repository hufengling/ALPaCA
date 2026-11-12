# ── Internal helpers ──────────────────────────────────────────────────────────

#' Load a NIfTI file as a float torch tensor with a leading channel dimension
#' @noRd
load_nifti_tensor <- function(path) {
  arr <- check_ants(path)
  torch_tensor(as.array(arr), dtype = torch_float())$unsqueeze(1L)
}


#' Build isolation and EPI masks for one candidate
#' @noRd
isolate_lesion <- function(lesion_mask_patch,
                           lesion_erode_patch,
                           candidate_id) {
  tmp_lesion_mask <- lesion_mask_patch$clone()
  tmp_erode       <- lesion_erode_patch$clone()

  lesion_ids <- unique(as.vector(as.array(tmp_lesion_mask)))
  lesion_ids <- lesion_ids[as.numeric(lesion_ids) != 0]

  tmp_mask   <- (tmp_lesion_mask == candidate_id)
  epi_tensor <- torch_clamp(tmp_mask + (tmp_erode == candidate_id), max = 1L)

  isolation_mask <- torch_ones_like(tmp_lesion_mask)
  if (length(lesion_ids) > 1L) {
    isolation_mask <- 0.1 + 0.9 * ((tmp_lesion_mask == 0) + tmp_mask)
  }
  list(isolation_mask = isolation_mask, epi_tensor = epi_tensor)
}


#' Pad a patch that clips the image boundary back to 24 × 24 × 24
#' @noRd

pad_patches <- function(patch, lesion_mask_patch, lesion_erode_patch,
                        start_ends, t1, n_channels = 4L) {

  patch_pad_tensor <- torch_zeros(n_channels, 24L, 24L, 24L)
  mask_pad_tensor  <- torch_zeros(1L, 24L, 24L, 24L)
  erode_pad_tensor <- torch_zeros(1L, 24L, 24L, 24L)

  starts <- start_ends[c(1L, 3L, 5L)]

  # Number of missing voxels on the low side
  start_patch <- ifelse(starts < 1L, 1L - starts, 0L)

  # Size of the actually extracted/clipped patch
  sx <- as.integer(patch$size(2L))
  sy <- as.integer(patch$size(3L))
  sz <- as.integer(patch$size(4L))

  # 1-based destination starting positions
  x0 <- as.integer(start_patch[1] + 1L)
  y0 <- as.integer(start_patch[2] + 1L)
  z0 <- as.integer(start_patch[3] + 1L)

  # narrow() creates a view into the destination tensor.
  # copy_() fills that view from the complete clipped source patch.
  patch_pad_tensor$
    narrow(2L, x0, sx)$
    narrow(3L, y0, sy)$
    narrow(4L, z0, sz)$
    copy_(patch)

  mask_pad_tensor$
    narrow(2L, x0, sx)$
    narrow(3L, y0, sy)$
    narrow(4L, z0, sz)$
    copy_(lesion_mask_patch)

  erode_pad_tensor$
    narrow(2L, x0, sx)$
    narrow(3L, y0, sy)$
    narrow(4L, z0, sz)$
    copy_(lesion_erode_patch)

  list(
    patch = patch_pad_tensor,
    lesion_mask_patch = mask_pad_tensor,
    lesion_erode_patch = erode_pad_tensor
  )
}

# ── Exported function ─────────────────────────────────────────────────────────

#' Make Predictions using ALPaCA
#'
#' This function runs the pre-trained Automated Lesion, PRL (Paramagnetic Rim Lesion), and CVS (Central Vein Sign) Analysis network. This model takes in multi-modal data in the form of pre-processed T1, FLAIR, EPI, and phase images to generate predictions of whether lesion candidates (identified via MIMoSA) are true lesions, PRLs, or CVS.
#'
#' @param ants_list An optional named list containing input images or file paths to images. Names must be: t1, flair, epi, phase, labeled_candidates, and eroded_candidates.
#' @param t1 antsImage or file path to .nii.gz representing T1-weighted MRI image.
#' @param flair antsImage or file path to .nii.gz representing the FLAIR MRI image.
#' @param epi antsImage or file path to .nii.gz representing the EPI MRI image.
#' @param phase antsImage or file path to .nii.gz representing the phase MRI image.
#' @param labeled_candidates antsImage or file path to .nii.gz representing labeled candidates for lesion regions.
#' @param eroded_candidates antsImage or file path to .nii.gz representing eroded candidates for lesion regions.
#' @param output_dir Directory where results will be saved.
#' @param lesion_priority Threshold priority for lesion prediction -- Youden's J, Specificity, Sensitivity.
#' @param prl_priority Threshold priority for PRL prediction.
#' @param cvs_priority Threshold priority for CVS prediction.
#' @param return_raw_probabilities Logical. Whether to return raw probability antsImages. (Default FALSE)
#' @param clear_discordant_predictions Logical. Whether to clear discordant predictions. (Default TRUE)
#' @param n_patches Integer. Number of patches to sample per candidate. (Default 25)
#' @param model_ids Vector of CV models to use. (Default 1:10)
#' @param rotate_patches Logical. Whether to rotate patches. (Default TRUE)
#' @param verbose Logical. Whether to print progress. (Default FALSE)
#'
#' @return A list containing the ALPaCA mask, predictions, and prediction uncertainties.
#'
#' @importFrom ANTsR antsImageRead antsImageWrite antsImageClone
#' @importFrom utils write.csv
#' @importFrom extrantsr check_ants
#' @import torch
#'
#' @export
make_predictions <- function(ants_list = NULL,
                             t1 = NULL, flair = NULL, epi = NULL, phase = NULL,
                             labeled_candidates = NULL,
                             eroded_candidates = NULL,
                             output_dir = NULL,
                             lesion_priority = c("Youden's J", "Specificity", "Sensitivity"),
                             prl_priority    = c("Youden's J", "Specificity", "Sensitivity"),
                             cvs_priority    = c("Youden's J", "Specificity", "Sensitivity"),
                             return_raw_probabilities     = FALSE,
                             clear_discordant_predictions = TRUE,
                             n_patches      = 25,
                             model_ids      = 1:10,
                             rotate_patches = TRUE,
                             verbose        = FALSE) {

  # ── Input validation ───────────────────────────────────────────────────────
  if (n_patches < 1) {
    stop("n_patches must be a positive integer.")
  }
  if (any(model_ids < 1) || any(model_ids > 10)) {
    stop("model_ids must be a vector of integers between 1 and 10, inclusive.")
  }
  if (is.null(output_dir) || !dir.exists(output_dir)) {
    stop("`output_dir` must be a valid, existing directory path.")
  }

  lesion_priority <- match.arg(lesion_priority, c("Youden's J", "Specificity", "Sensitivity"))
  prl_priority    <- match.arg(prl_priority,    c("Youden's J", "Specificity", "Sensitivity"))
  cvs_priority    <- match.arg(cvs_priority,    c("Youden's J", "Specificity", "Sensitivity"))

  device <- if (cuda_is_available()) "cuda" else "cpu"
  device <- torch_device(device)
  message("Using device: ", device)

  # ── Load images ────────────────────────────────────────────────────────────
  # extract_patch requires torch tensors ([1, X, Y, Z]); it indexes 4D and calls
  # t1$size().  We keep an ants clone of labeled_candidates only for building the
  # output images / writing NIfTIs.
  message("Loading images")

  if (!is.null(ants_list)) {
    if (!all(c("t1", "flair", "epi", "phase", "labeled_candidates",
               "eroded_candidates") %in% names(ants_list))) {
      stop("ants_list must be a named list with items: t1, flair, epi, phase, labeled_candidates, and eroded_candidates.")
    }
    t1_src <- ants_list$t1; flair_src <- ants_list$flair
    epi_src <- ants_list$epi; phase_src <- ants_list$phase
    lab_src <- ants_list$labeled_candidates
    ero_src <- ants_list$eroded_candidates
  } else {
    if (any(is.null(t1), is.null(flair), is.null(epi), is.null(phase),
            is.null(labeled_candidates), is.null(eroded_candidates))) {
      stop("Images must either be provided via ants_list, or individually for each of t1, flair, epi, phase, labeled_candidates, and eroded_candidates.")
    }
    t1_src <- t1; flair_src <- flair
    epi_src <- epi; phase_src <- phase
    lab_src <- labeled_candidates
    ero_src <- eroded_candidates
  }

  # Torch tensors for inference
  t1_t    <- load_nifti_tensor(t1_src)
  flair_t <- load_nifti_tensor(flair_src)
  epi_t   <- load_nifti_tensor(epi_src)
  phase_t <- load_nifti_tensor(phase_src)
  lab_t   <- load_nifti_tensor(lab_src)
  ero_t   <- load_nifti_tensor(ero_src)

  # Ants clone for output image assembly
  labeled_candidates_ants <- check_ants(lab_src)

  if (as.numeric(lab_t$sum()$item()) == 0) {
    warning("No lesion candidates detected.")
    return(NULL)
  }

  # ── Load models ────────────────────────────────────────────────────────────
  model_dir <- system.file("extdata", package = "ALPaCA")
  if (verbose) message("Loading models from: ", model_dir)

  autoencoder_paths <- lapply(model_ids, \(i)
    file.path(model_dir, paste0("autoencoder_", i, ".pt")))
  predictor_paths <- lapply(model_ids, \(i)
    file.path(model_dir, paste0("predictor_", i, ".pt")))

  # ── Per-candidate inference ────────────────────────────────────────────────
  n_candidates <- as.integer(lab_t$max()$item())
  n_classes    <- 3L

  if (verbose) message("Making predictions for ", n_candidates, " lesion candidates")

  output_tensor   <- torch_zeros(n_candidates, n_classes)
  variance_tensor <- torch_zeros(n_candidates, n_classes)

  for (candidate_id in seq_len(n_candidates)) {
    if (verbose) message("  Candidate ", candidate_id, " / ", n_candidates)

    # torch_nonzero on a [1, X, Y, Z] tensor returns 1-based [ch, x, y, z] rows,
    # which is exactly what extract_patch's `coord` argument expects.
    coords <- torch_nonzero(lab_t == candidate_id)
    n_avail <- coords$size(1L)
    if (n_avail == 0L) next

    n_actual <- min(n_patches, n_avail)
    sel      <- sample.int(n_avail, n_actual)

    all_patch <- torch_zeros(n_actual, 4L, 24L, 24L, 24L)

    for (i in seq_len(n_actual)) {
      coord <- coords[sel[i], ]

      # extract_patch handles bounds-clamping, padding, isolation and rotation
      # internally.  Note the arg order: coord, candidate_id, t1, flair, phase,
      # epi, lesion_mask, lesion_erode, rotate_patches.
      patch <- extract_patch(
        coord, candidate_id,
        t1_t, flair_t, phase_t, epi_t,
        lab_t, ero_t,
        rotate_patches = rotate_patches
      )

      all_patch[i, , , , ] <- patch
    }

    all_patch      <- all_patch$to(device = device)
    prediction_sum <- torch_zeros(n_actual, n_classes)

    for (model_ind in seq_along(model_ids)) {
      prl_autoencoder <- jit_load(autoencoder_paths[[model_ind]][1], device = "cpu")
      prl_autoencoder$to(device = device)

      prl_predictor <- jit_load(predictor_paths[[model_ind]][1], device = "cpu")
      prl_predictor$to(device = device)

      with_no_grad({
        encoded    <- prl_autoencoder$encoder(all_patch)
        raw_output <- prl_predictor(encoded)
      })

      prediction_sum <- prediction_sum + raw_output$to(device = "cpu")

      rm(prl_autoencoder, prl_predictor)
      gc()
    }

    prediction_mean <- prediction_sum / length(model_ids)

    if (n_actual > 1L) {
      output_tensor[candidate_id, ]   <- torch_mean(prediction_mean, dim = 1L)
      variance_tensor[candidate_id, ] <- torch_var(prediction_mean,  dim = 1L)
    } else {
      output_tensor[candidate_id, ]   <- prediction_mean$squeeze(1L)
      variance_tensor[candidate_id, ] <- torch_zeros(n_classes)
    }
  }

  # ── Build output data.frame ────────────────────────────────────────────────
  combined  <- torch_cat(list(output_tensor, variance_tensor), dim = 2L)
  output_df <- as.data.frame(as.matrix(combined$detach()))
  names(output_df) <- c("lesion_prob", "PRL_prob", "CVS_prob",
                        "lesion_prob_variance", "PRL_prob_variance", "CVS_prob_variance")

  if (verbose) message("Assembling predictions")

  binary_predictions <- matrix(
    nrow = nrow(output_df), ncol = 3,
    dimnames = list(NULL, c("lesion", "PRL", "CVS"))
  )

  # Lesion thresholding
  if (lesion_priority == "Youden's J")  binary_predictions[, 1] <- output_df[, 1] > 0.5517
  if (lesion_priority == "Specificity") binary_predictions[, 1] <- output_df[, 1] > 0.7243
  if (lesion_priority == "Sensitivity") binary_predictions[, 1] <- output_df[, 1] > 0.3787

  # PRL thresholding
  if (prl_priority == "Youden's J")  binary_predictions[, 2] <- output_df[, 2] > 0.0744
  if (prl_priority == "Specificity") binary_predictions[, 2] <- output_df[, 2] > 0.1135
  if (prl_priority == "Sensitivity") binary_predictions[, 2] <- output_df[, 2] > 0.0350

  # CVS thresholding
  if (cvs_priority == "Youden's J")  binary_predictions[, 3] <- output_df[, 3] > 0.2094
  if (cvs_priority == "Specificity") binary_predictions[, 3] <- output_df[, 3] > 0.3500
  if (cvs_priority == "Sensitivity") binary_predictions[, 3] <- output_df[, 3] > 0.1102

  # Encode columns for unique row-sum identification:
  # 0=(0,0,0), 1=(1,0,0), 3=(1,1,0), 5=(1,0,1), 7=(1,1,1)
  binary_predictions[, 2] <- binary_predictions[, 2] * 2
  binary_predictions[, 3] <- binary_predictions[, 3] * 4

  lesion_sums <- rowSums(binary_predictions)

  if (clear_discordant_predictions) {
    lesion_sums[lesion_sums %% 2 == 0] <- 0
    binary_predictions[lesion_sums == 0, ] <- 0
  }

  # ── Build output images ────────────────────────────────────────────────────
  alpaca_mask <- antsImageClone(labeled_candidates_ants) * 0

  if (return_raw_probabilities) {
    lesion_prob_image <- antsImageClone(alpaca_mask)
    prl_prob_image    <- antsImageClone(alpaca_mask)
    cvs_prob_image    <- antsImageClone(alpaca_mask)
  }

  n_lesions <- max(labeled_candidates_ants)
  for (j in seq_len(n_lesions)) {
    tmp_lesion_mask <- labeled_candidates_ants == j
    alpaca_mask     <- alpaca_mask + tmp_lesion_mask * lesion_sums[j]

    if (return_raw_probabilities) {
      lesion_prob_image <- lesion_prob_image + tmp_lesion_mask * output_df[j, 1]
      prl_prob_image    <- prl_prob_image    + tmp_lesion_mask * output_df[j, 2]
      cvs_prob_image    <- cvs_prob_image    + tmp_lesion_mask * output_df[j, 3]
    }
  }

  # ── Write outputs ──────────────────────────────────────────────────────────
  antsImageWrite(alpaca_mask, file.path(output_dir, "alpaca_mask.nii.gz"))

  if (return_raw_probabilities) {
    antsImageWrite(lesion_prob_image, file.path(output_dir, "lesion_prob.nii.gz"))
    antsImageWrite(prl_prob_image,    file.path(output_dir, "prl_prob.nii.gz"))
    antsImageWrite(cvs_prob_image,    file.path(output_dir, "cvs_prob.nii.gz"))
  }

  write.csv((binary_predictions != 0) * 1, file.path(output_dir, "predictions.csv"))
  write.csv(output_df,                     file.path(output_dir, "probabilities.csv"))

  # ── Return ─────────────────────────────────────────────────────────────────
  result <- list(
    alpaca_mask   = alpaca_mask,
    predictions   = (binary_predictions != 0) * 1,
    probabilities = output_df
  )

  if (return_raw_probabilities) {
    result$raw_probabilities <- list(
      lesion_probs = lesion_prob_image,
      prl_probs    = prl_prob_image,
      cvs_probs    = cvs_prob_image
    )
  }

  result
}
