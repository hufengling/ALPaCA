

#' Extract and optionally augment a 24 × 24 × 24 patch around one coordinate
#'
#' @param coord         `[4]` torch tensor: (channel, x, y, z) voxel indices.
#' @param candidate_id  Integer. Label of the target candidate.
#' @param t1,flair,phase,epi `[1, X, Y, Z]` torch tensors for each modality.
#' @param lesion_mask        `[1, X, Y, Z]` labeled candidate mask.
#' @param lesion_erode       `[1, X, Y, Z]` eroded candidate mask.
#' @param rotate_patches Logical. Apply random 3-D rotation augmentation?
#' @return `[4, 24, 24, 24]` torch tensor (masked multi-channel patch).
#' @noRd
extract_patch <- function(coord, candidate_id,
                          t1, flair, phase, epi,
                          lesion_mask, lesion_erode,
                          rotate_patches = TRUE) {
  # Convert to R integer (1-indexed); coord is [channel, x, y, z]
  cx <- as.integer(coord[2])
  cy <- as.integer(coord[3])
  cz <- as.integer(coord[4])

  start_ends <- c(
    cx - 12L, cx + 11L,
    cy - 12L, cy + 11L,
    cz - 12L, cz + 11L
  )

  x_start <- max(start_ends[1], 1L)
  x_end <- min(start_ends[2], t1$size(2L))
  y_start <- max(start_ends[3], 1L)
  y_end <- min(start_ends[4], t1$size(3L))
  z_start <- max(start_ends[5], 1L)
  z_end <- min(start_ends[6], t1$size(4L))

  t1_p <- t1[, x_start:x_end, y_start:y_end, z_start:z_end]
  flair_p <- flair[, x_start:x_end, y_start:y_end, z_start:z_end]
  phase_p <- phase[, x_start:x_end, y_start:y_end, z_start:z_end]
  epi_p <- epi[, x_start:x_end, y_start:y_end, z_start:z_end]

  patch <- torch::torch_cat(list(t1_p, flair_p, phase_p, epi_p), dim = 1L)
  lesion_mask_p <- lesion_mask[, x_start:x_end, y_start:y_end, z_start:z_end]
  lesion_erode_p <- lesion_erode[, x_start:x_end, y_start:y_end, z_start:z_end]

  # Pad if the patch clips the volume boundary
  if (!identical(as.integer(patch$size()), c(4L, 24L, 24L, 24L))) {
    padded <- pad_patches(
      patch, lesion_mask_p, lesion_erode_p,
      start_ends, t1
    )
    patch <- padded$patch
    lesion_mask_p <- padded$lesion_mask_patch
    lesion_erode_p <- padded$lesion_erode_patch
  }

  masks <- isolate_lesion(lesion_mask_p, lesion_erode_p, candidate_id)
  isolation_mask <- masks$isolation_mask
  epi_mask <- masks$epi_tensor


  if (rotate_patches) {
    invert <- sample.int(2L, 1L) - 1L # 0 or 1
    face <- sample.int(6L, 1L) # 1 – 6
    rotations <- sample.int(4L, 1L) - 1L # 0 – 3

    patch <- torch::torch_cat(list(
      rotate_patch(patch[1, , , ], invert, face, rotations)$unsqueeze(1L),
      rotate_patch(patch[2, , , ], invert, face, rotations)$unsqueeze(1L),
      rotate_patch(patch[3, , , ], invert, face, rotations)$unsqueeze(1L),
      rotate_patch(patch[4, , , ], invert, face, rotations)$unsqueeze(1L)
    ), dim = 1L)

    isolation_mask <- rotate_patch(isolation_mask[1, , , ], invert, face, rotations)$unsqueeze(1L)
    epi_mask <- rotate_patch(epi_mask[1, , , ], invert, face, rotations)$unsqueeze(1L)
  }

  combined_mask <- torch::torch_cat(
    list(isolation_mask$`repeat`(c(3L, 1L, 1L, 1L)), epi_mask),
    dim = 1L
  )

  patch * combined_mask
}
