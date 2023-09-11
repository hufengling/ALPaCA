#' Extract Patch
#'
#' This function extracts a patch from various MRI images based on specified parameters, 
#' optionally rotates the patches, and returns them as a concatenated tensor.
#'
#' @param candidate_id An identifier for the candidate lesion region.
#' @param patch_starts A vector specifying the starting coordinates of the patch in three dimensions.
#' @param patch_ends A vector specifying the ending coordinates of the patch in three dimensions.
#' @param t1 antsImage representing the T1-weighted MRI image.
#' @param flair antsImage representing the FLAIR MRI image.
#' @param epi antsImage representing the EPI MRI image.
#' @param phase antsImage representing the phase MRI image.
#' @param labeled_candidates antsImage representing labeled candidates for lesion regions.
#' @param rotate_patches A logical flag indicating whether to rotate the extracted patches.
#'
#' @return A concatenated tensor containing the extracted and optionally rotated patches.
#'
#' @export
#' 
#' @import torch
#'
#' @examples \dontrun{
#' # Extract a patch with no rotation.
#' extracted_patch <- extract_patch(candidate_id = 1,
#'                                  patch_starts = c(10, 20, 30),
#'                                  patch_ends = c(20, 30, 40),
#'                                  t1 = t1_image,
#'                                  flair = flair_image,
#'                                  epi = epi_image,
#'                                  phase = phase_image,
#'                                  labeled_candidates = labeled_candidates_image,
#'                                  rotate_patches = FALSE)
#'
#' # Extract a patch and apply random rotation.
#' extracted_patch <- extract_patch(candidate_id = 2,
#'                                  patch_starts = c(15, 25, 35),
#'                                  patch_ends = c(25, 35, 45),
#'                                  t1 = t1_image,
#'                                  flair = flair_image,
#'                                  epi = epi_image,
#'                                  phase = phase_image,
#'                                  labeled_candidates = labeled_candidates_image,
#'                                  rotate_patches = TRUE)
#' }

extract_patch <- function(candidate_id, patch_starts, patch_ends,
                          t1, flair, epi, phase,
                          labeled_candidates,
                          rotate_patches) {
  lesion_mask <- labeled_candidates[patch_starts[1]:patch_ends[1],
                                    patch_starts[2]:patch_ends[2],
                                    patch_starts[3]:patch_ends[3]]
  isolation_mask <- (lesion_mask == 0) + (lesion_mask == candidate_id)
  
  patches <-  list(t1_patch = t1[patch_starts[1]:patch_ends[1],
                                 patch_starts[2]:patch_ends[2],
                                 patch_starts[3]:patch_ends[3]],
                   flair_patch = flair[patch_starts[1]:patch_ends[1],
                                       patch_starts[2]:patch_ends[2],
                                       patch_starts[3]:patch_ends[3]],
                   epi_patch = epi[patch_starts[1]:patch_ends[1],
                                   patch_starts[2]:patch_ends[2],
                                   patch_starts[3]:patch_ends[3]],
                   phase_patch = phase[patch_starts[1]:patch_ends[1],
                                       patch_starts[2]:patch_ends[2],
                                       patch_starts[3]:patch_ends[3]])
  patches <- lapply(patches, function(patch) {
    torch_tensor(patch * isolation_mask)
  })
  
  if (rotate_patches) {
    invert <- sample(0:1, 1) # Mirror patch
    face <- sample(1:6, 1) # Which face of the tensor is "down"
    rotations <- sample(0:3, 1) # Rotate tensor radially once correct face is "down"
    
    patches <- lapply(patches, function(patch_tensor) {
      rotate_patch(patch_tensor, invert, face, rotations)
    })
  }
  
  patches <- lapply(patches, torch_unsqueeze, dim = 1)
  
  return(torch_cat(patches, dim = 1))
}