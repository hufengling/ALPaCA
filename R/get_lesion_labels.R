#' @title Lesion Labeling
#' @description This function is a helper function for lesion_identification(). It takes in a lesion segmentation mask and a NIfTI image with identified lesion centers.
#' @param lesion_mask Lesion mask.
#' @param centers Lesion center map. Provided by lesiontools::lesion_centers().
#'
#' @importFrom Rfast dista
#'
#' @return A NIfTI with each lesion assigned to its closest lesion center
get_lesion_labels <- function(lesion_mask, centers) {
  #### knn on mimosa segmentations ####
  inds_lab <- which(lesion_mask * centers > 0, arr.ind = TRUE) # labeled indices
  inds_cand <- which(lesion_mask == 1 & centers == 0, arr.ind = TRUE) # candidate indices

  pairwise_distances <- Rfast::dista(inds_cand, inds_lab, k = 1, index = TRUE)

  lesion_mask[inds_lab] <- centers[inds_lab]
  lesion_mask[inds_cand] <- centers[inds_lab[pairwise_distances, ]]
  lesion_mask
}
