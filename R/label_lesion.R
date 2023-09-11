## Label Lesions
#' This function labels lesions in an image based on the provided probability map and a binary mask.
#'
#' @param prob_map A probability map representing the likelihood of lesions.
#' @param bin_map A binary mask identifying lesions above a threshold.
#' @param mincluster The minimum cluster size for lesions (default is 100).
#'
#' @return A labeled image with each lesion identified. Confluent lesions are split
#' @export
#'
#' @examples \dontrun{
#' prob_map <- check_ants("prob_map.nii.gz")
#' mimosa_mask <- check_ants("mimosa_mask.nii.gz")
#' labeled_image <- label_lesion(prob_map, mimosa_mask, mincluster = 100)
#' }
#'
#' @importFrom ANTsRCore labelClusters
#' @importFrom extrantsr ants2oro oro2ants
#'
#' @seealso
#' \code{\link{lesion_centers}} for identifying lesion centers.
#' \code{\link{split_confluent}} for splitting confluent lesions.
#'
#' @description
#' This function takes a probability map and a binary mask as input and labels lesions in the mask
#' based on the probability map and minimum cluster size criteria. It returns a labeled image with
#' lesions identified and labeled.
#'
label_lesion <- function(prob_map, bin_map, mincluster = 100) {
  labeled_image <- labelClusters(bin_map,
                                 minClusterSize = 1,
                                 fullyConnected = TRUE
  )
  size_control <- table(labeled_image[labeled_image != 0])
  if (length(size_control != 0)) {
    size_control <- size_control[size_control > mincluster]
  }
  if (length(size_control) == 0) {
    zero_mask <- bin_map
    zero_mask[zero_mask == 1] <- 0
    return(zero_mask)
  }

  lesion_count <- 1:length(size_control)
  lesion_center_image <- lesion_centers(
    prob_map = prob_map, bin_map = bin_map,
    minCenterSize = mincluster / 10, radius = 1
  )$lesioncenters

  subimg <- lapply(lesion_count, split_confluent,
                   labeled_image = ants2oro(labeled_image),
                   lesion_center_image = ants2oro(lesion_center_image)
  )

  current_lesion_count <- max(subimg[[1]])
  sum_mask <- oro2ants(subimg[[1]])
  for (i in 2:length(subimg)) {
    mask <- subimg[[i]]
    non_zero <- mask[mask > 0]
    mask[mask > 0] <- non_zero + current_lesion_count
    current_lesion_count <- current_lesion_count + max(non_zero)
    sum_mask <- sum_mask + mask
  }
  return(sum_mask)
}
