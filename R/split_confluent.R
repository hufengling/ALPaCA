#' Split Confluent Lesions
#'
#' This function splits confluent lesions labeled as "i" in the labeled image
#' using the lesion center information and returns a labeled image with split lesions.
#'
#' @param i The label of the confluent lesions to be split.
#' @param labeled_image The labeled image containing lesions.
#' @param lesion_center_image The image containing lesion center information.
#'
#' @return A labeled image with split lesions.
#' @export
#'
#' @examples \dontrun{
#' labeled_image <- check_ants("labeled_image.nii.gz")
#' lesion_center_image <- check_ants("lesion_center_image.nii.gz")
#' split_lesion_image <- split_confluent(i = 1, labeled_image, lesion_center_image)
#' }
#'
#' @seealso
#' \code{\link{get_lesion_labels}} for obtaining lesion labels based on lesion centers.
#'
#' @description
#' This function takes a labeled image, lesion center image, and the label "i" of the confluent lesion
#' to be split. It uses the lesion center information to split the confluent lesions and returns a labeled
#' image with the split lesions labeled distinctly.
#'
#' @export
split_confluent <- function(i, labeled_image, lesion_center_image) {
  centers_in_label <- lesion_center_image[labeled_image == i]
  n_centers <- unique(centers_in_label[centers_in_label != 0])

  if (length(n_centers) == 0) {
    return(labeled_image == i)
  }

  lesion <- labeled_image == i
  split_lesion <- get_lesion_labels(
    lesion,
    lesion_center_image * lesion
  )
  s <- unique(split_lesion[split_lesion != 0])
  for (i in 1:length(s)) {
    split_lesion[split_lesion == s[i]] <- i
  }
  return(split_lesion)
}
