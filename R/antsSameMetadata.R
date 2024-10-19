#' Check if Two ANTs Images Have the Same Metadata
#'
#' This function checks if two ANTs images have the same metadata, including direction, origin, and spacing.
#'
#' @param ants1 An ANTs image object.
#' @param ants2 Another ANTs image object.
#' @param tolerance Tolerance for differences in metadata
#'
#' @return A logical value indicating whether the two ANTs images have the same metadata. Returns TRUE if they match and FALSE otherwise.
#'
#' @importFrom ANTsRCore antsGetDirection antsGetOrigin antsGetSpacing
#' @examples \dontrun{
#' antsImage1 <- antsImageRead("image1.nii")
#' antsImage2 <- antsImageRead("image2.nii")
#' antsSameMetadata(antsImage1, antsImage2)
#' }
#'
#' @seealso
#' \code{\link{antsGetDirection}}, \code{\link{antsGetOrigin}}, \code{\link{antsGetSpacing}}
#'
antsSameMetadata <- function(ants1, ants2, tolerance = 0.000001) {
  same_direction <- all(abs(antsGetDirection(ants1) - antsGetDirection(ants2)) > tolerance)
  same_origin <- all(antsGetOrigin(ants1) - antsGetOrigin(ants2) > tolerance)
  same_spacing <- all(antsGetSpacing(ants1) - antsGetSpacing(ants2) > tolerance)
  return(same_direction & same_origin & same_spacing)
}
