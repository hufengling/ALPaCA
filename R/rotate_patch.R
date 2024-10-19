#' Rotate Patch
#'
#' This function rotates a patch based on specified parameters.
#'
#' @param patch A tensor representing the patch to be rotated.
#' @param invert An integer flag indicating whether to invert the patch (1 for inversion, 0 for no inversion).
#' @param face An integer specifying the face of the cube to rotate towards (values from 1 to 6).
#' @param rotations An integer indicating the number of rotations (values from 0 to 3) to apply to the patch.
#'
#' @return A tensor representing the rotated patch.
#'
#' @import torch
#'
#' @examples \dontrun{
#' # Rotate the patch with no inversion, with the second face facedown, and performing one rotation.
#' rotated_patch <- rotate_patch(patch, invert = 0, face = 2, rotations = 1)
#'
#' # Rotate the patch with inversion, with the third face facedown, and performing three rotations.
#' rotated_patch <- rotate_patch(patch, invert = 1, face = 3, rotations = 3)
#' }

rotate_patch <- function(patch, invert, face, rotations) {
  if (invert == 1)
    patch <- torch_flip(patch, 1) # Reflect

  if (face == 2)
    patch <- torch_transpose(torch_flip(patch, 1), 1, 2) # Flip cube towards

  if (face == 3) {
    patch <- torch_transpose(torch_flip(patch, 1), 1, 2) # Flip cube towards twice
    patch <- torch_transpose(torch_flip(patch, 1), 1, 2)
  }

  if (face == 4)
    patch <- torch_transpose(torch_flip(patch, 2), 2, 1) # Flip cube away

  if (face == 5)
    patch <- torch_transpose(torch_flip(patch, 3), 1, 3) # Flip cube to left

  if (face == 6)
    patch <- torch_transpose(torch_flip(patch, 1), 3, 1) # Flip cube to right

  if (rotations == 1)
    patch <- torch_transpose(torch_flip(patch, 2), 2, 3)

  if (rotations == 2) {
    patch <- torch_transpose(torch_flip(patch, 2), 2, 3)
    patch <- torch_transpose(torch_flip(patch, 2), 2, 3)
  }

  if (rotations == 3)
    patch <- torch_transpose(torch_flip(patch, 3), 3, 2)

  return(patch)
}
