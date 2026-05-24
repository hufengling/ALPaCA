#' Apply a random rigid rotation / reflection to a 3-D patch
#'
#' Rotations are parameterised by three independent draws:
#' * `invert`    — 0 or 1; whether to reflect the patch.
#' * `face`      — 1 – 6; which face of the cube faces "down".
#' * `rotations` — 0 – 3; additional 90-degree rotations around the vertical.
#'
#' Dimension conventions (1-indexed, as required by `torch`):
#'   dim 1 = X, dim 2 = Y, dim 3 = Z  (for the 3-D input tensor).
#'
#' @param patch     `[X, Y, Z]` torch tensor (single-channel 3-D cube).
#' @param invert    Integer scalar (0 or 1).
#' @param face      Integer scalar (1 – 6).
#' @param rotations Integer scalar (0 – 3).
#' @return Rotated `[X, Y, Z]` torch tensor.
#' @noRd
rotate_patch <- function(patch, invert, face, rotations) {
  # Optional reflection along the Y axis
  if (invert == 1L) {
    patch <- torch::torch_flip(patch, dims = 2L)
  }

  # Reorient the cube so a different face points "down"
  if (face == 2L) {
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 2L), 2L, 3L)
  }
  if (face == 3L) {
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 2L), 2L, 3L)
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 2L), 2L, 3L)
  }
  if (face == 4L) {
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 3L), 3L, 2L)
  }
  if (face == 5L) {
    # "Flip cube to left": swap X and Z
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 3L), 1L, 3L)
  }
  if (face == 6L) {
    # "Flip cube to right": swap X and Z in the opposite direction
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 1L), 3L, 1L)
  }

  # In-plane rotation (multiples of 90 degrees in the Y-Z plane)
  if (rotations == 1L) {
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 3L), 3L, 2L) # was flip(2),transp(2,3)
  }
  if (rotations == 2L) {
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 3L), 3L, 2L)
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 3L), 3L, 2L)
  }
  if (rotations == 3L) {
    patch <- torch::torch_transpose(torch::torch_flip(patch, dims = 2L), 2L, 3L)
  }
  patch
}
