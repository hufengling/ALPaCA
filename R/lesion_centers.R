#' @title Distinct Lesion Centers
#' @description This function finds the centers of distinct lesions based on a lesion probability map. Each lesion center is labeled via a unique integer. The method is described in Dworkin et al., (2018).
#' @param prob_map an image of class \code{antsImage} containing the probability that each voxel is a lesion voxel
#' @param bin_map a binarized image of class \code{antsImage} in which voxels are classified as either lesion voxels or not lesion voxels.
#' @param minCenterSize an integer value representing the minimum number of connected voxels that can be considered a lesion center
#' @param radius an integer specifying radius of the neighborhood (in voxels) for which the hessian should be calculated.
#'
#' @importFrom ANTsRCore labelClusters
#' @importFrom extrantsr ants2oro oro2ants
#' @return A list containing lesioncenters (antsImage object with labeled lesion centers) and lesioncount (an integer value representing the number of distinct lesions)
#' @examples \dontrun{
#' library(neurobase)
#' lesion.probs <- readnii("path/to/probabilitymap")
#' centers <- lesioncenters(
#'   probmap = lesion.probs, binmap = lesion.probs > 0.30
#' )
#' }
#' @export
#' @references J.D. Dworkin, K.A. Linn, I. Oguz, G.M. Fleishman, R. Bakshi, G. Nair, P.A. Calabresi, R.G. Henry, J. Oh, N. Papinutto, D. Pelletier, W. Rooney, W. Stern, N.L. Sicotte, D.S. Reich, R.T. Shinohara. An automated statistical technique for counting distinct multiple sclerosis lesions. American Journal of Neuroradiology, 2018; 39, 626-633.
lesion_centers <- function(prob_map, bin_map,
                           minCenterSize = 10, radius = 1) {
  prob_map_oro <- ants2oro(prob_map)
  scale <- ceiling((1 / mean(prob_map_oro@pixdim[2:4]))^3)
  phes <- hessian3D(
    image = prob_map_oro, mask = ants2oro(bin_map),
    radius
  )
  phes1 <- phes$eigval1
  phes2 <- phes$eigval2
  phes3 <- phes$eigval3

  cluster_map <- labelClusters(bin_map, minClusterSize = 20 * scale)

  les <- antsImageClone(cluster_map)
  les[les != 0] <- 1
  les[phes1 > 0 | phes2 > 0 | phes3 > 0] <- 0
  les <- labelClusters(les, minClusterSize = minCenterSize * scale)

  return(list(lesioncenters = les, lesioncount = max(les)))
}
