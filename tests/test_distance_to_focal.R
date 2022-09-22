#-------------- test correctness of resistance distance calculation ----------------#

library(sibships)
library(gdistance)
library(raster)

is_equal <- function(x, y)
{
  abs(x - y) < sqrt(.Machine$double.eps)
}

toy_problem <- function(dim=20, n=5, seed=1)
{
  set.seed(seed)
  data(volcano)
  data <- volcano[1:min(dim, nrow(volcano)), 1:min(dim, ncol(volcano))]/100
  landscape <- raster(data)
  coords <- cbind(runif(n), runif(n))
  cells <- raster::cellFromXY(landscape, coords)
  while(length(unique(cells)) < n & max(cells) < ncell(landscape))
  {
    coords <- cbind(runif(n), runif(n))
    cells <- raster::cellFromXY(landscape, coords)
  }
  coords <- SpatialPoints(coords)
  list(landscape=landscape, coords=coords, cells=cells)
}

# ----------- test data ---------- #
test_data <- toy_problem()
landscape <- test_data$landscape
cells <- test_data$cells
coords <- test_data$coords

# --------------- gdistance --------------- #
# sort of annoying, can't get distances to last cell in raster this way
transition_operator <- gdistance::transition(landscape, function(x) 1/mean(x), 8)
all_coords <- raster::xyFromCell(landscape, 2:ncell(landscape)-1)
gdistance_distance <- as.matrix(gdistance::commuteDistance(transition_operator, all_coords))
gdistance_result <- gdistance_distance[cells,]

# ----------- radish + this package ------- #
surface <- radish::conductance_surface(
  raster::stack(list(landscape=landscape)),
  coords,
  directions=8
)
radish_distance <- sibships::distance_to_focal_raw(
  conductance=1./landscape[],
  s=surface, 
  cells_per_block=13,
  average_conductance=FALSE
)
radish_result <- radish_distance[,-c(400)]

# ------------ are they the same (up to constant factor)? --------- #
stopifnot(is_equal(1, cor(c(gdistance_result), c(radish_result))))

# ------------ check the last fucking cell that gdistance leaves out -------- #
# using the inefficient but easy calculation for resistance distance via eigendecomposition
Q <- surface$laplacian
Q@x[] <- -1/(landscape[surface$adj[1,]+1] + landscape[surface$adj[2,]+1])
Qd <- Matrix::Diagonal(nrow(Q), x = -Matrix::rowSums(Q))
Q <- Q + Qd
Qinv <- MASS::ginv(as.matrix(Q)) #the inefficient but easy way to calculate resistance distance
Rd <- outer(rep(1, nrow(Qinv)), diag(Qinv)) + outer(diag(Qinv), rep(1, nrow(Qinv))) - 2*Qinv
naive_distance <- Rd[cells,]
stopifnot(is_equal(1, cor(c(naive_distance), c(radish_distance))))

# ------------ check that unblocked calculation works -----------#
radish_distance_2 <- sibships::distance_to_focal_raw(
  conductance=1./landscape[],
  s=surface, 
  average_conductance=FALSE
)
stopifnot(is_equal(1, cor(c(radish_distance_2), c(radish_distance))))
