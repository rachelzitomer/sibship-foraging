distance_to_focal_deprecated <- function(conductance, s, cells_per_block=5000, average_conductance=TRUE){

  stopifnot(all(conductance > 0))

  symm <- function(X) (X + t(X))/2

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N     <- length(conductance)
  Q     <- s$laplacian
  if (average_conductance){
    Q@x[] <- -conductance[s$adj[1,]+1] - conductance[s$adj[2,]+1]
  } else {
    Q@x[] <- -1/(1/conductance[s$adj[1,]+1] + 1/conductance[s$adj[2,]+1])
  }

  # Eq. ??? in radish paper
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- Matrix::forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
  ones <- matrix(1, N, 1)
  v    <- sqrt(ones / N)
  Z    <- Matrix::Diagonal(N)[,s$demes]
  Zn   <- In %*% Z - (In %*% v) %*% (t(v) %*% Z)
  LQn  <- Matrix::update(s$choleski, Qn)
  G    <- Matrix::solve(LQn, Zn)
  tG   <- t(as.matrix(G))
  E    <- Matrix::t(Zn) %*% G 
  E2   <- Matrix::t(In) %*% G - (v %*% (t(v) %*% (Matrix::t(In) %*% G)))

  lE <- rep(0, N)
  nonfocal <- 1:(N-1)
  num_blocks <- ceiling((N-1)/cells_per_block)
  blocks <- as.numeric(cut(nonfocal, breaks=num_blocks))
  for(i in unique(blocks))
  {
    idx <- which(blocks == i)
    lZ  <- Matrix::Diagonal(N)[,idx,drop=FALSE]
    lZn <- In %*% lZ - (In %*% v) %*% (t(v) %*% lZ)
    lG  <- Matrix::solve(LQn, lZn)
    lE2 <- Matrix::t(In) %*% lG - (v %*% (t(v) %*% (Matrix::t(In) %*% lG)))
    llE    <- Matrix::t(lZn) %*% lG 
    lE[idx] <- diag(as.matrix(llE))
    lE[N] <- lE[N] - sum(lE2[N,])
  }

  diag(as.matrix(E)) + t(-2*as.matrix(E2) + lE)
}

distance_to_focal_raw <- function(conductance, s, cells_per_block=5000, average_conductance=TRUE){

  stopifnot(all(conductance > 0))

  symm <- function(X) (X + t(X))/2

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N     <- length(conductance)
  Q     <- s$laplacian
  if (average_conductance){
    Q@x[] <- -conductance[s$adj[1,]+1] - conductance[s$adj[2,]+1]
  } else {
    Q@x[] <- -1/(1/conductance[s$adj[1,]+1] + 1/conductance[s$adj[2,]+1])
  }

  # Eq. ??? in radish paper
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- Matrix::forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
  ones <- matrix(1, N, 1)
  v    <- sqrt(ones / N)
  Z    <- Matrix::Diagonal(N)[,s$demes]
  #speed this up using the same strategy
  Zn   <- In %*% Z - (In %*% v) %*% (t(v) %*% Z)
  LQn  <- Matrix::update(s$choleski, Qn)
  G    <- Matrix::solve(LQn, Zn)
  tG   <- t(as.matrix(G))
  E    <- Matrix::t(Zn) %*% G 
  E2   <- Matrix::t(In) %*% G - (v %*% (t(v) %*% (Matrix::t(In) %*% G)))

  lE <- rep(0, N)
  nonfocal <- 1:(N-1)
  num_blocks <- ceiling((N-1)/cells_per_block)
  if (num_blocks > 1)
  {
    blocks <- as.numeric(cut(nonfocal, breaks=num_blocks))
  } else {
    blocks <- rep(1, length(nonfocal))
  }

  B <- In %*% (v^2)
  Binv <- Matrix::solve(LQn, B)
  subproblem <- function(idx){
    lZ  <- Matrix::Diagonal(N)[,idx,drop=FALSE]
    A <- In %*% lZ
    Ainv <- Matrix::solve(LQn, A)
    a1 <- Matrix::diag(Matrix::t(In %*% lZ) %*% Ainv) - 1/N * Matrix::colSums(Ainv)
    #a1_orig <- Matrix::diag(Matrix::t(lZn) %*% Ainv)
    b1 <- Matrix::t(In %*% lZ) %*% Binv - 1/N * sum(Binv)
    #b1_orig <- Matrix::t(lZn) %*% Binv
    #max(abs( (a - b) - Matrix::diag(llE) ))

    a2 <- sum(-1/N * Matrix::colSums(Ainv))
    #a2_orig <- Matrix::t(In) %*% Ainv - (v %*% (t(v) %*% (Matrix::t(In) %*% Ainv)))
    #sum(a2_orig[N,])

    b2 <- (cells_per_block - ncol(Zn)) * sum(1/N * Binv)
    #Binv2 <- outer(Binv, rep(1, block_size - ncol(Zn)) )
    #b2_orig <- Matrix::t(In) %*% Binv2 - (v %*% (t(v) %*% (Matrix::t(In) %*% Binv2)))
    #sum(b2_orig[N,])

    #sum(a2_orig[N,]) - b2_orig[N]
    #a2 - b2

    list(llE = a1 - b1, lE2 = a2 - b2)
  }

  for(i in unique(blocks))
  {
    idx <- which(blocks == i)
#    lZ  <- Matrix::Diagonal(N)[,idx,drop=FALSE]
#    lZn <- In %*% lZ - (In %*% v) %*% (t(v) %*% lZ)
#    # lZn <- In %*% lZ - 1/N
#    # this is just In %*% lZ - ones / N because v is constant and lZ are indicators
#    lG  <- Matrix::solve(LQn, lZn)
#    # because of the above, we have Matrix::solve(LQn, In %*% lZ) - Matrix::solve(LQn, matrix(ones / N, nrow(In), ncol(lZ)))
#    # let these be A, B
#    lE2 <- Matrix::t(In) %*% lG - (v %*% (t(v) %*% (Matrix::t(In) %*% lG)))
#    # t(In) %*% A - t(In) %*% B - 1/N * t(In) %*% A + 1/N * t(In) %*% B
#    # and t(In) %*% B is rbind(matrix(1/N, ...), 0) 
#    llE    <- Matrix::t(lZn) %*% lG 
#    # t(In %*% lZ) %*% A - t(In %*% lZ) %*% B - 1/N * t(In %*% lZ) %*% A + 1/N * t(In %*% lZ) %*% B
#    # so need diag(t(In %*% lZ) %*% A) - diag(t(In %*% lZ) %*% B)
#    lE[idx] <- diag(as.matrix(llE))
#    lE[N] <- lE[N] - sum(lE2[N,])
    tmp <- subproblem(idx)
    lE[idx] <- tmp$llE
    lE[N] <- lE[N] - tmp$lE2

   # #this works. taken this, disassemble
   # Binv2 <- outer(Binv, rep(1, 975))
   # boof <- Matrix::t(In) %*% (Ainv - Binv2) - (v %*% (t(v) %*% (Matrix::t(In) %*% (Ainv - Binv2) )))
   # sum(boof[N,])
   # sum(lE2[N,])
   # #----

  }

  diag(as.matrix(E)) + t(-2*as.matrix(E2) + lE)
}

#TODO: clean this up, these are kinda redundant -- could just delete as we are now doing resistance model manually for flexibility
distance_to_focal <- function(theta, f, s, num_blocks=2, conductance=TRUE){

  symm <- function(X) (X + t(X))/2

  # conductance
  C <- f(theta)

  # Form the Laplacian. "adj" is assumed to contain at a minimum
  # the upper triangular part of the Laplacian (e.g. all edges [i,j]
  # where i < j). Duplicated edges are ignored.
  N     <- length(C$conductance)
  Q     <- s$laplacian
  if (conductance){
    Q@x[] <- -C$conductance[s$adj[1,]+1] - C$conductance[s$adj[2,]+1]
  } else {
    Q@x[] <- -1/(1/C$conductance[s$adj[1,]+1] + 1/C$conductance[s$adj[2,]+1])
  }

  # Eq. ??? in radish paper
  Qd   <- Matrix::Diagonal(N, x = -Matrix::rowSums(Q))
  In   <- Matrix::Diagonal(N)[-N,]
  Qn   <- Matrix::forceSymmetric(In %*% (Q + Qd) %*% Matrix::t(In))
  ones <- matrix(1, N, 1)
  v    <- sqrt(ones / N)
  Z    <- Matrix::Diagonal(N)[,s$demes]
  Zn   <- In %*% Z - (In %*% v) %*% (t(v) %*% Z)
  LQn  <- Matrix::update(s$choleski, Qn)
  G    <- Matrix::solve(LQn, Zn)
  tG   <- t(as.matrix(G))
  E    <- Matrix::t(Zn) %*% G 
  E2   <- t(In) %*% G - v %*% t(v) %*% t(In) %*% G #TODO replace E with this

  lE <- rep(0, N)
  nonfocal <- 1:(N-1)
  blocks <- as.numeric(cut(nonfocal, breaks=num_blocks))
  for(i in unique(blocks))
  {
    idx <- which(blocks == i)
    lZ  <- Matrix::Diagonal(N)[,idx,drop=FALSE]
    lZn <- In %*% lZ - (In %*% v) %*% (t(v) %*% lZ)
    lG  <- Matrix::solve(LQn, lZn)
    lE2 <- t(In) %*% lG - v %*% t(v) %*% t(In) %*% lG
    llE    <- Matrix::t(lZn) %*% lG 
    lE[idx] <- diag(llE)
    lE[N] <- lE[N] - sum(lE2[N,])
  }

  diag(E) + t(-2*as.matrix(E2) + lE)
}

radish_distance_to_focal <- function(theta,
                                     formula, 
                                     data,
                                     conductance_model = radish::loglinear_conductance, 
                                     conductance = TRUE)
{
  stopifnot(is.matrix(theta))
  # TODO: check that colnames of theta match formula

  # get response, remove lhs from formula
  terms    <- terms(formula)
  is_ibd   <- length(attr(terms, "factors")) == 0

  stopifnot(!is_ibd) #IBD; nothing to do

  formula  <- reformulate(attr(terms, "term.labels"))

  # "conductance_model" (a factory) is then responsible for parsing formula,
  # constructing design matrix, and returning actual "conductance_model"
  conductance_model <- conductance_model(formula, data$x) #TODO: is_ibd here and have factory modify accordingly
  default <- attr(conductance_model, "default")

  stopifnot(ncol(theta) == length(default))

  colnames(theta) <- names(default)
  
  output <- array(NA, c(length(surface$demes), nrow(surface$laplacian), nrow(theta)))
  for (i in 1:nrow(theta)){
    output[,,i] <- 
      distance_to_focal(theta = theta[i,], 
        f = conductance_model, 
        s = data, 
        conductance = conductance) 
  }
  list(theta = data.frame(theta), distance = output)
}

