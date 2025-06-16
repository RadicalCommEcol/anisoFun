#' Calculate the distance of each species to its extinction border
#'
#' @description 
#'
#' @param A_int the interaction matrix. NB COMPETITIVE INTERACTIONS ARE NEGATIVE
#' @param r0 the intrinsic growth rate
#' @param norm the norm to embed feasibility domain, either "l1" or "l2". WE USE "l2" NORM.
#' @param nsample the number of sampling points on the border of feasibility domain. Min 500.
#' @param just_min whether to return all resistance values from the samples or taking the minimal value as the estimated resistance
#'   all if TRUE, returns the minimum distance to the border of feasibility of the r0 vector.
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map_dbl
#' @importFrom utils combn
#'
#' @return an df with two numeric arrays = the distance of species i to its extinctin border. A binary array (extinct) indicating if the species is feasible (0) or not (1)
#' @export
Measure_Distances <- function(A_int, r0, norm = "l2", nsample = 1000, just_min = FALSE) {
    norm1 <- function(a) {
        if (is.matrix(a)) {
            a <- apply(a, 2, function(x) x / sum(x))
        } else {
            a <- a / sum(a)
        }
        return(a)
    }
    norm2 <- function(a) {
        if (is.matrix(a)) {
            a <- apply(a, 2, function(x) x / sqrt(sum(x^2)))
        } else {
            a <- a / sqrt(sum(a^2))
        }
        return(a)
    }

    feasibility_resistance_full <- function(matA, r, norm = "l2", nsample = 100, all = FALSE) {
        simplex_sampling <- function(m, n) {
            r <- list()
            for (j in 1:m) {
                dist <- c(sort(runif(n - 1, 0, 1)), 1)
                r[[j]] <- c(dist[1], diff(dist))
            }
            return(r)
        }
        euclidean_distance <- function(a, b) {
            sqrt(sum((a - b)^2))
        }
        arc_length <- function(a, b) {
            acos(sum(a * b))
        }
        matA <- -matA # feasibility_resistance_full originally was designed for competitive interactions as positive

        vertices <- combn(seq_len(nrow(matA)), nrow(matA) - 1, simplify = FALSE)
        distances_list <- list()
        for (i in seq_len(length(vertices))) {
        vertex <- vertices[[i]]
        if (norm == "l1") {
            r <- norm1(r); matA <- norm1(matA)
            distances_list[[i]] <- 1:nsample %>%
            map_dbl(function(x) {
                t <- unlist(simplex_sampling(1, nrow(matA) - 1))
                border_point <- c(matA[, vertex] %*% t)
                euclidean_distance(r, border_point)
            })
        }
        if (norm == "l2") {
            r <- norm2(r); matA <- norm2(matA)
            distances_list[[i]] <- 1:nsample %>%
            map_dbl(function(x) {
                t <- unlist(simplex_sampling(1, nrow(matA) - 1))
                border_point <- c(matA[, vertex] %*% t)
                border_point_norm <- border_point / sqrt(sum(border_point^2))
                arc_length(r, border_point_norm)
            })
        }
        }
        names(distances_list) <- unlist(lapply(vertices, paste, collapse = "_"))
        if (all) {
        return(distances_list)
        } else {
        return(sapply(distances_list, min))
        }
        }
  
    # 1.Checks 
    # if the input matrix is square
    if (nrow(A_int) != ncol(A_int)) {
        stop("The interaction matrix must be square.")
    }   
    # if the length of r0 matches the number of rows/columns in A_int
    if (length(r0) != nrow(A_int)) {
        stop("The length of r0 must match the number of rows/columns in A_int.")
    }
 
    nspecies <- nrow(A_int)

    colnames(A_int) <- paste0("sp_",1:nrow(A_int)) 
    rownames(A_int) <- paste0("sp_",1:nrow(A_int))
    
    # 2. Find the equilibrium abundance of the species. If the equilibrium abundance is negative, the species is not feasible.
    aux <- solve(-A_int,r0)

    # 3. Calculate the distances to the border of feasibility for each species
    distances <- feasibility_resistance_full(A_int, r0, nsample = nsample, all = FALSE)
    if(all(aux>0)){     
        extinct <- rep(0,nspecies) 
    }else {        
        extinct <- rep(0,nspecies)
        extinct[aux<=0] <- rep(1,sum(aux<=0)) # extinct = 1 if the species is not feasible
    }

    if(just_min) {
        names(extinct) <- paste0("sp",1:nspecies)
        print(distances)
        df <- list(
            distance = min(distances),
            extinct = unlist(extinct)
        )
    } else {
        names(distances) <- paste0("sp",1:nspecies)
        names(extinct) <- paste0("sp",1:nspecies)
        df <- data.frame(
            distance = unlist(distances),
            extinct = unlist(extinct)
        )
    }

    return(df)

}
