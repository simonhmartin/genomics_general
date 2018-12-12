
get_block_indices <- function(block_size, positions, chromosomes = NULL){
    if (is.null(chromosomes) == TRUE) {
        block_starts <- seq(min(positions), max(positions), block_size)
        
        block_ends <- block_starts + block_size - 1
        
        lapply(1:length(block_starts), function(x) which(positions >= block_starts[x] &
                                                         positions <= block_ends[x]))
        }
    else {
        chrom_names <- unique(chromosomes)
        
        block_starts <- lapply(chrom_names, function(chrom_name) seq(min(positions[chromosomes==chrom_name]),
                                                                     max(positions[chromosomes==chrom_name]), block_size)) 
        
        block_chroms <- unlist(lapply(1:length(block_starts), function(x) rep(chrom_names[x], length(block_starts[[x]]))))
        
        block_starts <- unlist(block_starts)
        
        block_ends <- block_starts + block_size - 1
        
        lapply(1:length(block_starts), function(x) which(chromosomes == block_chroms[x] &
                                                     positions >= block_starts[x] &
                                                     positions <= block_ends[x]))
        }
    }

#this function runs the jackknife procedure by calculating pseudovalues by removing one block at a time
#if the arguments specified by "..." are vectors, they will be indexed as they are.
#if they have two dimensions, they will be indexed along the first dimension
get_jackknife_sd <- function(block_indices, FUN, ...){
    n_blocks <- length(block_indices)
    args = list(...)
    overall_mean <- FUN(...)
    if (is.null(dim(args[1])) == TRUE){
        return(sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]]]))*(n_blocks-1))))
        }
    else{
        return(sd(sapply(1:n_blocks, function(i) overall_mean*n_blocks - do.call(FUN, lapply(args, function(a) a[-block_indices[[i]],]))*(n_blocks-1))))
        }
    }




