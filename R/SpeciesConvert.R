
#' Human-mouse ortholog pairs from Ensembl 87.
#'
#' @export
get_ortholog_table = function(){
  data( orthologs_ens87 )
  return( orthologs_ens87 )
}

# # Set up data on human-mouse orthologs as named vectors for fast access
ortholog_table = get_ortholog_table()
human_dupes = duplicated( ortholog_table$humansym )
mouse_dupes = duplicated( ortholog_table$mousesym )
h2m = setNames( ortholog_table$mousesym, nm = ortholog_table$humansym )[!human_dupes]
m2h = setNames( ortholog_table$humansym, nm = ortholog_table$mousesym )[!mouse_dupes]

#' Return the ortholog of a given gene or NA (if no match). Not vectorized, so use get_ortholog instead.
#'
get_ortholog_nonvec = function( gene, from, to ){
  get_if_present = function( gene, db ){
    if( gene %in% names(db) ){
      return( db[gene] )
    } else {
      return( NA )
    }
  }

  if       ( from == "human" && to == "mouse"){
    return( get_if_present( gene, h2m ) )
  } else if( from == "mouse" && to == "human"){
    return( get_if_present( gene, m2h ) )
  } else {
    warning(' The only working options are from = "human", to = "mouse" and from = "mouse", to = "human". Returning your gene unaltered. ')
    return( gene )
  }
}

#' Return the ortholog of a given gene or NA (if no match).
#' Human and mouse only so far.
#'
#' @export
get_ortholog = function(x, from, to, ... ) {
  x %<>% as.character
  y = sapply(x, get_ortholog_nonvec, from = from, to = to, ...)
  names(y) = x
  return(y)
}

#' Same as get_ortholog but returns just T or F.
#'
#' @export
has_ortholog = function( ... ){ !is.na( get_ortholog( ... ) ) }




#' Convert a gene expression matrix from one species to another.
#'
#' @details If two genes have the same ortholog, the molecule counts get added.
#' If a gene has no ortholog, it is omitted.
#' The input must be a matrix with genes stored in rownames( raw_dge ).
#'
#' @export
#'
convert_species_rownames = function( expr, from, to, ... ){
  cat( paste( "\nConverting from", from, "to", to, "...\n" ) )
  genes = rownames(expr)
  eligible_genes = genes[ has_ortholog( genes, from=from, to=to ) ]
  if(length(eligible_genes) == 0 ){
    stop("Can't convert this object: no genes have orthologs!\n")
  }
  orthologs_by_gene                = eligible_genes %>% get_ortholog( ., from, to )
  to_merge = intersect(rownames(expr),eligible_genes)
  expr %<>% extract(to_merge, )
  return( merge_dupe_rows( X = expr, row.names = orthologs_by_gene[to_merge], ... ) )
}

#' Merge duplicate rows of a Matrix (from the Matrix package).
#'
#' @param X Matrix to start with.
#' @param FUN Merging function.
#'
#' @export
#'
merge_dupe_rows = function( X, row.names, FUN = Matrix::colSums, verbose = F){
  X %<>% (Matrix::Matrix)
  assertthat::are_equal(length(row.names), nrow(X))
  rownames(X) = row.names
  idx_dup = row.names %>% duplicated
  dupes = row.names %>% extract(idx_dup) %>% unique
  not_dupes = row.names %>% setdiff(dupes)
  if( length(dupes) == 0) {return(X)}
  Y = X[not_dupes, ]; rownames(Y) = not_dupes
  Z = Matrix::Matrix(0, ncol = ncol(X), nrow = length(dupes), sparse = T)
  rownames(Z) = dupes
  colnames(Z) = colnames(X)
  if( verbose ){
    cat("There are ", length(dupes), " duplicate rows.\n" )
    cat("Original dimensions: ",               dim(X)[1], " by ", dim(X)[2], "\n" )
    cat("Dimension with singlets only: " ,     dim(Y)[1], " by ", dim(Y)[2], "\n" )
    cat("Dimension of merged row template: " , dim(Z)[1], " by ", dim(Z)[2], "\n" )
    cat("Merging...\n" )
  }

  for (d in dupes) {
    rows_use = X[row.names == d, , drop = F]
    if (verbose) {
      cat(d, ": ", dim(rows_use)[1], "occurrences \n")
    }
    Z[d, ] = FUN(rows_use)
  }

  W = rbind(Y, Z)
  if(verbose){ cat("Dimension of final result: " , dim(W)[1], " by ", dim(W)[2], "\n" ) }
  return(W)
}
