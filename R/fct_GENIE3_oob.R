#' @title GENIE3OOB
#'
#' @description \code{GENIE3OOB} Infers a gene regulatory network (in the form of a pvalues-filled adjacency matrix)
#' from expression data, using ensembles of regression trees.
#'
#' @param exprMatrix Expression matrix (genes x samples). Every row is a gene, every column is a sample.
#' @param regulators Subset of genes used as candidate regulators. Must be either a
#' vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g.
#' \code{c("at_12377", "at_10912")}. The default value NULL means that all the genes are used as candidate regulators.
#' @param targets Subset of genes to which potential regulators will be calculated. Must be either
#' a vector of indices, e.g. \code{c(1,5,6,7)}, or a vector of gene names, e.g. \code{c("at_12377", "at_10912")}.
#' If NULL (default), regulators will be calculated for all genes in the input matrix.
#' @param nTrees Number of trees in an ensemble for each target gene. Default: 1000.
#' @param nCores Number of cores to use for parallel computing. Default: 1.
#' @param verbose If set to TRUE, a feedback on the progress of the calculations is given. Default: TRUE
#'
#' @return matrix containing oob importances for all TF-gene pairs
GENIE3OOB <-
  function(exprMatrix,
           regulators = NULL,
           targets = NULL,
           nTrees = 1000,
           nCores = 1,
           verbose = FALSE)
  {
    .checkArguments(
      exprMatrix = exprMatrix,
      regulators = regulators,
      targets = targets,
      nTrees = nTrees,
      nCores = nCores,
      verbose = verbose
    )
    
    if (is.numeric(regulators))
      regulators <- rownames(exprMatrix)[regulators]
    
    ############################################################
    # transpose expression matrix to (samples x genes)
    exprMatrixT <- t(exprMatrix)
    rm(exprMatrix)
    num.samples <- nrow(exprMatrixT)
    allGeneNames <- colnames(exprMatrixT)
    
    # get names of input genes
    if (is.null(regulators))
    {
      regulatorNames <- allGeneNames
    } else
    {
      # input gene indices given as integers
      if (is.numeric(regulators))
      {
        regulatorNames <- allGeneNames[regulators]
        # input gene indices given as names
      } else
      {
        regulatorNames <- regulators
        # for security, abort if some input gene name is not in gene names
        missingGeneNames <- setdiff(regulatorNames, allGeneNames)
        if (length(missingGeneNames) != 0)
          stop(paste(
            "Regulator genes missing from the expression matrix:",
            paste(missingGeneNames, collapse =
                    ", ")
          ))
      }
    }
    regulatorNames <- sort(regulatorNames)
    rm(regulators)
    
    # get names of target genes
    if (is.null(targets))
    {
      targetNames <- allGeneNames
    } else
    {
      # input gene indices given as integers
      if (is.numeric(targets))
      {
        targetNames <- allGeneNames[targets]
        # input gene indices given as names
      } else
      {
        targetNames <- targets
        # for security, abort if some input gene name is not in gene names
        missingGeneNames <- setdiff(targetNames, allGeneNames)
        if (length(missingGeneNames) != 0)
          stop(paste(
            "Target genes missing from the expression matrix:",
            paste(missingGeneNames, collapse =
                    ", ")
          ))
      }
    }
    targetNames <- sort(targetNames)
    rm(targets)
    
    if (verbose)
      message(
        paste(
          "Starting GENIE3 network inference with MSEincrease on OOB as importance metric.",
          "\nK: ",
          "sqrt",
          "\nNumber of trees: ",
          nTrees,
          sep = ""
        )
      )
    
    # setup weight matrix
    Matrix <-
      matrix(0,
             nrow = length(regulatorNames),
             ncol = length(targetNames))
    rownames(Matrix) <- regulatorNames
    colnames(Matrix) <- targetNames
    
    #if (!foreach::getDoParRegistered()) {
    doParallel::registerDoParallel(cores  = nCores)
    #}
    
    
    if (verbose)
      message(paste("\nUsing", foreach::getDoParWorkers(), "cores."))
    "%dopar%" <- foreach::"%dopar%"
    suppressPackageStartupMessages(Matrix.reg <-
                                     doRNG::"%dorng%"(
                                       foreach::foreach(targetName = targetNames, .combine = cbind),
                                       {
                                         # remove target gene from input genes
                                         theseRegulatorNames <-
                                           setdiff(regulatorNames, targetName)
                                         numRegulators <-
                                           length(theseRegulatorNames)
                                         mtry <-
                                           round(sqrt(numRegulators))
                                         
                                         x <-
                                           exprMatrixT[, theseRegulatorNames]
                                         y <-
                                           exprMatrixT[, targetName]
                                         
                                         rf <-
                                           randomForest::randomForest(
                                             x,
                                             y,
                                             mtry = mtry,
                                             ntree = nTrees,
                                             importance = TRUE,
                                             nodesize = 1
                                           )
                                         
                                         im <-
                                           randomForest::importance(rf)[, "%IncMSE"]
                                         
                                         c(setNames(0, targetName), setNames(im, names(im)))[regulatorNames]
                                       }
                                     ))
    
    attr(Matrix.reg, "rng") <- NULL
    Matrix[regulatorNames, ] <- Matrix.reg[regulatorNames, ]
    
    if (verbose & sum(is.na(Matrix)) > 0) {
      warning(
        paste(
          sum(is.na(Matrix)),
          "Na importance values were returned,
      it may be caused by too few samples. You can run
      network inference with the node purity importance metric
      instance of OOB MSE increase, to avoid this issue."
        )
      )
    }
    
    
    return(Matrix)
  }


.checkArguments <-
  function(exprMatrix,
           regulators,
           targets,
           nTrees,
           nCores,
           verbose)
  {
    ############################################################
    # check input arguments
    if (!is.matrix(exprMatrix) && !is.array(exprMatrix)) {
      stop(
        "Parameter exprMatrix must be a two-dimensional matrix where each row corresponds to a
        gene and each column corresponds to a condition/sample/cell."
      )
    }
    
    if (length(dim(exprMatrix)) != 2) {
      stop(
        "Parameter exprMatrix must be a two-dimensional matrix where each row corresponds to a
        gene and each column corresponds to a condition/sample/cell."
      )
    }
    
    if (is.null(rownames(exprMatrix))) {
      stop("exprMatrix must contain the names of the genes as rownames.")
    }
    
    countGeneNames <- table(rownames(exprMatrix))
    nonUniqueGeneNames <- countGeneNames[countGeneNames > 1]
    if (length(nonUniqueGeneNames) > 0)
      stop("The following gene IDs (rownames) are not unique: ",
           paste(names(nonUniqueGeneNames), collapse = ", "))
    
    if (!is.numeric(nTrees) || nTrees < 1) {
      stop("Parameter nTrees should be a stricly positive integer.")
    }
    
    if (!is.null(regulators))
    {
      if (length(regulators) < 2)
        stop("Provide at least 2 potential regulators.")
      
      if (!is.vector(regulators)) {
        stop("Parameter 'regulators' must a vector (of indices or gene names).")
      }
      
      if (is.numeric(regulators)) {
        if (max(regulators) > nrow(exprMatrix))
          stop("At least one index in 'regulators' exceeds the number of genes.")
        if (min(regulators) <= 0)
          stop("The indexes in 'regulators' should be >=1.")
      }
      
      if (any(table(regulators) > 1))
        stop("Please, provide each regulator (name/ID) only once.")
      
      if (is.character(regulators)) {
        regulatorsInMatrix <- intersect(regulators, rownames(exprMatrix))
        if (length(regulatorsInMatrix) == 0)
          stop("The genes must contain at least one regulators")
        
        if (length(regulatorsInMatrix) < length(regulators))
          warning(
            "Only",
            length(regulatorsInMatrix),
            "out of",
            length(regulators),
            " candidate regulators (IDs/names) are in the expression matrix."
          )
      }
    }
    
    if (!is.null(targets))
    {
      if (!is.vector(targets)) {
        stop("Parameter 'targets' must a vector (of indices or gene names).")
      }
      
      if (is.numeric(targets)) {
        if (max(targets) > nrow(exprMatrix))
          stop("At least one index in 'targets' exceeds the number of genes.")
        if (min(targets) <= 0)
          stop("The indexes in 'targets' should be >=1.")
      }
      
      if (any(table(targets) > 1))
        stop("Please, provide each target (name/ID) only once.")
      
      if (is.character(targets)) {
        targetsInMatrix <- intersect(targets, rownames(exprMatrix))
        if (length(targetsInMatrix) == 0)
          stop("The genes must contain at least one targets.")
        
        
        if (length(targetsInMatrix) < length(targets))
          warning(
            "Only",
            length(targetsInMatrix),
            "out of",
            length(targets),
            "target genes (IDs/names) are in the expression matrix."
          )
      }
    }
    
    if (!is.numeric(nCores) || nCores < 1)
    {
      stop("Parameter nCores should be a stricly positive integer.")
    }
  }