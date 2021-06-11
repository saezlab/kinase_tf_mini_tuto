df_to_viper_regulon <- function(df)
{
  names(df) <- c("feature","pathway","sign")
  df <- df[complete.cases(df),]

  pathway_regulon <- list(0)
  i <- 1
  for(pathway in unique(df$pathway))
  {
    pathway_feature_list <- list(0)
    features <- df[df$pathway == pathway, 3]
    names(features) <- df[df$pathway == pathway, 1]
    pathway_feature_list[[1]] <- features
    pathway_feature_list[[2]] <- rep(1,length(features))
    names(pathway_feature_list) <- c("tfmode","likelihood")

    pathway_regulon[[i]] <- pathway_feature_list
    i <- i+1
  }
  names(pathway_regulon) <- unique(df$pathway)
  return(pathway_regulon)
}

plot_top_features <- function(data, n_top) {

  colnames(data) <- "value"

  arranged <- data %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "id") %>%
    arrange(desc(value))

  top_up <- slice_head(arranged, n = n_top)
  top_down <- slice_tail(arranged, n = n_top)

  p <- bind_rows(list(up = top_up, down = top_down), .id = "status") %>%
    mutate(id = fct_inorder(id)) %>%
    ggplot(aes(x = value, y = id, fill = status)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    theme_bw()

  return(p)

}
#' Translate matrix rownames
#'
#' Translates the rownames of the input matrix into the
#' desired ids using a translator data frame. When
#' input ids maps to several target ids, uses the summarise
#' function to resolve conflicts.
#'
#' @param mat Input matrix. Should have rownames.
#' @param df A 2-columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#'
#' @export
#'
#' @importFrom dplyr %>% group_by summarise_all sym select
#' @importFrom tibble rownames_to_column column_to_rownames
#'
translateMatrix <- function(mat, df, sourceKey, targetKey, summariseFun) {
  
  # if not function is passed, use first match as default function
  if(missing(summariseFun)) {
    message("------------------------------------------------")
    message("No input summarise function detected, using first match on multi-mapping situations.")
    summariseFun <- function(x) return(x[1])
  }
  # check matrix and translator df, printing mapping info
  inputMatrixCheck <- validateMatrix(mat)
  df <- validateTranslatorDf(df, sourceKey, targetKey)
  # transform mat into data frame and rownames into column
  xDf <- as.data.frame(mat) %>%
    tibble::rownames_to_column("mergingVariable")
  # merge with data matrix
  mergedDf <- merge(x = df, y = xDf, by.x = sourceKey, by.y = "mergingVariable")
  # get unique sources and keys
  uniqueSources <- df[, sourceKey][df[, sourceKey] %in% names(table(df[, sourceKey]))[table(df[, sourceKey]) == 1]]
  uniqueTargets <- df[, targetKey][df[, targetKey] %in% names(table(df[, targetKey]))[table(df[, targetKey]) == 1]]
  uniqueMatches <- mergedDf[, sourceKey] %in% uniqueSources & mergedDf[,targetKey] %in% uniqueTargets
  # remove source id, group and split matrixes
  if(sum(uniqueMatches) != 0) {
    mergedDf[, sourceKey] <- NULL
    uniqueMat <- mergedDf[uniqueMatches, ] %>%
      tibble::rownames_to_column(var = "toDiscard") %>%
      dplyr::select(-toDiscard) %>%
      column_to_rownames(var = targetKey) %>%
      as.matrix()
  } else {
    uniqueMat <- NULL
  }
  # tidy eval for targetKey (https://tidyeval.tidyverse.org/introduction.html)
  if(sum(!uniqueMatches) != 0) {
    duplicatedDf <- mergedDf[ !uniqueMatches, ]
    ids <- unique(duplicatedDf[, targetKey])
    duplicatedMat <- sapply(ids, function(id) {
      outVec <- subset(duplicatedDf, duplicatedDf[, targetKey] == id) %>%
        dplyr::select( -any_of(targetKey)) %>%
        as.matrix() %>%
        apply(., 2, FUN = summariseFun)
      return(outVec)
    }) %>%
      # transpose
      t()
    # and set new rownames
    rownames(duplicatedMat) <- ids
  } else {
    duplicatedMat <- NULL
  }
  # bind cols and return matrix
  outMat <- rbind(uniqueMat, duplicatedMat)
  return(outMat)
  
}

#' Translate matrix rownames using annotation package
#'
#' Uses an annotation package like
#' \href{https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html}{org.Hs.eg.db} or
#' \href{https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html}{org.Mm.eg.db}
#' to translate the rownames of the input matrix using \link[biokit]{translateMatrix}.
#'
#' @param mat Input matrix.
#' @param db Annotation package object to use for the translation.
#' @param sourceKey Source key in the DB.
#' @param targetKey Target key in the DB.
#' @param summariseFun Function used to resolve multi-mapping situations.
#'
#' @return A matrix with translated and summarised ids.
#'
#' @importFrom AnnotationDbi select
#'
#' @export
#'
translateMatrixWithDb <- function(mat, db, sourceKey, targetKey, summariseFun) {
  
  # get mapping df
  srcTarget <- AnnotationDbi::select(db, keys = rownames(mat), keytype = sourceKey, columns = targetKey)
  # translate matrix
  translatedMat <- translateMatrix(mat = mat, df = srcTarget, sourceKey = sourceKey,
                                   targetKey = targetKey, summariseFun = summariseFun)
  return(translatedMat)
  
}

#' Validate matrix
#'
#' Evaluates integrity of the input matrix.
#'
#' @param mat The matrix to validate.
#' @param checkRowNames Check rownames are present? Stop execution if TRUE and mat does not contain rownames.
#' @param allowNa Allow NA presence in mat? Stop execution if NAs are found in mat.
#'
#'
#' @return TRUE if no errors are found.
#'
validateMatrix <- function(mat, checkRowNames = TRUE, allowNa = TRUE) {
  
  # check class
  matClass <- class(mat)
  if(! all(matClass == c("matrix", "array"))) stop("Input is not a matrix.")
  # check rownames
  if(checkRowNames & is.null(rownames(mat))) stop("Matrix does not have rownames.")
  # check NAs
  naCount <- sum(is.na(mat))
  if(naCount != 0 & allowNa) message("Matrix contain NAs. This may produce unexpected results.")
  if(naCount != 0 & (! allowNa)) stop("Matrix contain NAs.")
  return(TRUE)
  
}

#' Validate sample information
#'
#' Evaluates integrity of the grouping variable at the sample information data frame.
#'
#' @param sampInfo The sample information data frame.
#' @param groupCol The column containing the grouping variable.
#' @param checkNames Check if grouping variable contains syntactically valid names? If TRUE, they are fixed in the output data frame.
#' @param checkSingleSample Check if any group contains only one sample? A warning message appears if true.
#' @param mat Matrix containing data. If provided, it will check number of columns in mat.
#'
#' @return The checked sampInfo data frame
#'
validateSampInfo <- function(sampInfo, groupCol, mat,
                             checkNames = TRUE, checkSingleSample = TRUE) {
  
  # check not syntactically valid names
  if(checkNames) {
    groupVar <- sampInfo[ , groupCol]
    goodNames <- make.names(groupVar)
    if(any(goodNames != groupVar)) {
      message("Some levels at the grouping variable are not syntactically valid names. They will be fixed.")
      sampInfo[ , groupCol] <- goodNames
    }
  }
  # check number of rows is same than matrix ncol
  if(!missing(mat)) {
    if(ncol(mat) != nrow(sampInfo)) {
      stop("Number of rows for sampInfo is not equal to the number of matrix columns.")
    }
  }
  # check presence of single sample groups
  if(checkSingleSample) {
    groupN <- table(sampInfo[, groupCol])
    if(any(groupN == 1)) {
      message("Some sample groups contains only 1 sample. This may produce unexpected results.")
    }
  }
  return(sampInfo)
  
}

#' Evaluate and print translator df information
#'
#' @param df A 2 columns translator data frame with source to target ids.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#'
messageMappingInfo <- function(df, sourceKey, targetKey) {
  
  # count source ids
  uniqueSource <- length(unique(df[, sourceKey]))
  # remove NAs
  notMapped <- is.na(df[, targetKey])
  df <- df[!notMapped, ]
  # count multimapping IDs
  sourceMultiMap <- sum(table(df[, sourceKey]) >= 2)
  targetMultiMap <- sum(table(df[, targetKey]) >= 2)
  # count target ids
  uniqueTarget <- length(unique(df[, targetKey]))
  # print info
  m <- paste("------------------------------------------------",
             paste0(sum(notMapped), " of ", uniqueSource, " input ids on the translator data frame could not be mapped."),
             paste0(sourceMultiMap, " of ", uniqueSource, " input ids on the translator data frame were mapped to 2 or more target ids."),
             paste0(targetMultiMap, " of ", uniqueTarget, " target ids on the translator data frame were mapped to 2 or more input ids."),
             "------------------------------------------------",
             paste0("Input keys were finally mapped to ", uniqueTarget, " target ids."),
             "------------------------------------------------",
             sep = "\n")
  message(m)
  
}

#' Validate translator data frame
#'
#' Evaluates integrity of the translator data frame, subsetting it
#' to selected columns and removing duplicated rows.
#'
#' @param df The translator data frame to validate.
#' @param sourceKey Character indicating the column of the df where the source ids are stored.
#' @param targetKey Character indicating the column of the df where the source ids are stored.
#'
#' @return The validated and ammended translator data frame if no problems are found.
#'
validateTranslatorDf <- function(df, sourceKey, targetKey) {
  
  # subset to desired columns
  translatorDf <- df[ , c(sourceKey, targetKey)]
  # remove duplicated rows in source to target df
  uniqueDf <- unique(translatorDf)
  if(nrow(uniqueDf) != nrow(translatorDf)) {
    message("Translator data frame contains duplicated rows and will be removed.")
  }
  # print mapping info
  messageMappingInfo(uniqueDf, sourceKey, targetKey)
  # remove not mapped ids
  mapped <- !is.na(uniqueDf[, targetKey])
  uniqueDf <- uniqueDf[mapped, ]
  return(uniqueDf)
  
}