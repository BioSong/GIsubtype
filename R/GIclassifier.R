#' description >> Function for the classification of GI cancer subtypes
#' argument
#' @param Expr a dataframe with log2_scaled Gene Expression Profiles data values, samples in columns, genes in rows
#' @param idType the type of geneid "SYMBOL", "ENSEMBL", "ENTREZID"
#' @param verbose TRUE/FALSE
#' @importFrom utils data
#' @export
GIclassifier <- function(Expr,idType = c("SYMBOL", "ENSEMBL", "ENTREZID"), verbose = TRUE) {
	  if (isTRUE(verbose)) {
        message("Checking input dataset and parameters......")
    }
    if (!is.matrix(Expr) & !is.data.frame(Expr)) {
        stop("Only gene expression profile in dataframe or matrix format is accepted.")
    }
    if (is.null(rownames(Expr)) | is.null(colnames(Expr))) {
        stop("Rownames and colnames are madatory in gene expression profile.")
    }
    if (sum(apply(Expr, 2, is.numeric)) != ncol(Expr)) {
        stop("Only numeric values in gene expression profile is accepted.")
    }
	  if (any(is.na(Expr))) {
        stop("Gene expression profile cannot contain any NA value(s).")
    }
    if (sum(c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ") %in%
        colnames(Expr)) != 0) {
        stop("Sample names in expression profile should not contain \"SYMBOL\", \"ENSEMBL\", \"ENTREZID\" and \"REFSEQ\".")
    }
    if (is.matrix(Expr)) {
        message("Transforming the input gene expression profile into dataframe format.")
        Expr <- as.data.frame(Expr)
    }
	  if (is.null(idType) | length(idType) != 1 | !(idType %in% c("SYMBOL", "ENSEMBL", "ENTREZID"))) {
        stop("idType should be one of \"SYMBOL\", \"ENSEMBL\", \"ENTREZID\".")
    }
    if (idType != "SYMBOL") {
        message("Converting gene ids to gene symbols......")
    }
    if (max(Expr)>50) {
        message("Performing log2 transformation......")
        Expr <- log2(Expr + 1)
    }
	  if (isTRUE(verbose)) {
        message("Making molecular subtype prediction......")
    }
	  utils::data(TrainExp)
	  genes <- intersect(rownames(TrainExp),rownames(Expr))
	  Expr <- Expr[match(genes,rownames(Expr)),]
	  TrainExp <- TrainExp[match(genes,rownames(TrainExp)),]
    res <- clusterRepro::IGP.clusterRepro(Expr, TrainExp)
	  return(res)
}


