#' Get promoter sequences for target genes
#'
#' @param target_gene A character vector of gene symbols
#' @param species A character string specifying species ("human" or "mouse")
#' @param promoter_start Distance upstream of TSS (default: 2000)
#' @param promoter_end Distance downstream of TSS (default: 200)
#' @return A DNAStringSet object containing promoter sequences
#' @export
#' @examples
#' promoter_seq <- get_promoter_sequences("GATA1", "human", 2000, 200)
get_promoter_sequences <- function(target_gene, species = c("human", "mouse"),
                                   promoter_start = 2000, promoter_end = 200) {

  species <- match.arg(species)

  # 检查并加载必要的包
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  # 设置物种特定的参数
  species_params <- switch(species,
                           human = list(
                             bsgenome = "BSgenome.Hsapiens.UCSC.hg19",
                             txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                             orgdb = "org.Hs.eg.db"
                           ),
                           mouse = list(
                             bsgenome = "BSgenome.Mmusculus.UCSC.mm10",
                             txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene",
                             orgdb = "org.Mm.eg.db"
                           )
  )

  # 安装必要的包（如果尚未安装）
  if (!requireNamespace(species_params$bsgenome, quietly = TRUE))
    BiocManager::install(species_params$bsgenome)
  if (!requireNamespace(species_params$txdb, quietly = TRUE))
    BiocManager::install(species_params$txdb)
  if (!requireNamespace(species_params$orgdb, quietly = TRUE))
    BiocManager::install(species_params$orgdb)

  # 加载包
  suppressPackageStartupMessages({
    library(species_params$bsgenome, character.only = TRUE)
    library(species_params$txdb, character.only = TRUE)
    library(species_params$orgdb, character.only = TRUE)
    library(GenomicFeatures)
    library(Biostrings)
  })

  # 获取基因组对象
  genome <- get(species_params$bsgenome)
  txdb <- get(species_params$txdb)
  orgdb <- get(species_params$orgdb)

  # 将基因符号转换为Entrez ID
  entrez_ids <- mapIds(orgdb,
                       keys = target_gene,
                       column = "ENTREZID",
                       keytype = "SYMBOL")

  if (any(is.na(entrez_ids))) {
    stop("Could not find Entrez IDs for: ",
         paste(target_gene[is.na(entrez_ids)], collapse = ", "))
  }

  # 获取基因的转录本信息
  transcripts <- transcriptsBy(txdb, by = "gene")[entrez_ids]

  # 提取每个基因的主要转录本（选择最长的转录本）
  primary_transcripts <- lapply(transcripts, function(tx) {
    tx[which.max(width(tx))]
  })

  # 创建空的DNAStringSet
  promoter_seqs <- DNAStringSet()

  for (i in seq_along(primary_transcripts)) {
    tx <- primary_transcripts[[i]]
    gene_symbol <- target_gene[i]

    # 确定链的方向
    if (as.character(strand(tx)) == "+") {
      prom_start <- max(1, start(tx) - promoter_start)  # 确保不超出染色体边界
      prom_end <- start(tx) + promoter_end
    } else {
      prom_start <- max(1, end(tx) - promoter_end)
      prom_end <- end(tx) + promoter_start
    }

    # 获取序列（直接返回DNAStringSet）
    seq <- getSeq(genome,
                  names = as.character(seqnames(tx)),
                  start = prom_start,
                  end = prom_end,
                  strand = "+")  # 统一使用正链获取序列

    # 转换为DNAStringSet并命名
    seq_set <- DNAStringSet(seq)
    names(seq_set) <- gene_symbol

    # 添加到结果集合
    promoter_seqs <- c(promoter_seqs, seq_set)
  }

  return(promoter_seqs)
}
