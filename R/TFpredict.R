#' Predict transcription factors binding to promoter regions
#'
#' @param target_gene A character vector of gene symbols
#' @param species A character string specifying species ("human" or "mouse")
#' @param promoter_start Distance upstream of TSS (default: 2000)
#' @param promoter_end Distance downstream of TSS (default: 200)
#' @param min_score Minimum score threshold for PWM matching (default: "85%")
#' @return A data.frame with predicted transcription factors and their binding sites
#' @export
#' @examples
#' TF <- TFpredict(target_gene = "GATA1", species = "human",promoter_start = 2000, promoter_end = 200)
TFpredict <- function(target_gene, species = c("human", "mouse"),
                      promoter_start = 2000, promoter_end = 200,
                      min_score = "85%") {

  # 获取启动子序列（返回DNAStringSet）
  promoter_seqs <- get_promoter_sequences(target_gene, species,
                                          promoter_start, promoter_end)

  # 转换为单个DNAString（取第一个序列）
  if (length(promoter_seqs) == 0) {
    stop("No promoter sequences found for the target genes")
  }
  promoter_seq <- promoter_seqs[[1]]  # 提取第一个序列作为DNAString

  # 获取JASPAR数据库矩阵
  species_code <- switch(species,
                         human = 9606,
                         mouse = 10090)

  opts <- list(species = species_code, collection = "CORE")
  pfm_list <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts)

  # 转换为PWM
  pwm_list <- lapply(pfm_list, function(pfm) TFBSTools::toPWM(pfm))

  # 查找TF结合位点
  hits <- lapply(pwm_list, function(pwm) {
    if (methods::is(pwm, "PWMatrix")) {
      tryCatch({
        Biostrings::matchPWM(pwm@profileMatrix,
                             promoter_seq,  # 使用DNAString而非DNAStringSet
                             min.score = min_score)
      }, error = function(e) NULL)
    } else {
      return(NULL)
    }
  })

  # 筛选非空结果
  hits <- hits[!sapply(hits, is.null)]
  hits <- hits[sapply(hits, length) > 0]

  if (length(hits) == 0) {
    message("No transcription factor binding sites found with the current threshold")
    return(NULL)
  }

  # 获取TF信息
  tf_names <- sapply(pfm_list, function(pfm) pfm@name)
  matched_tf_names <- tf_names[names(hits)]

  # 构建结果数据框
  result <- do.call(rbind, lapply(names(hits), function(tf_id) {
    tf_views <- hits[[tf_id]]
    sequences <- as.character(tf_views)
    data.frame(
      JASPAR_ID = tf_id,
      Gene_Symbol = matched_tf_names[tf_id],
      Count = length(sequences),
      Sequences = paste(sequences, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }))
  result <- result[order(-result$Count), , drop = FALSE]

  rownames(result) <- NULL

  # 添加序列信息
  attr(result, "promoter_sequence") <- as.character(promoter_seq)
  attr(result, "promoter_range") <- paste0("TSS -", promoter_start, " to +", promoter_end)

  return(result)
}
