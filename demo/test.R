library(TFpredict)
target_gene<-c("GATA1")
species<-c("human")
promoter_start<-2000
promoter_end<-200

TF<-TFpredict(target_gene=target_gene, species=species,romoter_start=promoter_start, promoter_end=promoter_end)
