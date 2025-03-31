# TFpredict.R
## 中文描述：给一个基因，然后提取对应的启动子序列，使用JARSPR进行预测上游的转录因子

预测基因启动子区域转录因子结合位点的R包

## 安装

```r
# 安装依赖
if (!require("remotes")) install.packages("remotes")

# 从GitHub安装
remotes::install_github("yudonglin506311858/TFpredict.R", dependencies = TRUE)
#软件安装
library(devtools)
install_github("yourusername/TFpredict.R")
library(TFpredict.R)

```

## 使用示例
```r
library(TFpredict)
result <- TFpredict("GATA1", "human", 2000, 200)

#示例
target_gene <- c("GATA1")
species <- c("human")
promoter_start <- "2000"
promoter_end <- "200"
min_score = "85%"
TF <- TFpredict(target_gene = target_gene, species = species,
                promoter_start = promoter_start, promoter_end = promoter_end,min_score = "85%")


```

## 结果示例
```r
   JASPAR_ID Gene_Symbol Count                                  Sequences
1   MA0130.1     ZNF354C     5 CTCCAC | ATCCAC | CTCCAC | CTCCAC | GTCCAC
2   MA1536.1       NR2C2     3             AAGGTCAC | TGGGTCAT | CAGGTCAG
3   MA0528.2      ZNF263     3 CTGGGGAGGAGG | GGAGGGAGGAGG | GGAGGGAGGAAG
4   MA0719.1      RHOXF1     2                        CTGATCCA | GTGATCCA
5   MA0766.2       GATA5     2                    GAGATAAGAA | CTGATAAGAC
6   MA0597.2       THAP1     2                GGGCAGGGCAAG | CCCCAGGGCAGG
7   MA0479.1       FOXH1     1                                ATCAATCCACA
8   MA0599.1        KLF5     1                                 GCCCCACCCA
9   MA0672.1      NKX2-3     1                                 GCCACTTAAG
10  MA0673.1      NKX2-8     1                                  CCACTTAAG
11  MA0674.1      NKX6-1     1                                   GTCATTAA
12  MA0675.1      NKX6-2     1                                   GTCATTAA
13  MA0132.2        PDX1     1                                   GTCATTAA
14  MA0738.1        HIC2     1                                  ATGCCCAGG
15  MA0741.1       KLF16     1                                GACCCGCCCCC
16  MA0801.1         MGA     1                                   GGGTGTGA
17  MA0810.1      TFAP2A     1                               TGCCCCAGGGCA
18  MA0811.1      TFAP2B     1                               TGCCCCAGGGCA
19  MA0524.2      TFAP2C     1                               TGCCCCAGGGCA
20  MA0914.1        ISL2     1                                   CCACTTAA
21  MA0103.3        ZEB1     1                                CTCACCTGCAA
22  MA1496.1       HOXA4     1                                   GTCATTAA
23  MA1499.1       HOXB4     1                                   GTCATTAA
24  MA1502.1       HOXB8     1                                   GTCATTAA
25  MA1504.1       HOXC4     1                                   GTCATTAA
26  MA1507.1       HOXD4     1                                   GTCATTAA
27  MA1522.1         MAZ     1                                GCCCCCTCCCC
28  MA1530.1      NKX6-3     1                                  GTCATTAAT
29  MA1535.1       NR2C1     1                                  CAAGGTCAC
30  MA0893.2        GSX2     1                                   GTCATTAA
31  MA0910.2       HOXD8     1                                   GTCATTAA
32  MA1635.1     BHLHE22     1                                 CACAGCTGGT
33  MA1648.1       TCF12     1                                CTCACCTGCAA
34  MA1728.1      ZNF549     1                               TCTGCTGCCCCA
35  MA0592.3       ESRRA     1                              TCCCAAGGTCACC
36  MA0039.4        KLF4     1                               TGCCCCACCCAC
37  MA1108.2        MXI1     1                                 ACCACATGCC
38  MA0508.3       PRDM1     1                                CTCTTTCTCCC
39  MA0829.2      SREBF1     1                             AGGTCACCTGATCC
40  MA1959.1        KLF7     1                                  GGGGCGTGG
41  MA1965.1         SP5     1                                 CCCCTCCCCT
42  MA2003.1      NKX2-4     1                                 CCACTCCACC
43  MA0493.2        KLF1     1                                  GGGGCGTGG
44  MA0723.2        VAX2     1                                   GTCATTAA
45  MA0740.2       KLF14     1                                  GGGGCGTGG
46  MA0742.2       KLF12     1                                  GGGGCGTGG
47  MA1498.2       HOXA7     1                                   GTCATTAA
48  MA1511.2       KLF10     1                                  GGGGCGTGG
```