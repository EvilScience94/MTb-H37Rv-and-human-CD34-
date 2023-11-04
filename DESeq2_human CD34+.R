> dds <- DESeqDataSetFromMatrix(countData = CD_pro, colData = CD_col, design = ~ Condition + Time)
> dds <- DESeq(dds)
estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 12 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
> show(dds)
class: DESeqDataSet 
dim: 21504 81 
metadata(1): version
assays(6): counts mu ... replaceCounts replaceCooks
rownames(21504): HUM_01 HUM_02 ... HUM_21503 HUM_21504
rowData names(31): baseMean baseVar ... maxCooks replace
colnames(81): GSM3745998 GSM3745998.1 ... GSM3746014.1 GSM3746014.2
colData names(5): X Condition Time sizeFactor replaceable
> res <- results(dds)
> summary(res)
out of 18653 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4376, 23%
LFC < 0 (down)     : 5434, 29%
outliers [1]       : 0, 0%
low counts [2]     : 1794, 9.6%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
##LFC shrinkage analysis
> res_lfc <- lfcShrink(dds, coef="Condition_Uninfected_vs_Infected", type="apeglm")

#Multi-Factor Analysis
> colData(dds)
DataFrame with 81 rows and 5 columns
X  Condition     Time        sizeFactor replaceable
<factor>   <factor> <factor>         <numeric>   <logical>
  GSM3745998   GSM3745998 Uninfected      24h  0.86498230524255        TRUE
GSM3745998.1 GSM3745998 Uninfected      24h 0.851073064203406        TRUE
GSM3745998.2 GSM3745998 Uninfected      24h 0.853090935747772        TRUE
GSM3745998.3 GSM3745998 Uninfected      24h 0.871906684525629        TRUE
GSM3745998.4 GSM3745998 Uninfected      24h 0.875610399224274        TRUE
...                 ...        ...      ...               ...         ...
GSM3746013.1 GSM3746013 Uninfected     120h  1.79953557984354        TRUE
GSM3746013.2 GSM3746013 Uninfected     120h  1.81479959978726        TRUE
GSM3746014   GSM3746014   Infected     120h  1.40697628366602        TRUE
GSM3746014.1 GSM3746014   Infected     120h  1.34820963413772        TRUE
GSM3746014.2 GSM3746014   Infected     120h  1.35891863191325        TRUE
> ddsMF <- dds
> levels(ddsMF$Time)
[1] "120h" "24h"  "72h" 
> design(ddsMF) <- formula(~ Time + Condition)
> ddsMF <- DESeq(ddsMF)
using pre-existing size factors
estimating dispersions
found already estimated dispersions, replacing these
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
-- replacing outliers and refitting for 12 genes
-- DESeq argument 'minReplicatesForReplace' = 7 
-- original counts are preserved in counts(dds)
estimating dispersions
fitting model and testing
> resMF <- results(ddsMF)
> head(resMF)
log2 fold change (MLE): Condition Uninfected vs Infected 
Wald test p-value: Condition Uninfected vs Infected 
DataFrame with 6 rows and 6 columns
baseMean     log2FoldChange              lfcSE              stat               pvalue                padj
<numeric>          <numeric>          <numeric>         <numeric>            <numeric>           <numeric>
  HUM_01                  0                 NA                 NA                NA                   NA                  NA
HUM_02 0.0375790376895239  0.416071596556706   2.91711257724158 0.142631312827201    0.886581366250173                  NA
HUM_03 0.0375790376895239  0.416071596556706   2.91711257724158 0.142631312827201    0.886581366250173                  NA
HUM_04   9.24884522531863  0.133094829612702  0.178139335819907 0.747138912358226    0.454979713686539   0.573563694859862
HUM_05   422.396309151892 -0.325442804882407 0.0914909051963687 -3.55710553069622 0.000374963552245771 0.00123801629379244
HUM_06    37.964107677336  0.142445805423864 0.0846169320310048  1.68341964196563    0.092293883494193   0.157152728624817
> summary(resMF)

out of 18653 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 4835, 26%
LFC < 0 (down)     : 3947, 21%
outliers [1]       : 0, 0%
low counts [2]     : 2511, 13%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

##Multi-Factor analysis at LFC=1.5
> resMFT1_lfc1.5 <- results(ddsMF, 
                            contrast = c("Time", "24h", "72h"),
                            lfcThreshold = 1.5,
                            altHypothesis = "greaterAbs", alpha = 0.1)                          
> summary(resMFT1_lfc1.5)

out of 18653 with nonzero total read count
adjusted p-value < 0.1
LFC > 1.50 (up)    : 0, 0%
LFC < -1.50 (down) : 23, 0.12%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> resMFT2_lfc1.5 <- results(ddsMF,
                            contrast = c("Time", "24h", "120h"), 
                            lfcThreshold = 1.5,
                            altHypothesis = "greaterAbs", alpha = 0.1)
> summary(resMFT2_lfc1.5)

out of 18653 with nonzero total read count
adjusted p-value < 0.1
LFC > 1.50 (up)    : 25, 0.13%
LFC < -1.50 (down) : 228, 1.2%
outliers [1]       : 0, 0%
low counts [2]     : 1077, 5.8%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> resMFT3_lfc1.5 <- results(ddsMF, contrast = c("Time", "72h", "120h"),
                            lfcThreshold = 1.5,
                            altHypothesis = "greaterAbs", alpha = 0.1)
> summary(resMFT3_lfc1.5)

out of 18653 with nonzero total read count
adjusted p-value < 0.1
LFC > 1.50 (up)    : 16, 0.086%
LFC < -1.50 (down) : 23, 0.12%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

> resMFcond_lfc1.5 <- results(ddsMF, 
                              contrast = c("Condition", "Uninfected",
                                           "Infected"),
                              lfcThreshold = 1.5,
                              altHypothesis = "greaterAbs", alpha = 0.1)
> summary(resMFcond_lfc1.5)

out of 18653 with nonzero total read count
adjusted p-value < 0.1
LFC > 1.50 (up)    : 6, 0.032%
LFC < -1.50 (down) : 61, 0.33%
outliers [1]       : 0, 0%
low counts [2]     : 1794, 9.6%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

##Plotting of MF DE analysis at LFC=1.5
> par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
> ylim <- c(-5.0,5.0)
> plotMA(resMFcond_lfc1.5, ylim=ylim,
         main= "LFC = 1.5, Uninfected vs Infected"); drawLines1()
> plotMA(resMFT1_lfc1.5, ylim=ylim,
         main= "LFC = 1.5, 24h vs 72h"); drawLines1()
> plotMA(resMFT2_lfc1.5, ylim=ylim,
         main= "LFC = 1.5, 24h vs 120h"); drawLines1()
> plotMA(resMFT3_lfc1.5, ylim=ylim,
         main= "LFC = 1.5, 72h vs 120h"); drawLines1()
