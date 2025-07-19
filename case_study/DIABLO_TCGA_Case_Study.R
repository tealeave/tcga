## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(123) # for reproducibility, remove for normal use


## -------------------------------------------------------------------------------------------------------------------
data(breast.TCGA) # load in the data

data = list(miRNA = breast.TCGA$data.train$mirna, # set a list of all the X dataframes
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)

lapply(data, dim) # check their dimensions

Y = breast.TCGA$data.train$subtype # set the response variable as the Y dataframe
summary(Y)


## ---- fig.show = "hold", out.width = "33%", fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the breast TCGA data. Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5. "----
list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

pls1 <- spls(data[["miRNA"]], data[["mRNA"]], keepX = list.keepX, keepY = list.keepY) # generate three pairwise PLS models
pls2 <- spls(data[["miRNA"]], data[["proteomics"]], keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["mRNA"]], data[["proteomics"]], keepX = list.keepX, keepY = list.keepY)

plotVar(pls1, cutoff = 0.5, title = "(a) miRNA vs mRNA", legend = c("miRNA", "mRNA"), # plot features of first PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(pls2, cutoff = 0.5, title = "(a) miRNA vs proteomics", legend = c("miRNA", "proteomics"), # plot features of second PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(pls3, cutoff = 0.5, title = "(a) mRNA vs proteomics", legend = c("mRNA", "proteomics"), # plot features of third PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))


## ---- echo = FALSE--------------------------------------------------------------------------------------------------
pls1 <- spls(data[["miRNA"]], data[["mRNA"]], ncomp = 1, keepX = 25, keepY = 25)
pls2 <- spls(data[["miRNA"]], data[["proteomics"]], ncomp = 1, keepX = 25, keepY = 25)
pls3 <- spls(data[["mRNA"]], data[["proteomics"]], ncomp = 1, keepX = 25, keepY = 25)


## -------------------------------------------------------------------------------------------------------------------
cor(pls1$variates$X, pls1$variates$Y) # calculate correlation of miRNA and mRNA
cor(pls2$variates$X, pls2$variates$Y) # calculate correlation of miRNA and proteins
cor(pls3$variates$X, pls3$variates$Y) # calculate correlation of mRNA and proteins


## -------------------------------------------------------------------------------------------------------------------
design = matrix(0.1, ncol = length(data), nrow = length(data), # for square matrix filled with 0.1s
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

design


## -------------------------------------------------------------------------------------------------------------------
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) # form basic DIABLO model


## ---- fig.cap = "FIGURE 2: Choosing the number of components in `block.plsda` using `perf()` with 10 × 10-fold CV function in the `breast.TCGA` study. Classification error rates (overall and balanced, see Section 7.3) are represented on the y-axis with respect to the number of components on the x-axis for each prediction distance presented in PLS-DA"----
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', folds = 10, nrepeat = 10) # run component number tuning with repeated CV

plot(perf.diablo) # plot output of tuning


## -------------------------------------------------------------------------------------------------------------------
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # set the optimal ncomp value
perf.diablo$choice.ncomp$WeightedVote # show the optimal choice for ncomp for each dist metric


## ---- eval=FALSE, include = FALSE-----------------------------------------------------------------------------------
## test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
## 
## t1 = proc.time()
## tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
##                               test.keepX = test.keepX, design = design,
##                               validation = 'Mfold', folds = 10, nrepeat = 1,
##                               cpus = 2, dist = "centroids.dist")
## t2 = proc.time()
## running_time = t2 - t1; running_time
## 
## list.keepX = tune.TCGA$choice.keepX
## list.keepX
## 
## save(tune.TCGA,list.keepX, file = 'RData/result-TCGA-diablo_design0.1.RData')


## ---- echo = FALSE--------------------------------------------------------------------------------------------------
load('RData/result-TCGA-diablo_design0.1.RData')


## ---- include =TRUE, eval=FALSE-------------------------------------------------------------------------------------
## test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)), # set grid of values for each component to test
##                    miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
##                    proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))
## 
## tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, # run the feature selection tuning
##                               test.keepX = test.keepX, design = design,
##                               validation = 'Mfold', folds = 10, nrepeat = 1,
##                               dist = "centroids.dist")


## -------------------------------------------------------------------------------------------------------------------
list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
list.keepX


## -------------------------------------------------------------------------------------------------------------------
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, # set the optimised DIABLO model
                          keepX = list.keepX, design = design)


## -------------------------------------------------------------------------------------------------------------------
final.diablo.model$design # design matrix for the final model


## -------------------------------------------------------------------------------------------------------------------
selectVar(final.diablo.model, block = 'mRNA', comp = 1)$mRNA$name # the features selected to form the first component


## ---- fig.cap = "FIGURE 3: Diagnostic plot from multiblock sPLS-DA applied on the `breast.TCGA` study. Samples are represented based on the specified component (here `ncomp = 1`) for each data set (mRNA, miRNA and protein). Samples are coloured by breast cancer subtype and 95% confidence ellipse plots are represented."----
plotDiablo(final.diablo.model, ncomp = 1)


## ---- fig.cap = "FIGURE 4: Sample plot from multiblock sPLS-DA performed on the `breast.TCGA` study. The samples are plotted according to their scores on the first 2 components for each data set. Samples are coloured by cancer subtype"----
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, title = 'DIABLO Sample Plots')


## ---- fig.cap = "FIGURE 5: Arrow plot from multiblock sPLS-DA performed on the `breast.TCGA` study. The samples are projected into the space spanned by the first two components for each data set then overlaid across data sets."----
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, title = 'DIABLO')


## ---- fig.cap = "FIGURE 6: Correlation circle plot from multiblock sPLS-DA performed on the `breast.TCGA` study. Variable types are indicated with different symbols and colours, and are overlaid on the same plot."----
plotVar(final.diablo.model, var.names = FALSE, style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), col = c('darkorchid', 'brown1', 'lightgreen'))


## ---- fig.cap = "FIGURE 7: Circos plot from multiblock sPLS-DA performed on the `breast.TCGA` study. The plot represents the correlations greater than 0.7 between variables of different types, represented on the side quadrants"----
circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)


## ---- eval = TRUE, fig.cap = "FIGURE 8: Relevance network for the variables selected by multiblock sPLS-DA performed on the `breast.TCGA` study on component 1. Each node represents a selected with colours indicating their type. The colour of the edges represent positive or negative correlations"----
network(final.diablo.model, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)


## ----eval = FALSE---------------------------------------------------------------------------------------------------
## library(igraph)
## my.network = network(final.diablo.model, blocks = c(1,2,3),
##         color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
## write.graph(my.network$gR, file = "myNetwork.gml", format = "gml")


## ---- fig.cap = "FIGURE 9: Loading plot for the variables selected by multiblock sPLS-DA performed on the `breast.TCGA` study on component 1. The most important variables (according to the absolute value of their coefficients) are ordered from bottom to top. As this is a supervised analysis, colours indicate the class for which the median expression value is the highest for each feature"----
plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median')


## ---- eval = TRUE, fig.cap = "FIGURE 10: Clustered Image Map for the variables selected by multiblock sPLS-DA performed on the `breast.TCGA` study on component 1. By default, Euclidean distance and Complete linkage methods are used. The CIM represents samples in rows (indicated by their breast cancer subtype on the left hand side of the plot) and selected features in columns (indicated by their data type at the top of the plot)."----
cimDiablo(final.diablo.model)


## -------------------------------------------------------------------------------------------------------------------
perf.diablo = perf(final.diablo.model, validation = 'Mfold', M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') # run repeated CV performance evaluation

perf.diablo$MajorityVote.error.rate

perf.diablo$WeightedVote.error.rate


## ---- fig.cap = "FIGURE 11: ROC and AUC based on multiblock sPLS-DA performed on the `breast.TCGA` study for the miRNA data set after 2 components. The function calculates the ROC curve and AUC for one class vs. the others."----
auc.splsda = auroc(final.diablo.model, roc.block = "miRNA", roc.comp = 2, print = FALSE)


## -------------------------------------------------------------------------------------------------------------------
data.test.TCGA = list(mRNA = breast.TCGA$data.test$mrna,
                      miRNA = breast.TCGA$data.test$mirna)

predict.diablo = predict(final.diablo.model, newdata = data.test.TCGA)


## -------------------------------------------------------------------------------------------------------------------
confusion.mat = get.confusion_matrix(truth = breast.TCGA$data.test$subtype,
                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat


## -------------------------------------------------------------------------------------------------------------------
get.BER(confusion.mat)


## -------------------------------------------------------------------------------------------------------------------
#sessionInfo()


## ---- include = FALSE-----------------------------------------------------------------------------------------------
# extract R code
#purl("DIABLO_TCGA.Rmd")

