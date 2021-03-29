
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CROCS (Changepoints for a Range of ComplexitieS)

CROCS is an extension of the CROPS [(Haynes et
al. 2017)](https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1116445)
and sequential search [(Hocking et
al. 2018)](https://arxiv.org/abs/1810.00117) algorithm which allows to
compute all optimal changepoint segmentations of data sequences for all
penalty values accross a peak range. This package implements the CROPS
algorithm as well as segmentation models for peak calling. They have
been described in our study [(Liehrmann et
al. 2020)](https://arxiv.org/abs/2012.06848) in a genomic context. On
top of the gfpop method [(Runge et
al. 2020)](https://arxiv.org/abs/2002.03646), it proposes a framework
to design your Peak Caller. The user can choose among several noise
assumptions `CROCS::lossFactory`, transformations
`CROCS::transformationFactory`, peak shape assumptions
`CROCS::graphFactory` and peak start/end post-processing rules
`CROCS::postProcessingRuleFactory` to build a Peak Caller
`CROCS::peakCallerFactory`. The `CROCS::CROCS` output is the first step
to compute an objective function that can be optimized in a supervised
learning procedure. We give two examples of its use in the following
section.

<!-- badges: start -->

<!-- badges: end -->

## Quick Start

We install the package from Github:

``` r
#devtools::install_github("aLiehrmann/CROCS")
library(CROCS)
```

## CROCS for supervised segmentation procedure on ChIP-Seq data

``` r
library(CROCS, quietly =TRUE)
library(data.table, quietly =TRUE)
library(purrr, quietly =TRUE)
library(furrr, quietly =TRUE)
library(PeakError, quietly =TRUE)
library(penaltyLearning, quietly =TRUE)
plan(multisession, workers=20)
```

We prepare the data.

``` r
head(regions <- as.data.table(CROCS::regions))
#>             cell.type  sample.id chrom chromStart chromEnd annotation chunk
#> 1: skeletalMuscleCtrl McGill0019  chr3   10016881 10018028  peakStart     1
#> 2: skeletalMuscleCtrl McGill0036  chr3   10016881 10018028  peakStart     1
#> 3: skeletalMuscleCtrl McGill0037  chr3   10016881 10018028  peakStart     1
#> 4:   skeletalMuscleMD McGill0012  chr3   10016881 10018028  peakStart     1
#> 5:   skeletalMuscleMD McGill0013  chr3   10016881 10018028  peakStart     1
#> 6: skeletalMuscleCtrl McGill0019  chr3   10018488 10019581    peakEnd     1
head(counts <- as.data.table(CROCS::counts))
#>             cell.type  sample.id chromStart chromEnd coverage chunk
#> 1: skeletalMuscleCtrl McGill0019   10011731 10011809        1     1
#> 2: skeletalMuscleCtrl McGill0019   10011809 10011971        0     1
#> 3: skeletalMuscleCtrl McGill0019   10011971 10011972        2     1
#> 4: skeletalMuscleCtrl McGill0019   10011972 10011977        3     1
#> 5: skeletalMuscleCtrl McGill0019   10011977 10012071        6     1
#> 6: skeletalMuscleCtrl McGill0019   10012071 10012075        3     1
training_chunks <- 1:3
test_chunks <- 4
training_counts <- counts[chunk %in% training_chunks,]
test_counts <- counts[chunk %in% test_chunks,]
training_regions <- regions[chunk %in% training_chunks,]
test_regions <- regions[chunk %in% test_chunks,]
```

We compute machine learning features on the problems. These features
will be use latter in the suprevised learning procedure.

``` r
all_features <- as.data.table(featureMatrix(
    data.sequences = counts, 
    problem.vars = c("sample.id","chunk"), 
    data.var = "coverage"
), keep.rownames = T)
all_features <- data.table(cbind(all_features, unique(counts[,.(sample.id,chunk)])))
setkey(all_features,chunk,sample.id)
training_features <- all_features[unique(training_counts[,.(chunk,sample.id)])]
test_features <- all_features[unique(test_counts[,.(chunk,sample.id)])]
```

### Example 1: the unconstrained segmentation model with gaussian transformed noise and max-jump post-processing rule as Peak Caller

We define the Peak Caller.

``` r
std_g <- graphFactory(graph="std")
loss <- lossFactory(type="mean")
tr <- transformationFactory(transformation="anscombe_poisson")
rule <- postProcessingRuleFactory(rule="maxjump")

my_peak_caller <- peakCallerFactory(
  mygraph=std_g, 
  transformation_f=tr, 
  loss_f=loss,
  postProcessingRule_f=rule
) 
```

For each problem in the training set we compute all optimal models
between 1 and 9 peaks with the help of `CROCS::CROCS`. The choice of the
lower and upper peak bound is problem-specific.

``` r
training_models <- future_map(
  split(training_counts, by=c("chunk", "sample.id"), drop = TRUE),
  function(counts_p){
    y = counts_p$coverage
    w = counts_p$chromEnd - counts_p$chromStart
    models <- CROCS(
      data=y, 
      weights=w, 
      lower_bound_peak=1, 
      upper_bound_peak=9,
      solver=my_peak_caller
    )
    map(models, function(model_p){
      if (is.na(model_p$start[[1]])){
        model_p$chromStart <- 0
        model_p$chromEnd <- 0.1
      } else {
        model_p$chromStart <- counts_p$chromEnd[model_p$start]
        model_p$chromEnd <- counts_p$chromEnd[model_p$end]
      }
      model_p
    })
  }, 
  .options = furrr_options(seed=NULL)
)
print(cbind(
  chromStart = training_models[[1]][[8]]$chromStart,
  chromEnd = training_models[[1]][[8]]$chromEnd
))
#>      chromStart chromEnd
#> [1,]   10017758 10018745
#> [2,]   10027697 10029200
#> [3,]   10067192 10068519
```

Given the labels of each problem we compute the errors associated to
each model previously computed.

``` r
training_errors <- future_map(
  split(training_regions, by=c("chunk", "sample.id"), drop=TRUE), 
  function(region_p){
    map(
      training_models[[paste0(region_p$chunk[[1]],".",region_p$sample.id[[1]])]], 
      function(model_p){
        peaks <- data.frame(
          chromStart=model_p$chromStart, 
          chromEnd=model_p$chromEnd
        )
        errors <- PeakErrorChrom(peaks=peaks, regions=region_p)
        errors$changepoints <- model_p$changepoints
        errors
      }
    )
  },
  .options = furrr_options(seed=NULL)
)
```

We merge the computed models and associated errors.

``` r
training_models_f <- rbindlist(pmap(
  unique(counts[,.(sample.id,chunk)]), 
  function(sample.id, chunk){
    rbindlist(map(
      training_models[[paste0(chunk,".", sample.id)]], 
      function(model_p){
        data.table(
          chunk=chunk,
          sample.id=sample.id,
          changepoints=model_p$changepoints,
          loss=model_p$loss
        )
      }
    ))
  }
))
head(training_models_f)
#>    chunk  sample.id changepoints      loss
#> 1:     1 McGill0019           39  47644.08
#> 2:     1 McGill0019            0 518730.26
#> 3:     1 McGill0019            8 130621.41
#> 4:     1 McGill0019           20  66944.04
#> 5:     1 McGill0019            4 203731.92
#> 6:     1 McGill0019           27  55586.47
training_errors_f <- rbindlist(pmap(
  unique(counts[,.(sample.id,chunk)]), 
  function(sample.id, chunk){
    rbindlist(map(
      training_errors[[paste0(chunk,".", sample.id)]], 
      function(errors_p){
        data.table(
          chunk=chunk,
          sample.id=sample.id,
          changepoints=errors_p$changepoints,
          errors=sum(errors_p$fn+errors_p$fp),
          labels=nrow(errors_p)
        )
      }
    ))
  }
))
head(training_errors_f)
#>    chunk  sample.id changepoints errors labels
#> 1:     1 McGill0019           39      0      7
#> 2:     1 McGill0019           39      0      7
#> 3:     1 McGill0019           39      0      7
#> 4:     1 McGill0019           39      0      7
#> 5:     1 McGill0019           39      0      7
#> 6:     1 McGill0019           39      0      7
dt_model_error <- merge(
  x = training_models_f, 
  y = training_errors_f, 
  by=c("chunk", "sample.id", "changepoints")
)
head(dt_model_error)
#>    chunk  sample.id changepoints     loss errors labels
#> 1:     1 McGill0019            0 518730.3      6      7
#> 2:     1 McGill0019            0 518730.3      6      7
#> 3:     1 McGill0019            0 518730.3      6      7
#> 4:     1 McGill0019            0 518730.3      6      7
#> 5:     1 McGill0019            0 518730.3      6      7
#> 6:     1 McGill0019            0 518730.3      6      7
```

We compute the piecewise constant objective function \(E(\lambda)\) for
each problem.

``` r
training_modelSelection <- dt_model_error[,modelSelection(
  models = .SD,
  loss = "loss",
  complexity = "changepoints"
), by=.(chunk, sample.id)]
head(training_modelSelection)
#>    chunk  sample.id min.lambda max.lambda min.log.lambda max.log.lambda
#> 1:     1 McGill0019     0.0000   461.2570           -Inf       6.133955
#> 2:     1 McGill0019   461.2570   474.3379       6.133955       6.161920
#> 3:     1 McGill0019   474.3379   571.9229       6.161920       6.349004
#> 4:     1 McGill0019   571.9229   572.4606       6.349004       6.349944
#> 5:     1 McGill0019   572.4606   673.6850       6.349944       6.512763
#> 6:     1 McGill0019   673.6850   777.7755       6.512763       6.656438
#>    cum.iterations changepoints     loss errors labels
#> 1:             28           39 47644.08      0      7
#> 2:             27           37 48566.59      0      7
#> 3:             26           36 49040.93      0      7
#> 4:             25           34 50184.77      0      7
#> 5:             24           33 50757.23      0      7
#> 6:             23           32 51430.92      0      7
```

We search the problem-specific interval of \(\lambda\) associated with
minimum error.

``` r
training_targets <- targetIntervals(
  models = training_modelSelection, 
  problem.vars = c("chunk", "sample.id")
)
head(training_targets)
#>    chunk  sample.id min.log.lambda max.log.lambda errors
#> 1:     1 McGill0019           -Inf       8.302645      0
#> 2:     1 McGill0036           -Inf       8.661698      0
#> 3:     1 McGill0037           -Inf       8.575239      0
#> 4:     1 McGill0012           -Inf       8.677474      0
#> 5:     1 McGill0013           -Inf       9.262656      0
#> 6:     2 McGill0019       7.417167      10.824984      0
```

We learn a function that predicts problem-specific \(\log(\lambda)\)
values.

``` r
training_targets <- training_targets[order(chunk, sample.id)]
training_features <- training_features[order(chunk, sample.id)]
training_targets_mat <- as.matrix(
  training_targets[ ,.(min.log.lambda, max.log.lambda)]
)
training_features_mat <- as.matrix(
  training_features[,!c("rn", "chunk", "sample.id")]
)
training_features_mat <- training_features_mat[
  is.finite(training_targets_mat[,1])|is.finite(training_targets_mat[,2]),
]
training_targets_mat <- training_targets_mat[
  is.finite(training_targets_mat[,1])|is.finite(training_targets_mat[,2]),
]
set.seed(2020)
fold.vec = sample(rep(1:4, l = nrow(training_features_mat)))
lambda_function <- IntervalRegressionCV(
  target.mat = training_targets_mat,
  feature.mat = training_features_mat,
  LAPPLY = future_map,
  fold.vec = fold.vec
)
#> Loading required namespace: directlabels
#> Loading required namespace: directlabels
#> Loading required namespace: directlabels
#> Loading required namespace: directlabels
```

We predict \(\lambda\) values for each problem from the validation set
or unlabeled problems.

``` r
predicted_log_lambda <- lambda_function$predict(
  as.matrix(test_features[,!c("rn", "chunk", "sample.id")])
)
predicted_lambda <- data.table(
  chunk=test_features$chunk,
  sample.id=test_features$sample.id,
  lambda=exp(predicted_log_lambda[,1])
)
head(predicted_lambda)
#>    chunk  sample.id   lambda
#> 1:     4 McGill0023 1863.090
#> 2:     4 McGill0022 2338.691
#> 3:     4 McGill0267 1507.355
#> 4:     4 McGill0019 2198.152
#> 5:     4 McGill0036 2429.329
```

Using previously computed problem-specific \(\lambda\) values, we
retreive peaks associated to each problem from the validation set or
unlabeled problems.

``` r
split_test_counts <- split(test_counts, by=c("chunk", "sample.id"),drop=TRUE)
predicted_models <- pmap(predicted_lambda, function(chunk, sample.id, lambda){
  current_pb <- split_test_counts[[paste0(chunk,".",sample.id)]]
  y = current_pb$coverage
  w = current_pb$chromEnd - current_pb$chromStart
  predicted_model <- my_peak_caller(
    data = y,
    weights = w,
    lambda = lambda
  )
  if (is.na(predicted_model$start[[1]])){
    predicted_model$chromStart <- 0
    predicted_model$chromEnd <- 0.1
  } else {
    predicted_model$chromStart <- current_pb$chromEnd[predicted_model$start]
    predicted_model$chromEnd <- current_pb$chromEnd[predicted_model$end]
  }
  predicted_model
})
names(predicted_models) <- paste0(
  predicted_lambda$chunk,
  ".",
  predicted_lambda$sample.id
)
names(predicted_models)[[3]]
#> [1] "4.McGill0267"
print(cbind(
  chromStart = predicted_models[[3]]$chromStart,
  chromEnd = predicted_models[[3]]$chromEnd
))
#>      chromStart chromEnd
#> [1,]   38060802 38061488
#> [2,]   38155074 38156798
```

Since we have labels for the validation set, we can compute the
accuracy.

``` r
predicted_errors <- rbindlist(future_map(
  split(test_regions, by=c("chunk", "sample.id"), drop=TRUE), 
  function(region_p){
    model <- predicted_models[[paste0(region_p$chunk[[1]],".",region_p$sample.id[[1]])]]
    peaks <- data.frame(
      chromStart=model$chromStart, 
      chromEnd=model$chromEnd
    )
    errors <- PeakErrorChrom(peaks=peaks, regions=region_p)
    errors
  },
  .options = furrr_options(seed=NULL)
))
print(paste0(1-predicted_errors[,sum(fp+fn)/.N]," accuracy"))
#> [1] "0.95 accuracy"
```

``` r
plan(sequential)
```

### Example 2: the unconstrained segmentation model with negative binomial noise and max-jump post-processing rule as Peak Caller

``` r
plan(multisession, workers=20)
```

We define the Peak Caller.

``` r
std_g <- graphFactory(graph="std")
loss <- lossFactory(type="negbin")
tr <- transformationFactory(transformation="raw")
rule <- postProcessingRuleFactory(rule="maxjump")

my_peak_caller <- peakCallerFactory(
  mygraph=std_g, 
  transformation_f=tr, 
  loss_f=loss,
  postProcessingRule_f=rule
) 
```

For each problem in the training set we compute all optimal models
between 1 and 9 peaks with the help of `CROCS::CROCS`. The choice of the
lower and upper peak bound is problem-specific. In order to calibrate
the dispersion \(\phi\) parameter from the negative binomial
distribution, we repeate this step for a range of \(\phi\) values.

``` r
phi <- exp(seq(log(1),log(1000),length=8))
training_models <- future_map(
  split(training_counts, by=c("chunk", "sample.id"), drop=TRUE),
  function(counts_p){
    map(phi, function(phi_p){
      y = counts_p$coverage
      w = counts_p$chromEnd - counts_p$chromStart
      models <- CROCS(
        data=y, 
        weights=w, 
        lower_bound_peak=1, 
        upper_bound_peak=9,
        solver=my_peak_caller,
        phi=phi_p
      )
      map(models, function(model_p){
        if (is.na(model_p$start[[1]])){
          model_p$chromStart <- 0
          model_p$chromEnd <- 0.1
        } else {
          model_p$chromStart <- counts_p$chromEnd[model_p$start]
          model_p$chromEnd <- counts_p$chromEnd[model_p$end]
        }
        model_p
      })
    })
  }, 
  .options = furrr_options(seed=NULL)
)
```

Given the labels of each problem we compute the errors associated to
each model previously computed.

``` r
training_errors <- future_map(
  split(training_regions, by=c("chunk", "sample.id"), drop=TRUE), 
  function(region_p){
    map(
      training_models[[paste0(region_p$chunk[[1]],".",region_p$sample.id[[1]])]], 
      function(phi_p){
        map(phi_p, function(model_p){
          peaks <- data.frame(
            chromStart=model_p$chromStart, 
            chromEnd=model_p$chromEnd
          )
          errors <- PeakErrorChrom(peaks=peaks, regions=region_p)
          errors$changepoints <- model_p$changepoints
          errors$phi <- model_p$phi
          errors
        })
      }
    )
  },
  .options = furrr_options(seed=NULL)
)
```

We merge the computed models and associated errors.

``` r
# /// compute lambda function
training_models_f <- rbindlist(future_pmap(
  unique(counts[,.(sample.id,chunk)]), 
  function(sample.id, chunk){
    rbindlist(map(
      training_models[[paste0(chunk,".", sample.id)]],
      function(phi_p){
        rbindlist(map(phi_p, 
          function(model_p){
            data.table(
              chunk=chunk,
              sample.id=sample.id,
              changepoints=model_p$changepoints,
              loss=model_p$loss,
              phi=model_p$phi
            )
          }
        ))
      }
    ))
  }
))
head(training_models_f)
#>    chunk  sample.id changepoints     loss phi
#> 1:     1 McGill0019           26 162576.2   1
#> 2:     1 McGill0019            0 229269.5   1
#> 3:     1 McGill0019            8 170538.6   1
#> 4:     1 McGill0019           14 166832.7   1
#> 5:     1 McGill0019            4 189953.9   1
#> 6:     1 McGill0019           19 164908.2   1
training_errors_f <- rbindlist(future_pmap(
  unique(counts[,.(sample.id,chunk)]), 
  function(sample.id, chunk){
    rbindlist(map(
      training_errors[[paste0(chunk,".", sample.id)]],
      function(phi_p){
        rbindlist(map(phi_p, 
          function(errors_p){
            data.table(
              chunk=chunk,
              sample.id=sample.id,
              changepoints=errors_p$changepoints,
              phi=errors_p$phi,
              errors=sum(errors_p$fn+errors_p$fp),
              labels=nrow(errors_p)
            )
          }
        ))
      }
    ))
  }
))
head(training_errors_f)
#>    chunk  sample.id changepoints phi errors labels
#> 1:     1 McGill0019           26   1      5      7
#> 2:     1 McGill0019           26   1      5      7
#> 3:     1 McGill0019           26   1      5      7
#> 4:     1 McGill0019           26   1      5      7
#> 5:     1 McGill0019           26   1      5      7
#> 6:     1 McGill0019           26   1      5      7
dt_model_error <- merge(
  x = training_models_f, 
  y = training_errors_f, 
  by=c("chunk", "sample.id", "changepoints","phi")
)
head(dt_model_error)
#>    chunk  sample.id changepoints phi     loss errors labels
#> 1:     1 McGill0019            0   1 229269.5      6      7
#> 2:     1 McGill0019            0   1 229269.5      6      7
#> 3:     1 McGill0019            0   1 229269.5      6      7
#> 4:     1 McGill0019            0   1 229269.5      6      7
#> 5:     1 McGill0019            0   1 229269.5      6      7
#> 6:     1 McGill0019            0   1 229269.5      6      7
```

We compute the piecewise constant objective function \(E(\lambda)\) for
each problem and each \(\phi\) value.

``` r
training_modelSelection <- dt_model_error[,modelSelection(
  models = .SD,
  loss = "loss",
  complexity = "changepoints"
), by=.(chunk, sample.id, phi)]
head(training_modelSelection)
#>    chunk  sample.id phi min.lambda max.lambda min.log.lambda max.log.lambda
#> 1:     1 McGill0019   1     0.0000   325.9183           -Inf       5.786647
#> 2:     1 McGill0019   1   325.9183   332.1626       5.786647       5.805624
#> 3:     1 McGill0019   1   332.1626   340.2916       5.805624       5.829803
#> 4:     1 McGill0019   1   340.2916   349.7040       5.829803       5.857087
#> 5:     1 McGill0019   1   349.7040   358.9196       5.857087       5.883098
#> 6:     1 McGill0019   1   358.9196   392.1469       5.883098       5.971637
#>    cum.iterations changepoints     loss errors labels
#> 1:             16           26 162576.2      5      7
#> 2:             15           23 163553.9      5      7
#> 3:             14           21 164218.2      5      7
#> 4:             13           20 164558.5      4      7
#> 5:             12           19 164908.2      4      7
#> 6:             11           17 165626.1      4      7
```

We search the \(\phi\) value which minimizes the errors over all
problems of the training set.

``` r
training_modelSelection_by_phi <- split(training_modelSelection, by="phi", drop=TRUE)
dt <- rbindlist(map(training_modelSelection_by_phi, function(dt_by_phi) {
  sorted_lambda <- sort(c(dt_by_phi$min.log.lambda, dt_by_phi$max.log.lambda))
  grid <- rle(sorted_lambda)$values
  errors <- vector(mode="integer", length= length(grid)-1)
  walk(1:nrow(dt_by_phi), function(i){
    a <- which(grid>=dt_by_phi$min.log.lambda[[i]] & grid<dt_by_phi$max.log.lambda[[i]])
    errors[a] <<- errors[a] + dt_by_phi$errors[[i]]
  })
  i_min <- which.min(errors)
  predicted.lambda <- sum(grid[i_min:(i_min+1)])/2
  data.table(phi = dt_by_phi$phi[[1]], errors = min(errors))
}))
print(predicted_phi <- dt$phi[which.min(dt$errors)])
#> [1] 51.79475
```

After fixing \(\phi\), we search the problem-specific interval of
\(\lambda\) associated with minimum error.

``` r
training_targets <- targetIntervals(
  models = training_modelSelection[phi==predicted_phi], 
  problem.vars = c("chunk", "sample.id")
)
head(training_targets)
#>    chunk  sample.id min.log.lambda max.log.lambda errors
#> 1:     1 McGill0019           -Inf       4.712877      0
#> 2:     1 McGill0036       2.774878       4.845226      0
#> 3:     1 McGill0037           -Inf       5.021139      0
#> 4:     1 McGill0012       2.824390       4.889388      0
#> 5:     1 McGill0013           -Inf       5.303372      0
#> 6:     2 McGill0019       4.158716       7.968949      0
```

Learn a function that predicts problem-specific (()) values.

``` r
training_targets <- training_targets[order(chunk, sample.id)]
training_features <- training_features[order(chunk, sample.id)]
training_targets_mat <- as.matrix(
  training_targets[ ,.(min.log.lambda, max.log.lambda)]
)
training_features_mat <- as.matrix(
  training_features[,!c("rn", "chunk", "sample.id")]
)
training_features_mat <- training_features_mat[
  is.finite(training_targets_mat[,1])|is.finite(training_targets_mat[,2]),
]
training_targets_mat <- training_targets_mat[
  is.finite(training_targets_mat[,1])|is.finite(training_targets_mat[,2]),
]
set.seed(2020)
fold.vec = sample(rep(1:4, l = nrow(training_features_mat)))
lambda_function <- IntervalRegressionCV(
  target.mat = training_targets_mat,
  feature.mat = training_features_mat,
  LAPPLY = future_map,
  fold.vec = fold.vec
)
#> Loading required namespace: directlabels
#> Loading required namespace: directlabels
#> Loading required namespace: directlabels
#> Loading required namespace: directlabels
```

We predict \(\lambda\) values for each problem from the validation set
or unlabeled problems.

``` r
predicted_log_lambda <- lambda_function$predict(
  as.matrix(test_features[,!c("rn", "chunk", "sample.id")])
)
predicted_lambda <- data.table(
  chunk=test_features$chunk,
  sample.id=test_features$sample.id,
  lambda=exp(predicted_log_lambda[,1])
)
head(predicted_lambda)
#>    chunk  sample.id  lambda
#> 1:     4 McGill0023 75.2662
#> 2:     4 McGill0022 75.2662
#> 3:     4 McGill0267 75.2662
#> 4:     4 McGill0019 75.2662
#> 5:     4 McGill0036 75.2662
```

Using previously computed problem-specific \(\lambda\) values and the
consensus \(\phi\) value, we retreive peaks associated to each problem
from the validation set or unlabeled problems.

``` r
split_test_counts <- split(test_counts, by=c("chunk", "sample.id"), drop=TRUE)
predicted_models <- pmap(predicted_lambda, function(chunk, sample.id, lambda){
  current_pb <- split_test_counts[[paste0(chunk,".",sample.id)]]
  y = current_pb$coverage
  w = current_pb$chromEnd - current_pb$chromStart
  predicted_model <- my_peak_caller(
    data = y,
    weights = w,
    lambda = lambda,
    phi = predicted_phi
  )
  if (is.na(predicted_model$start[[1]])){
    predicted_model$chromStart <- 0
    predicted_model$chromEnd <- 0.1
  } else {
    predicted_model$chromStart <- current_pb$chromEnd[predicted_model$start]
    predicted_model$chromEnd <- current_pb$chromEnd[predicted_model$end]
  }
  predicted_model
})
names(predicted_models) <- paste0(
  predicted_lambda$chunk,
  ".",
  predicted_lambda$sample.id
)
```

Since we have labels for the validation set, we can compute the
accuracy.

``` r
predicted_errors <- rbindlist(future_map(
  split(test_regions, by=c("chunk", "sample.id"), drop=TRUE), 
  function(region_p){
    model <- predicted_models[[paste0(region_p$chunk[[1]],".",region_p$sample.id[[1]])]]
    peaks <- data.frame(
      chromStart=model$chromStart, 
      chromEnd=model$chromEnd
    )
    errors <- PeakErrorChrom(peaks=peaks, regions=region_p)
    errors
  },
  .options = furrr_options(seed=NULL)
))

print(paste0(1-predicted_errors[,sum(fp+fn)/.N]," accuracy"))
#> [1] "0.85 accuracy"
```

``` r
plan(sequential)
```
