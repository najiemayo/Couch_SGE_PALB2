precisionBasedThreshold = function(score, label, weights, precisions) {
  # use precisions to choose the threshold for both sides
  
  # input
  # score: n*1 vector, the functional score
  # label: n*1 factor vector the label of the variant effects, 
  #   assume only two classes, the levels should be specifed in the order
  #   of increasing scores
  # weights: 2*1 vector weights for each class
  # precisions: a m*1 vector for different precision cutoffs
  
  # output
  # cutoffs: m*2 matrix corresponding to specified precisions for both classes  
  
  # sort scores
  indexSort = order(score, decreasing = F)
  score = score[indexSort]
  label = label[indexSort]
  
  # assign the precision values for each score for both classes
  uniqueValues = unique(score)
  k = length(uniqueValues)
  precision1 = numeric(length(uniqueValues)) 
  precision2 = numeric(length(uniqueValues))
  score1 = score[label == levels(label)[1]]
  score2 = score[label == levels(label)[2]]
  
  for (i in 1:k) { # not efficient, but easy to write
    # for class 1
    TP1 = sum(score1 <= uniqueValues[i]) * weights[1]
    FP1 = sum(score2 <= uniqueValues[i]) * weights[2]
    
    # for class 2
    TP2 = sum(score2 >= uniqueValues[i]) * weights[2]
    FP2 = sum(score1 >= uniqueValues[i]) * weights[1]
    
    precision1[i] = TP1 / (TP1 + FP1)
    precision2[i] = TP2 / (TP2 + FP2)
  }
  
  # find cutoffs for precisions
  cutoffs = matrix(NA, length(precisions), 2)
  rownames(cutoffs) = precisions
  colnames(cutoffs) = levels(label)
  for (i in 1:length(precisions)) {
    # for class1, choose the rightmost
    cutoffs[i, 1] = uniqueValues[max(which(precision1 >= precisions[i]))]
    # for class2, choose the leftmost
    cutoffs[i, 2] = uniqueValues[min(which(precision2 >= precisions[i]))]
  }
  
  return(cutoffs)
}

ROCOut = function(ROCResult, truth) {
  levelName = levels(truth)                 
  confusion = c(sum(ROCResult[[2]] == levelName[1] & truth == levelName[1]),
                sum(ROCResult[[2]] == levelName[2] & truth == levelName[1]),
                sum(ROCResult[[2]] == levelName[1] & truth == levelName[2]),
                sum(ROCResult[[2]] == levelName[2] & truth == levelName[2]))
  confusion = paste(confusion, collapse = " ")
  out = data.frame(ROCResult[[1]]$auc, ROCResult[[3]], ROCResult[[4]], 
                   ROCResult[[5]], confusion)
  colnames(out) = c("AUC", "cutoff", "Sensitivity", "Specificity", 
                    "ConfusionMatrix")
  return(out)
}

ROCBasedPrediction = function(score, truth) {
  # calculate AUC and the predicted labels based on the optimal cut
  # assume a binary classification
  
  # input
  # score: the score
  # truth: the class label, a factor with lower level corresponds to lower score
  
  roc_score = roc(predictor = score, response = truth,
                  levels = rev(levels(truth)))
  out = data.frame(functional.score = roc_score$"threshold", 
                   sens =  roc_score$"sensitivities", 
                   specs = roc_score$"specificities")
  
  out$avg = out$sens+out$specs
  index = which(out$avg == max(out$avg))
  if (length(index) > 1) {
    out$F1 = 1 / (1 / out$sens + 1 / out$specs)
    F1Scores = out$F1[index]
    index = index[which.max(F1Scores)] # only keep the first with the max F1
    
    if (length(index) > 1) {
      sensitivity = out$sens[index]
      index = index[which.max(sensitivity)]
    }
  }
  
  cutoff0 = out$functional.score[index]
  sen0 = out$sens[index]
  spec0 = out$specs[index]
  
  labels = levels(truth)
  predicted = ifelse(score < cutoff0, labels[1], labels[2])
  predicted = factor(predicted, levels = levels(truth))
  return(list(roc_score, predicted, cutoff0, sen0, spec0))
}

medianIRQ = function(data, referenceType) {
  # normalize the data using the minus median and divided by interquartile
  
  # input
  # data: a data frame
  # referenceType: a logical index selecting the rows that are used as the 
  #  reference points
  
  median = apply(X = data[referenceType, , drop = F], MARGIN = 2, FUN = median)
  quartiles = apply(X = data[referenceType, , drop = F], MARGIN = 2, 
                    FUN = quantile, probs = c(0.25, 0.75))
  IRQ = quartiles[2, ] - quartiles[1, ]
  normalized = data
  for (i in 1:dim(data)[2]) {
    normalized[, i] = (data[, i] - median[i]) / IRQ[i]
  }
  
  return(normalized)
}

median2Normalize = function(data, referenceType1, referenceType2) {
  # normalize the data using the median of two reference types
  
  # input
  # data: a data frame
  # referenceType1, referenceType2: 
  #    a logical index selecting the rows that are used as the 
  #    reference points
  
  D14libColumns = grepl("D14lib", colnames(data))
  median1 = apply(X = data[referenceType1, D14libColumns, drop = F], MARGIN = 2, 
                  FUN = median)
  median2 = apply(X = data[referenceType2, D14libColumns, drop = F], MARGIN = 2, 
                  FUN = median)
  medianMean1 = mean(median1)
  medianMean2 = mean(median2)
  
  nReplicate = length(median1)
  a = numeric(nReplicate)
  b = numeric(nReplicate)
  
  for (i in 1:nReplicate) {
    a[i] = (medianMean2 - medianMean1) / (median2[i] - median1[i])
    b[i] = medianMean2 - a * median2[i]
  }
  
  normalized = data
  
  replicateID = sub("_.*$", "", colnames(data))
  replicateNumber = as.integer(sub("R", "", replicateID))
  
  for (i in 1:dim(data)[2]) {
    normalized[, i] = a[replicateNumber[i]] * data[, i] + b[replicateNumber[i]]
  }
  
  return(normalized)
}

plotAcrossReplicates = function(data, nReplicate, featureID, levelsOfEvent) {
  dataForPlot = NULL
  for (i in 1:nReplicate) {
    names = paste0("R", i, "_", featureID)
    dataPerReplicate = data[, c(names, "EventType")]
    colnames(dataPerReplicate)[1] = "log2Ratio"
    replicate = paste0("R", i)
    dataPerReplicate = cbind(dataPerReplicate, replicate)
    dataForPlot = rbind(dataForPlot, dataPerReplicate)
  }
  g1 = ggplot(dataForPlot, aes(replicate, log2Ratio, 
                    fill = factor(EventType, levels = levelsOfEvent))) + 
    geom_violin()
  return(g1)
}

functionalClass = function(score, pathogenicT, benignT) {
  FSClass = character(length(score))
  FSClass[score <= pathogenicT] = "abnormal"
  FSClass[score >= benignT] = "normal"
  FSClass[score < benignT & score > pathogenicT] = 
    "I/U"
  
  FSClass = factor(FSClass, levels = c("normal", "I/U", "abnormal"))
  return(FSClass)
}

BayesDelNoAFClass = function(score) {
  FSClass = character(length(score))
  FSClass[] = "NA"

  FSClass[score >= 0.5] = "P_Strong"
  FSClass[score < 0.5 & score >= 0.27] = 
    "P_Moderate"
  FSClass[score < 0.27 & score >= 0.13] = 
    "P_Supporting"
  FSClass[score <= -0.18 & score > -0.36] = 
    "B_Supporting"
  FSClass[score <= -0.36] = 
    "B_Moderate"

  FSClass = factor(FSClass, levels = 
    c("P_Strong", "P_Moderate", "P_Supporting", "B_Supporting", "B_Moderate"))
  return(FSClass)
}
