#' @export
getTP = function(estimate, truth) {
    estimate[estimate!=0] = 1
    truth[truth!=0] = 1
    TP = estimate * truth
    sum(TP)
}

#' @export
getAdjTP = function(estimate, truth) {
    estimate[estimate!=0] = 1
    truth[truth!=0] = 1
    truth[t(truth)!=0] = 1
    TP = estimate * truth
    sum(TP)
}

#' @export
getFP = function(estimate, truth) {
    estimate[estimate!=0] = 1
    tmp = truth
    tmp[truth!=0] = 0
    tmp[truth==0] = 1
    FP = estimate * tmp
    sum(FP)
}

#' @export
getAdjFP = function(estimate, truth) {
    estimate[estimate!=0] = 1
    truth[truth!=0] = 1
    truth[t(truth)!=0] = 1
    tmp = truth
    tmp[truth!=0] = 0
    tmp[truth==0] = 1
    FP = estimate * tmp
    sum(FP)
}

#' @export
getFN = function(estimate, truth) {
    truth[truth!=0] = 1
    tmp = estimate
    tmp[estimate!=0] = 0
    tmp[estimate==0] = 1
    FN = truth * tmp
    sum(FN)
}

#' @export
getAdjFN = function(estimate, truth) {
    truth[truth!=0] = 1
    estimate[estimate!=0] = 1
    estimate[t(estimate)!=0] = 1
    tmp = estimate
    tmp[estimate!=0] = 0
    tmp[estimate==0] = 1
    FN = truth * tmp
    sum(FN)
}

#' @export
getTN = function(estimate, truth) {
    tmp1 = estimate
    tmp2 = truth

    tmp1[estimate!=0] = 0
    tmp1[estimate==0] = 1

    tmp2[truth!=0] = 0
    tmp2[truth==0] = 1

    TN = tmp1 * tmp2
    sum(TN)
}

#' @export
getPrecision = function(estimate, truth) {
    TP=getTP(estimate, truth)
    FP=getFP(estimate, truth)
    FN=getFN(estimate, truth)
    precision=TP*1.0/(TP+FP)
    precision*100
}

#' @export
getRecall = function(estimate, truth) {
    TP=getTP(estimate, truth)
    FP=getFP(estimate, truth)
    FN=getFN(estimate, truth)
    recall=TP*1.0/(TP+FN)
    recall*100
}

#' @export
getF1 = function(estimate, truth) {
    recall = getRecall(estimate, truth)
    precision = getPrecision(estimate, truth)
    F1 = 2 * precision * recall / (precision + recall)
    F1
}

#' @export
getAdjPrecision = function(estimate, truth) {
    TP=getAdjTP(estimate, truth)
    FP=getAdjFP(estimate, truth)
    FN=getAdjFN(estimate, truth)
    precision=TP*1.0/(TP+FP)
    precision*100
}

#' @export
getAdjRecall = function(estimate, truth) {
    TP=getAdjTP(estimate, truth)
    FP=getAdjFP(estimate, truth)
    FN=getAdjFN(estimate, truth)
    recall=TP*1.0/(TP+FN)
    recall
}

#' @export
getAdjF1 = function(estimate, truth) {
  recall = getAdjRecall(estimate, truth)
  precision = getAdjPrecision(estimate, truth)
  F1 = 2 * precision * recall / (precision + recall)
  F1*100
}

#' @export
getOrientAccuracy = function(estimate, truth) {
    OrientTP = getTP(estimate, truth)
    AdjTP = getAdjTP(estimate, truth)
    accuracy = 0
    if (AdjTP != 0) {
        accuracy=OrientTP/AdjTP
    }
    accuracy*100
}

#' @export
myShd = function (m1, m2) {
    m1 = as.matrix(m1)
    m2 = as.matrix(m2)
    m1[m1 != 0] <- 1
    m2[m2 != 0] <- 1
    shd <- 0
    s1 <- m1 + t(m1)
    s2 <- m2 + t(m2)
    s1[s1 == 2] <- 1
    s2[s2 == 2] <- 1
    ds <- s1 - s2
    ind <- which(ds > 0)
    m1[ind] <- 0
    shd <- shd + length(ind)/2
    ind <- which(ds < 0)
    m1[ind] <- m2[ind]
    shd <- shd + length(ind)/2
    d <- abs(m1 - m2)
    shd + sum((d + t(d)) > 0)/2
}