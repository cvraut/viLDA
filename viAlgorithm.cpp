#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <math.h>
#include <chrono>
#include <thread>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List svi(IntegerVector words, IntegerVector docIDs,
         IntegerVector topics, int lengthVocab,
         int numDocuments, int lengthDocuments, int maxIter,
         double alphaWords, double alphaTopics, double rho) {


  IntegerVector docStartIndices(numDocuments);
  for (int i = 0; i < numDocuments; i++) {
    docStartIndices[i] = i*lengthDocuments;
  }

  IntegerVector listDocIDs(numDocuments);
  for (int i = 0; i < numDocuments; i++) {
    listDocIDs[i] = i;
  }

  int numWords = words.length();
  int numTopics = topics.length();

  NumericVector alphaWordsVec(lengthVocab);
  NumericVector alphaTopicsVec(numTopics);

  for (int i = 0; i < lengthVocab; i++) {
    alphaWordsVec[i] = alphaWords;
  }

  for (int i = 0; i < numTopics; i++) {
    alphaTopicsVec[i] = alphaTopics;
  }

  NumericMatrix lambdaMat(numTopics, lengthVocab);
  for (int i = 0; i < numTopics; i++) {
    lambdaMat.row(i) = runif(lengthVocab);
    lambdaMat.row(i) = lambdaMat.row(i) / sum(lambdaMat.row(i));
  }

  NumericMatrix phiMat(numWords, numTopics);
  NumericMatrix gammaMat(numDocuments, numTopics);
  NumericMatrix intermediateTopics;
  NumericMatrix lambdaMatHat(numTopics, lengthVocab);

  NumericVector phiMatRowSum;
  NumericVector phiMatWeightedRowSum;
  NumericVector digammaSumgammaMat(numDocuments);
  NumericVector digammaSumlambdaMat(numTopics);

  IntegerVector singleSample;

  int sampledDoc;
  int docStartIndex;
  double phiMatEV;
  double phiMatNormalConst;

  NumericMatrix wordIndicatorMatrix(numWords, lengthVocab);
  for (int word = 0; word < numWords; word++) {
    for(int token = 0; token < lengthVocab; token++) {
      wordIndicatorMatrix(word, token) = 0;
    }
    wordIndicatorMatrix(word, words[word]) = 1;
  }

  Rcout << "starting SVI" << "\n";

  // Stochastic Variational Inference
  for (int iter = 0; iter < maxIter; iter++) {

    singleSample = sample(listDocIDs, 1, true);
    sampledDoc = singleSample[0];
    docStartIndex = docStartIndices[sampledDoc];

    // Initialize topic mixture variational parameters to one
    for (int doc = 0; doc < numDocuments; doc++) {
      for (int topic = 0; topic < numTopics; topic++) {
        gammaMat(doc, topic) = 1;
      }
    }

    // Calculating row sums of digamma function for topic mixture variational parameters
    for (int doc = 0; doc < numDocuments; doc++) {
      digammaSumgammaMat[doc] = 0;
      for (int topic = 0; topic < numTopics; topic++) {
        digammaSumgammaMat[doc] = digammaSumgammaMat[doc] + R::digamma(gammaMat(doc,topic));
      }
    }

    // calculating row sums of digamma function for word variational parameters
    for (int topic = 0; topic < numTopics; topic++) {
      digammaSumlambdaMat[topic] = 0;
      for (int token = 0; token < lengthVocab; token++) {
        digammaSumlambdaMat[topic] = digammaSumlambdaMat[topic] + R::digamma(lambdaMat(topic,token));
      }
    }

    // VI for parameters for a specific document
    for (int vbIter = 0; vbIter < maxIter; vbIter++) {

      // VI for topic distribution parameters for the sampled document
      // TODO: phi parameters are infinite because the arguments to the
      // exponential are too large. This is easy to see from the behavior
      // of the digamma function and the equations given in the SVI algorithm
      // for LDA
      for (int word = 0; word < lengthDocuments; word++) {
        for (int topic = 0; topic < numTopics; topic++) {
          phiMatEV = -1*digammaSumgammaMat[sampledDoc] - digammaSumlambdaMat[topic];
          phiMatEV = phiMatEV + R::digamma(gammaMat(sampledDoc,topic));
          phiMatEV = phiMatEV + R::digamma(lambdaMat(topic, words[docStartIndex + word]));
          phiMat(docStartIndex + word, topic) = exp(phiMatEV);
          Rcout << "current exp(phiMatEV) value: " << exp(phiMatEV) << "\n";
        }
        phiMatNormalConst = sum(phiMat.row(word));
        phiMat.row(word) = phiMat.row(word)/phiMatNormalConst;
      }

      // VI for topic mixture parameters for this document
      phiMatRowSum = phiMat.row(docStartIndex);
      for (int doc = 1; doc < lengthDocuments; doc++) {
        phiMatRowSum = phiMatRowSum + phiMat.row(docStartIndex + doc);
      }
      gammaMat.row(sampledDoc) = alphaTopicsVec + phiMatRowSum;
    }

    // VI for word distribution parameters (global)
    for (int topic = 0; topic < numTopics; topic++) {
      phiMatWeightedRowSum = phiMat(docStartIndex, topic) * wordIndicatorMatrix.row(docStartIndex);
      // Rcout << "started weighted row sum" << "\n";
      for (int word = 0; word < numWords; word++) {
        phiMatWeightedRowSum = phiMatWeightedRowSum + phiMat.row(word) * wordIndicatorMatrix.row(word);
      }
      lambdaMatHat.row(topic) = alphaWordsVec + numDocuments*phiMatWeightedRowSum;
    }

    // update lambda parameters
    for (int i = 0; i < numTopics; i++) {
      lambdaMat.row(i) = (1 - rho)*lambdaMat.row(i) + rho*lambdaMatHat.row(i);
    }

    Rcout << "gammaMat " << gammaMat << "\n";
    Rcout << "phiMat " << phiMat << "\n";
    Rcout << "lambdaMat " << lambdaMat << "\n";
    std::this_thread::sleep_for (std::chrono::seconds(8));


    if (iter % 1000 == 0) {
      rho = rho/2;
    }

  }

  List res(3);
  res[0] = lambdaMat;
  res[1] = phiMat;
  res[2] = gammaMat;
  return res;
}















