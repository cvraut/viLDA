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
         int numDocuments, int lengthDocuments, int maxIter, int maxVBiter,
         double alphaWords, double alphaTopics, double rho, double tol) {


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
  NumericVector digammaSumlambdaMat(numTopics);
  NumericVector phiMatArgumentArray(numTopics);

  IntegerVector singleSample;

  int sampledDoc;
  int docStartIndex;
  double phiMatArgument;
  double phiMatMaxArgument;
  double digammaSumGammaMat;

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

    // reset vector used to assess convergence in topic variational parameters
    NumericVector oldGammaVector(numTopics);
    for (int topic = 0; topic < numTopics; topic++) {
      oldGammaVector[topic] = 0;
    }

    // sample a document and update its local parameters in this iteration
    singleSample = sample(listDocIDs, 1, true);
    sampledDoc = singleSample[0];
    docStartIndex = docStartIndices[sampledDoc];

    // Initialize topic mixture variational parameters to one
    for (int topic = 0; topic < numTopics; topic++) {
      gammaMat(sampledDoc, topic) = 1;
    }

    // VI for parameters for a specific document
    for (int vbIter = 0; vbIter < maxVBiter; vbIter++) {

      // Calculating row sums of digamma function for topic mixture variational parameters
      digammaSumGammaMat = 0;
      for (int topic = 0; topic < numTopics; topic++) {
        digammaSumGammaMat = digammaSumGammaMat + R::digamma(gammaMat(sampledDoc,topic));
      }

      // calculating row sums of digamma function for word variational parameters
      for (int topic = 0; topic < numTopics; topic++) {
        digammaSumlambdaMat[topic] = 0;
        for (int token = 0; token < lengthVocab; token++) {
          digammaSumlambdaMat[topic] = digammaSumlambdaMat[topic] + R::digamma(lambdaMat(topic,token));
        }
      }

      // VI for topic distribution parameters for the sampled document
      for (int word = 0; word < lengthDocuments; word++) {

        // calculate EVs used in the local variational solution for topic probabilities
        for (int topic = 0; topic < numTopics; topic++) {
          phiMatArgument = -1*digammaSumGammaMat - digammaSumlambdaMat[topic];
          phiMatArgument = phiMatArgument + R::digamma(gammaMat(sampledDoc,topic));
          phiMatArgument = phiMatArgument + R::digamma(lambdaMat(topic, words[docStartIndex + word]));
          phiMatArgumentArray[topic] = phiMatArgument;
        }

        // subtract max exponential argument from all values to avoid numerical errors
        phiMatMaxArgument = max(phiMatArgumentArray);
        for (int topic = 0; topic < numTopics; topic++) {
          phiMat(docStartIndex + word, topic) = exp(phiMatArgumentArray[topic] - phiMatMaxArgument);
        }

        // normalize topic probabilities to one
        phiMat.row(docStartIndex + word) = phiMat.row(docStartIndex + word)/sum(phiMat.row(docStartIndex + word));
      }

      // VI for topic mixture parameters for this document
      phiMatRowSum = phiMat.row(docStartIndex);
      for (int word = 1; word < lengthDocuments; word++) {
        phiMatRowSum = phiMatRowSum + phiMat.row(docStartIndex + word);
      }
      gammaMat.row(sampledDoc) = alphaTopicsVec + phiMatRowSum;

      // break if convergence is reached
      if (sum(pow(oldGammaVector - gammaMat.row(sampledDoc), 2)) < tol) {
        break;
      } else {
        oldGammaVector = gammaMat.row(sampledDoc);
      }

    }

    // VI for word distribution parameters (global)
    for (int topic = 0; topic < numTopics; topic++) {
      phiMatWeightedRowSum = phiMat(docStartIndex, topic) * wordIndicatorMatrix.row(docStartIndex);
      for (int word = 1; word < lengthDocuments; word++) {
        phiMatWeightedRowSum = phiMatWeightedRowSum + phiMat(docStartIndex + word, topic) * wordIndicatorMatrix.row(docStartIndex + word);
      }
      lambdaMatHat.row(topic) = alphaWordsVec + numDocuments*phiMatWeightedRowSum;
    }

    // update lambda parameters
    for (int i = 0; i < numTopics; i++) {
      lambdaMat.row(i) = (1 - rho)*lambdaMat.row(i) + rho*lambdaMatHat.row(i);
    }

    // decrease learn rate parameter
    if (iter % 1000 == 0) {
      rho = rho/2;
    }

  }

  // Predict the topic assigned to each document
  NumericVector predictedTopics(numDocuments);
  NumericVector sumTopicProbs(lengthDocuments);
  for (int doc = 0; doc < numDocuments; doc++) {
    sumTopicProbs = phiMat.row(doc*lengthDocuments);
    for (int word = 1; word < lengthDocuments; word++) {
      sumTopicProbs = sumTopicProbs + phiMat.row(doc*lengthDocuments + word);
    }
    predictedTopics[doc] = which_max(sumTopicProbs);
  }

  List res(4);
  res[0] = lambdaMat;
  res[1] = phiMat;
  res[2] = gammaMat;
  res[3] = predictedTopics;
  return res;
}















