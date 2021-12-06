#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <math.h>
#include <chrono>
#include <thread>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List svi(IntegerMatrix data, IntegerVector numDistinctWordVec,
         IntegerVector topics, int lengthVocab, int numDocuments, int maxIterConst,
         int maxVBiterConst, double alphaWords, double alphaTopics, double rho, double tol) {


  IntegerVector docStartIndices(numDocuments);
  int cumulativeNumDocs = 0;
  for (int i = 0; i < numDocuments; i++) {
    docStartIndices[i] = cumulativeNumDocs;
    cumulativeNumDocs = cumulativeNumDocs + numDistinctWordVec[i];
  }

  IntegerVector listDocIDs(numDocuments);
  for (int i = 0; i < numDocuments; i++) {
    listDocIDs[i] = i;
  }

  int numWords = sum(data.column(2));
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

  NumericMatrix phiMat(data.nrow(), numTopics);
  NumericMatrix gammaMat(numDocuments, numTopics);
  NumericMatrix intermediateTopics;
  NumericMatrix lambdaMatHat(numTopics, lengthVocab);
  NumericMatrix oldLambdaMatrix(numTopics, lengthVocab);

  for (int topic = 0; topic < numTopics; topic++) {
    for (int token = 0; token < lengthVocab; token++) {
      oldLambdaMatrix[topic, token] = 0;
    }
  }

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

  NumericMatrix wordIndicatorMatrix(data.nrow(), lengthVocab);
  for (int word = 0; word < data.nrow(); word++) {
    for(int token = 0; token < lengthVocab; token++) {
      wordIndicatorMatrix(word, token) = 0;
    }
    wordIndicatorMatrix(word, data(word, 1)) = 1;
  }

  int maxEpoch = maxIterConst;
  int maxVBiter = maxVBiterConst;

  Rcout << "starting SVI" << "\n";

  // Stochastic Variational Inference
  for (int epoch = 0; epoch < maxEpoch; epoch++) {

    for (int iter = 0; iter < numDocuments; iter++) {
      // reset vector used to assess convergence in topic variational parameters
      NumericVector oldGammaVector(numTopics);
      for (int topic = 0; topic < numTopics; topic++) {
        oldGammaVector[topic] = 0;
      }

      // sample a document and update its local parameters in this iteration
      sampledDoc = iter;
      docStartIndex = docStartIndices[sampledDoc];

      // Initialize topic mixture variational parameters to one
      for (int topic = 0; topic < numTopics; topic++) {
        gammaMat(sampledDoc, topic) = 1;
      }

      // VI for parameters for a specific document
      for (int vbIter = 0; vbIter < maxVBiter; vbIter++) {

        // Calculating row sums of digamma function for topic mixture variational parameters
        digammaSumGammaMat = sum(digamma(gammaMat.row(sampledDoc)));

        // calculating row sums of digamma function for word variational parameters
        for (int topic = 0; topic < numTopics; topic++) {
          digammaSumlambdaMat[topic] = sum(digamma(lambdaMat.row(topic)));
        }

        // VI for topic distribution parameters for the sampled document
        for (int word = 0; word < numDistinctWordVec[sampledDoc]; word++) {

          // calculate EVs used in the local variational solution for topic probabilities
          for (int topic = 0; topic < numTopics; topic++) {
            phiMatArgument = -1*digammaSumGammaMat - digammaSumlambdaMat[topic];
            phiMatArgument = phiMatArgument + R::digamma(gammaMat(sampledDoc,topic));
            phiMatArgument = phiMatArgument + R::digamma(lambdaMat(topic, data(docStartIndex + word, 1)));
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
        phiMatRowSum = data(docStartIndex, 2) * phiMat.row(docStartIndex);
        for (int word = 1; word < numDistinctWordVec[sampledDoc]; word++) {
          phiMatRowSum = phiMatRowSum + data(docStartIndex + word, 2) * phiMat.row(docStartIndex + word);
        }
        gammaMat.row(sampledDoc) = alphaTopicsVec + phiMatRowSum;

      }

      // VI for word distribution parameters (global)
      for (int topic = 0; topic < numTopics; topic++) {
        phiMatWeightedRowSum = data(docStartIndex, 2) * phiMat(docStartIndex, topic) * wordIndicatorMatrix.row(docStartIndex);
        for (int word = 1; word < numDistinctWordVec[sampledDoc]; word++) {
          phiMatWeightedRowSum = phiMatWeightedRowSum + data(docStartIndex + word, 2) * phiMat(docStartIndex + word, topic) * wordIndicatorMatrix.row(docStartIndex + word);
        }
        lambdaMatHat.row(topic) = alphaWordsVec + numDocuments*phiMatWeightedRowSum;
      }

      // update lambda parameters
      for (int i = 0; i < numTopics; i++) {
        lambdaMat.row(i) = (1 - rho)*lambdaMat.row(i) + rho*lambdaMatHat.row(i);
      }
    }

  }

  // normalization
  for (int i = 0; i < numTopics; i++) {
    lambdaMat.row(i) = lambdaMat.row(i)/sum(lambdaMat.row(i));
  }

  for (int i = 0; i < numDocuments; i++) {
    gammaMat.row(i) = gammaMat.row(i)/sum(gammaMat.row(i));
  }



  List res(3);
  res[0] = lambdaMat;
  res[1] = phiMat;
  res[2] = gammaMat;
  return res;
}













