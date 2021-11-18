#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <chrono>
#include <thread>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List gibbsSampler(IntegerVector words, IntegerVector docIDs,
                  IntegerVector topics, int lengthVocab, int numDocuments,
                  double alphaWords, double alphaTopics,
                  int numEpochs, int warmUp, int lag) {




  arma::ivec wordVec = as<arma::ivec>(words);
  arma::ivec docIDvec = as<arma::ivec>(docIDs);
  arma::ivec topicVec = as<arma::ivec>(topics);

  int numWords = wordVec.n_elem;
  int numTopics = topicVec.n_elem;

  IntegerVector sampledTopicsIV = sample(topics, numWords, true);
  arma::ivec sampledTopics = as<arma::ivec>(sampledTopicsIV);
  arma::Mat<int> Cwt = arma::imat(lengthVocab, numTopics);
  arma::Mat<int> Cdt = arma::imat(numDocuments, numTopics);

  arma::ivec wordsInTopic;
  arma::ivec topicsInDoc;

  arma::uvec word_ID_vec;
  arma::uvec topic_ID_vec;
  arma::uvec doc_ID_vec;

  Rcout << "initializing count matrices" << "\n";

  for (int w = 0; w < lengthVocab; w++) {
    for (int k = 0; k < numTopics; k++) {
      word_ID_vec = arma::find(sampledTopics == k);
      wordsInTopic = wordVec.elem(word_ID_vec);
      topic_ID_vec = arma::find(wordsInTopic == w);
      Cwt(w,k) = topic_ID_vec.n_elem;
    }
  }


  for (int d = 0; d < numDocuments; d++) {
    for (int k = 0; k < numTopics; k++) {
      doc_ID_vec = arma::find(docIDvec == d);
      topicsInDoc = sampledTopics.elem(doc_ID_vec);
      // Rcout << "accessed sampled topics" << "\n";
      topic_ID_vec = arma::find(topicsInDoc == k);
      Cdt(d,k) = topic_ID_vec.n_elem;
    }
  }

  int numSamples = (numEpochs - warmUp)/lag - 1;
  List wordParamSamples(numSamples);
  List topicParamSamples(numSamples);

  arma::Row<double> probWord = arma::rowvec(numTopics);
  arma::Row<double> probDoc = arma::rowvec(numTopics);
  arma::Row<double> probTopic = arma::rowvec(numTopics);

  int currentWord;
  int currentDocID;
  int sampleIndex = 0;

  Rcout << "starting Gibbs sampler" << "\n";

  // Gibbs Sampler
  for (int epoch = 0; epoch < numEpochs; epoch++) {
    for (int i = 0; i < numWords; i++) {
      currentWord = words[i];
      currentDocID = docIDs[i];
      Cwt(currentWord, sampledTopics[i]) = Cwt(currentWord, sampledTopics[i]) - 1;
      Cdt(currentDocID, sampledTopics[i]) = Cdt(currentDocID, sampledTopics[i]) - 1;

      for (int topic = 0; topic < numTopics; topic++) {
        probWord[topic] = (Cwt(currentWord,topic) + alphaWords) / (arma::sum(Cwt.col(topic)) + lengthVocab*alphaWords);
      }

      for (int topic = 0; topic < numTopics; topic++) {
        probDoc[topic] = (Cdt(currentDocID,topic) + alphaTopics) / (arma::sum(Cdt.row(currentDocID)) + numTopics*alphaTopics);
      }

      probTopic = probWord % probDoc;
      IntegerVector singleSample = sample(topics, 1, true, wrap(probTopic/arma::sum(probTopic)));
      sampledTopics[i] = singleSample[0];

      Cwt(currentWord, sampledTopics[i]) = Cwt(currentWord, sampledTopics[i]) + 1;
      Cdt(currentDocID, sampledTopics[i]) = Cdt(currentDocID, sampledTopics[i]) + 1;

    }

    if ((epoch > warmUp) && (epoch % lag == 0)) {

      arma::Mat<double> wordParameters = arma::mat(lengthVocab, numTopics);
      arma::Mat<double> topicParameters = arma::mat(numDocuments, numTopics);
      for (int topic = 0; topic < numTopics; topic++) {
        for (int word = 0; word < lengthVocab; word++) {
          wordParameters(word,topic) = (Cwt(word,topic) + alphaWords) / (arma::sum(Cwt.col(topic)) + lengthVocab*alphaWords);
        }
      }
      for (int doc = 0; doc < numDocuments; doc++) {
        for (int topic = 0; topic < numTopics; topic++) {
          topicParameters(doc,topic) = (Cdt(doc,topic) + alphaTopics) / (arma::sum(Cdt.row(doc)) + numTopics*alphaTopics);
        }
      }

      for (int doc = 0; doc < numDocuments; doc++) {
        topicParameters.row(doc) = topicParameters.row(doc)/arma::sum(topicParameters.row(doc));
      }

      wordParamSamples[sampleIndex] = wrap(wordParameters);
      topicParamSamples[sampleIndex] = wrap(topicParameters);
      sampleIndex = sampleIndex + 1;
    }

  }

  List res(2);
  res[0] = wordParamSamples;
  res[1] = topicParamSamples;

  return res;

}

















