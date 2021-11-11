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

  // wordVec.print();
  // docIDvec.print();
  // topicVec.print();

  int numWords = wordVec.n_elem;
  int numTopics = topicVec.n_elem;

  // arma::ivec sampledTopics = RcppArmadillo::sample(topics, numWords, true, NumericVector::create());
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
      // Rcout << "topic" << k << "\n";
      word_ID_vec = arma::find(sampledTopics == k);
      // Rcout << "length words " << wordVec.n_elem << "\n";
      // Rcout << "found indices" << "\n";
      // word_ID_vec.print();
      wordsInTopic = wordVec.elem(word_ID_vec);
      // Rcout << "found words" << "\n";
      // wordsInTopic.print();
      topic_ID_vec = arma::find(wordsInTopic == w);
      // Rcout << "w" << w  << "\n";
      // Rcout << "about to access Cwt" << "\n";
      // Rcout << "length of topic ID vec" << topic_ID_vec.n_elem << "\n";
      Cwt(w,k) = topic_ID_vec.n_elem;
      // Cwt.print();
      // Rcout << "accessed Cwt" << "\n";
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

  // Cwt.print();

  // Rcout << "now printing Cdt" << "\n";

  // Cdt.print();

  // Rcout << "finished  initializing count matrices" << "\n";

  int numSamples = (numEpochs - warmUp)/lag - 1;

  // arma::Mat<double> wordParamSamples = arma::mat(numSamples, numWords*numTopics);
  // arma::Mat<double> topicParamSamples = arma::mat(numSamples, numDocuments*numTopics);

  List wordParamSamples(numSamples);
  List topicParamSamples(numSamples);

  // Rcout << "Instantiated lists" << "\n";

  arma::Row<double> probWord = arma::rowvec(numTopics);
  arma::Row<double> probDoc = arma::rowvec(numTopics);
  arma::Row<double> probTopic = arma::rowvec(numTopics);

  // Rcout << "instantiating probability vectors" << "\n";

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

      // Rcout << "calculating probability vectors" << "\n";

      for (int topic = 0; topic < numTopics; topic++) {
        probWord[topic] = (Cwt(currentWord,topic) + alphaWords) / (arma::sum(Cwt.col(topic)) + lengthVocab*alphaWords);
      }

      // Cdt.print();

      for (int topic = 0; topic < numTopics; topic++) {
        probDoc[topic] = (Cdt(currentDocID,topic) + alphaTopics) / (arma::sum(Cdt.row(currentDocID)) + numTopics*alphaTopics);
        // Rcout << "printing prob doc " << probDoc[topic] << "\n";
        // Rcout << "printing numerator " << Cdt(currentDocID,topic) + alphaTopics << "\n";
      }

      // Rcout << "trying to sample next topic" << "\n";

      probTopic = probWord % probDoc;

      // probWord.print();
      // probDoc.print();
      // probTopic.print();

      // arma::ivec sample = RcppArmadillo::sample(topicVec, 1, true, probTopic/arma::sum(probTopic));
      //
      // Rcout << "ran sample command" << "\n";
      //
      // sample.print();
      // sampledTopics[i] = sample(topics, 1, true, probTopic/arma::sum(probTopic))[0];

      IntegerVector singleSample = sample(topics, 1, true, wrap(probTopic/arma::sum(probTopic)));
      // Rcout << "ran sample command" << "\n";
      sampledTopics[i] = singleSample[0];

      // Rcout << "finished sampling next topic" << "\n";

      Cwt(currentWord, sampledTopics[i]) = Cwt(currentWord, sampledTopics[i]) + 1;
      Cdt(currentDocID, sampledTopics[i]) = Cdt(currentDocID, sampledTopics[i]) + 1;

    }

    if ((epoch > warmUp) && (epoch % lag == 0)) {

      // Rcout << "taking sample" << "\n";

      // Cdt.print();

      arma::Mat<double> wordParameters = arma::mat(lengthVocab, numTopics);
      arma::Mat<double> topicParameters = arma::mat(numDocuments, numTopics);
      for (int topic = 0; topic < numTopics; topic++) {
        for (int word = 0; word < lengthVocab; word++) {
          wordParameters(word,topic) = (Cwt(word,topic) + alphaWords) / (arma::sum(Cwt.col(topic)) + lengthVocab*alphaWords);
        }
      }
      for (int doc = 0; doc < numDocuments; doc++) {
        for (int topic = 0; topic < numTopics; topic++) {
          // Rcout << "print doc, topic " << doc << "," << topic << "\n";
          topicParameters(doc,topic) = (Cdt(doc,topic) + alphaTopics) / (arma::sum(Cdt.row(doc)) + numTopics*alphaTopics);
          // Rcout << "topicParam " << topicParameters[doc,topic] << "\n";
          // topicParameters.print();
        }
      }

      // Rcout << "topicParams " << "\n";
      // topicParameters.print();

      for (int doc = 0; doc < numDocuments; doc++) {
        topicParameters.row(doc) = topicParameters.row(doc)/arma::sum(topicParameters.row(doc));
      }

      // wordParamSamples.row(sampleIndex) = reshape(wordParameters, 1, numWords*numTopics);
      // topicParamSamples.row(sampleIndex) = reshape(topicParameters, 1, numDocuments*numTopics);
      // sampleIndex = sampleIndex + 1;

      // wordParameters.print();
      // topicParameters.print();

      // List res(1);
      // return res;

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

















