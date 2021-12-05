# example script to generate the data

source("R/data_gen.R")

n_docs = c(500,1000,2000)
n_vocabs = c(100,500,1000)
n_tops = c(2,5,10)

for(d in n_docs){
  for(v in n_vocabs){
    for(t in n_tops){
      print(sprintf("making data_file_d%d_v%d_t%d.Rda",d,v,t))
      data_obj <- data_gen(d,v,t)
      save(data_obj,file=sprintf("data_file_d%d_v%d_t%d.Rda",d,v,t))
      print(sprintf("wrote data_file_d%d_v%d_t%d.Rda",d,v,t))
    }
  }
}





# sourceCpp("viAlgorithm.cpp")
#
# tic()
# resList <- svi(words = words,
#                docIDs = docIDs,
#                lengthDocumentVec = lengthDocuments,
#                topics = topics,
#                lengthVocab = lengthVocab,
#                numDocuments = numDocuments,
#                maxIterConst = 100,
#                maxVBiterConst = 100,
#                alphaWords = alphaWords,
#                alphaTopics = alphaTopics,
#                rho = 0.1,
#                tol = 0.0001)
# toc()
# print("printing estimated variational parameters for the word distribution")
# print(resList[[1]])
# print("printing estimated variational parameters for the topic distribution")
# print(resList[[3]])
#
# print("printing true values")
# print(wordDistributions)
# print(generatedTopics)
# # print(resList)
