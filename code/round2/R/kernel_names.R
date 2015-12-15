#!/usr/bin/env Rscript
mutations_angular_similarity <- c('../../data/round2/kernels/angular_similarity_mutations_driver',
                                  '../../data/round2/kernels/angular_similarity_mutations_nonsilent',
                                  '../../data/round2/kernels/angular_similarity_mutations_original')

mutations_dot_product <- c('../../data/round2/kernels/dot_product_mutations_driver',
                                  '../../data/round2/kernels/dot_product_mutations_nonsilent',
                                  '../../data/round2/kernels/dot_product_mutations_original')

genex_angular_similarity <- c('../../data/round2/kernels/angular_similarity_genex_original',
                                  '../../data/round2/kernels/angular_similarity_genex_filtered',NA)

genex_dot_product <- c('../../data/round2/kernels/dot_product_genex_original',
                           '../../data/round2/kernels/dot_product_genex_filtered',NA)

cnv_angular_similarity <- c('../../data/round2/kernels/angular_similarity_cnv',NA,NA)

cnv_dot_product <- c('../../data/round2/kernels/dot_product_cnv',NA,NA)
methylation_angular_similarity <- c('../../data/round2/kernels/angular_similarity_methylation_shores_filtered',
                                  '../../data/round2/kernels/angular_similarity_methylation_shores_original',
                                  '../../data/round2/kernels/angular_similarity_methylation_islands')

methylation_dot_product <- c('../../data/round2/kernels/dot_product_methylation_shores_filtered',
                             '../../data/round2/kernels/dot_product_methylation_shores_original',
                             '../../data/round2/kernels/dot_product_methylation_islands')

drug_angular_similarity <- c('../../data/round2/kernels/angular_similarity_drug_pathway',
                              '../../data/round2/kernels/angular_similarity_drug_target',NA)

drug_dot_product <- c('../../data/round2/kernels/dot_product_drug_pathway',
                       '../../data/round2/kernels/dot_product_drug_target',NA)

kernel_names <- data.frame(mutations_angular_similarity,mutations_dot_product,genex_angular_similarity
                           ,genex_dot_product,cnv_angular_similarity,cnv_dot_product,methylation_angular_similarity
                           ,methylation_dot_product,drug_angular_similarity,drug_dot_product)
names(kernel_names) <- c("mutations_angular_similarity","mutations_dot_product","genex_angular_similarity"
                         ,"genex_dot_product","cnv_angular_similarity","cnv_dot_product","methylation_angular_similarity"
                         ,"methylation_dot_product","drug_angular_similarity","drug_dot_product")

write.table(kernel_names,"../../data/round2/kernels/kernel_names.txt",sep="\t",row.names = F,col.names = T,quote = F)
