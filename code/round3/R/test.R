
test_data<-read.table("../../../data/originals/ch1_leaderBoard_monoTherapy.csv",sep=",",header=T)

submission_target_1a <- as.data.frame(cbind(test_data[,1],test_data[,14],as.data.frame(predicted_glm_target_a)))
colnames(submission_target_1a)=c("CELL_LINE","COMBINATION_ID","PREDICTION")
write.table(submission_target_1a,file="../../../data/round2/submission/submission_target_1a.txt",quote=F,sep="\t",col.names=T,row.names=F)

submission_target_1b <- as.data.frame(cbind(test_data[,1],test_data[,14],as.data.frame(predicted_glm_target_b)))
colnames(submission_target_1b)=c("CELL_LINE","COMBINATION_ID","PREDICTION")
write.table(submission_target_1b,file="../../../data/round2/submission/submission_target_1b.txt",quote=F,sep="\t",col.names=T,row.names=F)

submission_pathway_1a <- as.data.frame(cbind(test_data[,1],test_data[,14],as.data.frame(predicted_glm_pathway_a)))
colnames(submission_pathway_1a)=c("CELL_LINE","COMBINATION_ID","PREDICTION")
write.table(submission_pathway_1a,file="../../../data/round2/submission/submission_pathway_1a.txt",quote=F,sep="\t",col.names=T,row.names=F)

submission_pathway_1b <- as.data.frame(cbind(test_data[,1],test_data[,14],as.data.frame(predicted_glm_pathway_b)))
colnames(submission_pathway_1b)=c("CELL_LINE","COMBINATION_ID","PREDICTION")
write.table(submission_pathway_1b,file="../../../data/round2/submission/submission_pathway_1b.txt",quote=F,sep="\t",col.names=T,row.names=F)

submission_all_1a <- as.data.frame(cbind(test_data[,1],test_data[,14],as.data.frame(predicted_glm_all_a)))
colnames(submission_all_1a)=c("CELL_LINE","COMBINATION_ID","PREDICTION")
write.table(submission_all_1a,file="../../../data/round2/submission/submission_all_1a.txt",quote=F,sep="\t",col.names=T,row.names=F)

submission_all_1b <- as.data.frame(cbind(test_data[,1],test_data[,14],as.data.frame(predicted_glm_all_b)))
colnames(submission_all_1b)=c("CELL_LINE","COMBINATION_ID","PREDICTION")
write.table(submission_all_1b,file="../../../data/round2/submission/submission_all_1b.txt",quote=F,sep="\t",col.names=T,row.names=F)


