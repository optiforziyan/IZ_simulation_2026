suppressPackageStartupMessages(library(terra))
args<-commandArgs(); script<-sub('^--file=','',args[grepl('^--file=',args)][1])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/')
files<-list.files(file.path(root,'results/checkpoints'),full.names=TRUE,pattern='rds$')
report<-list()
for(f in files) {
 o<-readRDS(f); x<-o$v; a<-o$input
 stopifnot(all(x %in% c(0,1)),all(x[a$lock==1]==1),all(colSums(a$features*x)>=colSums(a$features)*a$targets-1e-6),o$log$gap<1e-8)
 if(!is.null(a$budget)) stopifnot(sum(x)<=a$budget+1e-6)
 if(!is.null(a$equal_area)) stopifnot(sum(x)==a$equal_area)
 objective<-if(a$objective=='minset') sum(a$cost*x) else 0
 if(!is.null(a$penalty)) objective<-objective+a$strength*sum(a$penalty*x)
 stopifnot(abs(objective-o$log$objective)<1e-5)
 report[[length(report)+1]]<-o$log
}
out<-do.call(rbind,report); write.csv(out,file.path(root,'results/solver_audit.csv'),row.names=FALSE)
cat('Verified',nrow(out),'cached optimizations: hard constraints, lock-in, budgets, objectives, zero solver gaps.\n')
