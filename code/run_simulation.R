# Reproducible synthetic SCP experiment. No empirical SDM or genetic inference.
# Run: Rscript simulation_v2/code/run_simulation.R [main|replicates|all]
suppressPackageStartupMessages({library(terra); library(prioritizr)})
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[1] else 'all'
script <- sub('^--file=', '', commandArgs()[grepl('^--file=',commandArgs())][1])
root <- normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
for(d in c('results','figures','data')) dir.create(file.path(root,d),showWarnings=FALSE)
RNGkind('Mersenne-Twister','Inversion','Rejection')
norm01 <- function(x) if(diff(range(x))==0) x*0 else (x-min(x))/diff(range(x))
as_layer <- function(template,v,nm='value') {r<-template; values(r)<-v; names(r)<-nm; r}
csv <- function(x,name) write.csv(x,file.path(root,'results',paste0(name,'.csv')),row.names=FALSE)

generate_landscape <- function(seed, export=FALSE) {
 set.seed(seed)
 r<-rast(nrows=24,ncols=24,xmin=0,xmax=24000,ymin=0,ymax=24000,crs='EPSG:3857')
 xy<-crds(r)/24000; x<-xy[,1]; y<-xy[,2]; n<-ncell(r)
 # Three broad, three intermediate, two endemic species; smooth connected ranges.
 pars<-data.frame(species=paste0('Sp',1:8),class=c(rep('widespread',3),rep('intermediate',3),rep('endemic',2)),
 cx=c(.40,.63,.49,.24,.75,.56,.13,.87),cy=c(.52,.60,.34,.73,.31,.80,.16,.88),
 sx=c(.30,.27,.29,.14,.13,.12,.055,.06),sy=c(.25,.28,.20,.12,.15,.10,.07,.055))
 S<-H<-matrix(0,n,8)
 for(j in 1:8) {
  cx<-pars$cx[j]+runif(1,-.025,.025); cy<-pars$cy[j]+runif(1,-.025,.025)
  a<-runif(1,-.5,.5); dx<-x-cx; dy<-y-cy
  u<-cos(a)*dx+sin(a)*dy; v<- -sin(a)*dx+cos(a)*dy
  smooth<- .10*sin(2*pi*x+runif(1,0,2*pi))*cos(2*pi*y+runif(1,0,2*pi))
  H[,j]<-exp(-.5*((u/pars$sx[j])^2+(v/pars$sy[j])^2)+smooth)
  H[,j]<-H[,j]/max(H[,j]); S[,j]<-as.numeric(H[,j]>=.30)
 }
 colnames(S)<-colnames(H)<-pars$species
 ranges<-colSums(S); frac<-ranges/n
 alpha<-ifelse(frac>=.25,.20,ifelse(frac>=.10,.30,.50))
 # Toy rooted tree: 2 four-species clades, 4 sister pairs, 8 terminal branches.
 B<-rbind(c(1,1,1,1,0,0,0,0),c(0,0,0,0,1,1,1,1),
          c(1,1,0,0,0,0,0,0),c(0,0,1,1,0,0,0,0),c(0,0,0,0,1,1,0,0),c(0,0,0,0,0,0,1,1),diag(8))
 bl<-c(2,2,1,1,1,1,rep(1,8)); branch<-(S%*%t(B))>0
 # Four discrete, overlapping functional entities; FD = number locally present.
 T<-rbind(c(1,0,0,0),c(1,0,1,0),c(0,1,0,0),c(0,1,1,0),
          c(0,0,1,0),c(0,1,0,1),c(1,0,0,1),c(0,0,0,1))
 F<-(S%*%T)>0; colnames(F)<-paste0('Function',1:4)
 # Species-specific spatial population strata, not measured genetic lineages.
 P<-matrix(0,n,0); popmeta<-data.frame()
 for(j in 1:8) {
  ids<-which(S[,j]>0); np<-if(j<=3) 3 else 2
  angle<-runif(1,0,pi); gradient<-x[ids]*cos(angle)+y[ids]*sin(angle)
  ord<-ids[order(gradient)]; groups<-cut(seq_along(ord),breaks=np,labels=FALSE)
  for(k in 1:np) {
   z<-numeric(n); z[ord[groups==k]]<-1; P<-cbind(P,z)
   popmeta<-rbind(popmeta,data.frame(population=paste0('Sp',j,'_Pop',k),species=paste0('Sp',j),pi=runif(1,.002,.020),range_cells=sum(z)))
  }
 }
 colnames(P)<-popmeta$population
 D<-cbind(TD=rowSums(S),PD=as.numeric(branch%*%bl),FD=rowSums(F),GD=as.numeric(P%*%popmeta$pi))
 Z<-apply(D[,c('PD','FD','GD')],2,norm01)
 # Two compact pre-existing reserves shared by all scenarios within each seed.
 lock<-as.numeric(((x-.43)^2+(y-.48)^2<.085^2)|((x-.69)^2+(y-.60)^2<.055^2))
 pars$range_cells<-ranges; pars$range_fraction<-frac; pars$target_fraction<-alpha
 pars$target_cells<-ceiling(ranges*alpha-1e-9)
 landscape<-list(seed=seed,r=r,S=S,H=H,D=D,Z=Z,P=P,F=F,branch=branch,bl=bl,lock=lock,
                  alpha=alpha,pars=pars,popmeta=popmeta)
 if(export) {
  dr<-file.path(root,'data'); dir.create(file.path(dr,'SDM'),showWarnings=FALSE)
  for(j in 1:8) {
   writeRaster(as_layer(r,S[,j]),file.path(dr,'SDM',paste0('Sp',j,'_presence.tif')),overwrite=TRUE)
   writeRaster(as_layer(r,S[,j]),file.path(dr,'SDM',paste0('Sp',j,'_SDM.tif')),overwrite=TRUE)
   writeRaster(as_layer(r,H[,j]),file.path(dr,'SDM',paste0('Sp',j,'_suitability.tif')),overwrite=TRUE)
  }
  for(d in colnames(D)) writeRaster(as_layer(r,D[,d]),file.path(dr,paste0(d,'.tif')),overwrite=TRUE)
  writeRaster(as_layer(r,lock),file.path(dr,'Locked_in.tif'),overwrite=TRUE)
  write.csv(pars,file.path(dr,'species_parameters.csv'),row.names=FALSE)
  write.csv(popmeta,file.path(dr,'population_parameters.csv'),row.names=FALSE)
  write.csv(T,file.path(dr,'species_functional_entities.csv'),row.names=FALSE)
  write.csv(data.frame(branch=paste0('b',1:nrow(B)),length=bl,B),file.path(dr,'branch_species_incidence.csv'),row.names=FALSE)
  writeLines('(((Sp1:1,Sp2:1):1,(Sp3:1,Sp4:1):1):2,((Sp5:1,Sp6:1):1,(Sp7:1,Sp8:1):1):2);',file.path(dr,'toy_tree.tre'))
  write.csv(S,file.path(dr,'species_presence_matrix.csv'),row.names=FALSE)
  write.csv(P,file.path(dr,'population_presence_matrix.csv'),row.names=FALSE)
  write.csv(F*1,file.path(dr,'functional_presence_matrix.csv'),row.names=FALSE)
  write.csv(branch*1,file.path(dr,'branch_presence_matrix.csv'),row.names=FALSE)
  portable<-landscape; portable$r<-wrap(as_layer(r,rep(0,n)))
  saveRDS(portable,file.path(dr,'landscape.rds'))
 }
 landscape
}

solver_log<-data.frame()
optimize <- function(L,label,features=L$S,targets=L$alpha,cost=rep(1,nrow(L$S)),
                     penalty=NULL,strength=1,budget=NULL,equal_area=NULL,objective='minset',extra=NULL) {
 gc(FALSE)
 cache_dir<-file.path(root,'results/checkpoints'); dir.create(cache_dir,showWarnings=FALSE)
 cache<-file.path(cache_dir,paste0(L$seed,'_',label,'.rds'))
 fingerprint<-list(features=features,targets=targets,cost=cost,penalty=penalty,strength=strength,budget=budget,equal_area=equal_area,objective=objective,extra=extra,lock=L$lock)
 if(file.exists(cache)) {
  old<-readRDS(cache)
  if(identical(old$input,fingerprint)) {solver_log<<-rbind(solver_log,old$log);return(old$v)}
 }
 fr<-rast(lapply(seq_len(ncol(features)),function(j) as_layer(L$r,features[,j],paste0('f',j))))
 p<-problem(as_layer(L$r,cost),fr) |> add_relative_targets(targets)
 if(objective=='minset') p<-add_min_set_objective(p) else p<-add_min_penalties_objective(p,budget=budget)
 if(!is.null(penalty)) p<-add_linear_penalties(p,strength,as_layer(L$r,penalty))
 if(!is.null(budget)&&objective=='minset') p<-add_linear_constraints(p,budget,'<=',rep(1,nrow(features)))
 if(!is.null(equal_area)) p<-add_linear_constraints(p,equal_area,'=',rep(1,nrow(features)))
 if(!is.null(extra)) p<-add_linear_constraints(p,extra$threshold,extra$sense,extra$data)
 p<-p |> add_locked_in_constraints(as_layer(L$r,L$lock)) |> add_binary_decisions() |>
  add_highs_solver(gap=0,time_limit=120,threads=1,verbose=FALSE,control=list(random_seed=1))
 t<-proc.time()[3]; sol<-solve(p); elapsed<-proc.time()[3]-t
 v<-as.numeric(values(sol)); if(anyNA(v)) stop('No solution: ',label)
 stopifnot(all(abs(v-round(v))<1e-6)); v<-round(v)
 held<-colSums(features*v); required<-colSums(features)*targets
 stopifnot(all(held>=required-1e-6),all(v[L$lock==1]==1))
 if(!is.null(budget)) stopifnot(sum(v)<=budget+1e-6)
 if(!is.null(equal_area)) stopifnot(sum(v)==equal_area)
 # prioritizr attaches exact solver status and gaps to solution output.
 info<-attributes(sol)
 stopifnot(length(info$gap)==1,is.finite(info$gap),info$gap<1e-8)
 logrow<-data.frame(seed=L$seed,scenario=label,elapsed_seconds=elapsed,area=sum(v),status=info$status,gap=info$gap,objective=info$objective)
 solver_log<<-rbind(solver_log,logrow)
 saveRDS(list(input=fingerprint,v=v,log=logrow),cache)
 cat(L$seed,label,'area',sum(v),'seconds',round(elapsed,2),'\n'); flush.console()
 v
}
evaluate <- function(L,v,label,weights=c(1/3,1/3,1/3)) {
 spheld<-colSums(L$S*v); popheld<-colSums(L$P*v)
 sums<-colSums(L$D*v); q<-as.numeric(L$Z%*%weights)
 data.frame(seed=L$seed,scenario=label,area=sum(v),area_percent=100*mean(v),
  species_targets_met=sum(spheld>=colSums(L$S)*L$alpha-1e-6),species_present=sum(spheld>0),
  min_species_target_ratio=min(spheld/(colSums(L$S)*L$alpha)),
  local_PD_sum=sums['PD'],local_FD_sum=sums['FD'],local_GD_sum=sums['GD'],Q_sum=sum(v*q),
  network_PD=sum(L$bl[colSums(L$branch*v)>0]),network_FD=sum(colSums(L$F*v)>0),
  populations_present=sum(popheld>0),population_targets_met=sum(popheld>=.2*colSums(L$P)-1e-6),
  total_populations=ncol(L$P),row.names=NULL)
}
export_solution <- function(L,v,label) writeRaster(as_layer(L$r,v,label),file.path(root,'results',paste0(label,'.tif')),overwrite=TRUE)

main <- function() {
 L<-generate_landscape(20260925,TRUE); q<-rowMeans(L$Z); penalty<-1-q
 sols<-list(); metrics<-data.frame()
 add<-function(v,label) {sols[[label]]<<-v; metrics<<-rbind(metrics,evaluate(L,v,label)); export_solution(L,v,label)}
 base<-optimize(L,'Species_min_area'); K<-sum(base); add(base,'Species_min_area')
 # Envelope over ALL species-feasible minimum-area alternatives: avoids lucky baseline.
 lo<-optimize(L,'Species_Q_min_at_Cstar',penalty=q,equal_area=K,objective='penalties')
 hi<-optimize(L,'Species_Q_max_at_Cstar',penalty=-q,equal_area=K,objective='penalties')
 add(lo,'Species_Q_min_at_Cstar'); add(hi,'Species_Q_max_at_Cstar')
 # True aggregated-feature conservation targets: explicit indices, valid but different requirements.
 av<-optimize(L,'Aggregated_features_20pct',features=L$D,targets=rep(.2,4)); add(av,'Aggregated_features_20pct')
 # Fixed-area diagnostic isolates loss of identity; no species constraints.
 hot<-optimize(L,'Aggregated_score_fixed_area',features=L$D,targets=rep(0,4),penalty=-q,equal_area=K,objective='penalties')
 add(hot,'Aggregated_score_fixed_area')
 for(rho in c(.25,1,4)) {
  a<-optimize(L,paste0('Effective_cost_rho',rho),cost=1+rho*penalty)
  b<-optimize(L,paste0('Soft_penalty_rho',rho),penalty=penalty,strength=rho)
  oa<-sum((1+rho*penalty)*a); ob<-sum((1+rho*penalty)*b)
  stopifnot(abs(oa-ob)<1e-5)
  add(a,paste0('Effective_cost_rho',rho)); add(b,paste0('Soft_penalty_rho',rho))
 }
 sensitivity<-data.frame()
 wsets<-list(equal=c(1,1,1)/3,PD_focus=c(.6,.2,.2),FD_focus=c(.2,.6,.2),GD_focus=c(.2,.2,.6),previous=c(.3,.3,.4))
 for(wn in names(wsets)) for(eps in c(0,.05,.10,.20)) {
  w<-wsets[[wn]]; qw<-as.numeric(L$Z%*%w); cap<-floor(K*(1+eps)+1e-8)
  nm<-paste0('Hier_',wn,'_eps',eps)
  v<-optimize(L,nm,penalty=-qw,budget=cap,objective='penalties')
  m<-evaluate(L,v,nm,w); m$weight_set<-wn; m$epsilon<-eps; m$area_cap<-cap
  sensitivity<-rbind(sensitivity,m)
  if(wn=='equal') add(v,nm)
 }
 # Additional explicit population representation illustrates what extra identity constraints do.
 feat<-cbind(L$S,L$P); targ<-c(L$alpha,rep(.2,ncol(L$P)))
 pop<-optimize(L,'Species_and_populations_min_area',features=feat,targets=targ)
 add(pop,'Species_and_populations_min_area')
 kp<-sum(pop)
 popq<-optimize(L,'Species_and_populations_Q',features=feat,targets=targ,penalty=-q,budget=kp,objective='penalties')
 add(popq,'Species_and_populations_Q')
 # Policy sensitivity, explicitly illustrative and independent of observed outcomes.
 targetcheck<-data.frame()
 for(pol in c('range_based','uniform_30pct','range_plus_10pct')) {
  L2<-L; L2$alpha<-switch(pol,range_based=L$alpha,uniform_30pct=rep(.3,8),range_plus_10pct=pmin(L$alpha+.1,1))
  bv<-optimize(L2,paste0('Target_',pol)); kk<-sum(bv)
  hv<-optimize(L2,paste0('TargetHier_',pol),penalty=-q,budget=floor(kk*1.1),objective='penalties')
  mm<-rbind(evaluate(L2,bv,'min_area'),evaluate(L2,hv,'hier_epsilon_0.1')); mm$target_policy<-pol
  targetcheck<-rbind(targetcheck,mm)
 }
 csv(metrics,'main_comparison'); csv(sensitivity,'weight_epsilon_sensitivity'); csv(targetcheck,'target_sensitivity')
 details<-do.call(rbind,lapply(names(sols),function(nm) data.frame(scenario=nm,species=L$pars$species,
    range_cells=L$pars$range_cells,target_fraction=L$alpha,required_cells=L$pars$target_cells,
    held_cells=colSums(L$S*sols[[nm]]),held_fraction=colSums(L$S*sols[[nm]])/colSums(L$S))))
 csv(details,'species_representation'); saveRDS(sols,file.path(root,'results/solutions.rds'))
 saveRDS(list(Cstar=K,metrics=metrics,sensitivity=sensitivity),file.path(root,'results/main_analysis.rds'))
 plot_main(L,sols,sensitivity)
}

plot_main <- function(L,sols,sens) {
 draw<-function(v,title,col) {plot(as_layer(L$r,v),main=title,col=col,axes=FALSE,legend=FALSE,mar=c(1,1,3,1))}
 png(file.path(root,'figures/01_species_ranges.png'),width=2400,height=1300,res=180)
 par(mfrow=c(2,4),oma=c(0,0,2,0))
 for(j in 1:8) draw(L$S[,j],sprintf('%s  %s\nRange %.1f%% | target %.0f%%',L$pars$species[j],L$pars$class[j],100*L$pars$range_fraction[j],100*L$alpha[j]),c('#f0f3f5','#187e89'))
 mtext('Synthetic clustered ranges; 24 x 24 equal-area planning units',outer=TRUE,cex=1.2); dev.off()
 png(file.path(root,'figures/02_local_indices.png'),width=2000,height=1700,res=180)
 par(mfrow=c(2,2))
 for(d in colnames(L$D)) plot(as_layer(L$r,L$D[,d]),main=paste('Local',if(d=='GD') 'GD proxy' else d),col=hcl.colors(50,'YlGnBu',rev=TRUE),axes=FALSE)
 dev.off()
 keys<-c('Species_min_area','Aggregated_features_20pct','Aggregated_score_fixed_area','Species_Q_max_at_Cstar','Hier_equal_eps0.1','Species_and_populations_Q')
 png(file.path(root,'figures/03_network_comparison.png'),width=2400,height=1600,res=180)
 par(mfrow=c(2,3),oma=c(2,0,0,0))
 titles<-c('Species minimum area','Aggregated features: 20% targets','Aggregated score: same area','Species + best score at C*','Species + score: 10% cost tolerance','Species + population hard targets')
 for(i in seq_along(keys)) {v<-sols[[keys[i]]]; draw(v+L$lock,paste0(titles[i],'\n',sum(v),' cells'),c('#f0f3f5','#2195a6','#192d45'))}
 mtext('Light grey: unselected   |   Teal: newly selected   |   Navy: existing reserves (locked in)',side=1,outer=TRUE,cex=.9)
 dev.off()
 png(file.path(root,'figures/04_tradeoffs.png'),width=1900,height=900,res=180)
 par(mfrow=c(1,2),mar=c(5,5,3,1)); ss<-sens[sens$weight_set=='equal',]
 plot(ss$area,ss$Q_sum,type='b',pch=19,col='#187e89',xlab='Selected area (cells)',ylab='Sum of normalized local index score',main='Local-score gain with cost tolerance')
 plot(ss$area,ss$population_targets_met,type='b',pch=19,col='#b96b34',ylim=c(0,ncol(L$P)),xlab='Selected area (cells)',ylab='Population targets met (20% each)',main='Population representation is not guaranteed')
 dev.off()
}

replicates <- function() {
 out<-data.frame()
 # Predeclared seed sequence; no screening or exclusion of unfavourable landscapes.
 for(seed in 20260925:20260934) {
  L<-generate_landscape(seed); q<-rowMeans(L$Z)
  b<-optimize(L,'rep_species_min'); K<-sum(b)
  h0<-optimize(L,'rep_hier_eps0',penalty=-q,budget=K,objective='penalties')
  h1<-optimize(L,'rep_hier_eps0.1',penalty=-q,budget=floor(K*1.1),objective='penalties')
  ag<-optimize(L,'rep_aggregate_same_area',features=L$D,targets=rep(0,4),penalty=-q,equal_area=K,objective='penalties')
  for(nm in c('b','h0','h1','ag')) {m<-evaluate(L,get(nm),nm);m$Cstar<-K;out<-rbind(out,m)}
  csv(out,'replicate_results')
 }
}
if(mode %in% c('main','all')) main()
if(mode %in% c('replicates','all')) replicates()
csv(solver_log,'solver_log')
writeLines(c(capture.output(sessionInfo()),paste('HiGHS R package:',packageVersion('highs')),
 paste('RNG:',paste(RNGkind(),collapse=', '))),file.path(root,'results/sessionInfo.txt'))
cat('COMPLETED',mode,'\n')
