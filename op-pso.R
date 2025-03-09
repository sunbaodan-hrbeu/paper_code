
#directed acyclic detection
legal2<-function(A,nVar){
  C<-A
  for (i in 1:nVar){
    if (sum(diag(C))==0){
      C<-A%*%C
      if ((i==nVar)&&(sum(diag(C))==0)){
        return(A)
      }
    }
    else{
      location=which(C!=0,arr.ind = TRUE)
      a<-sample(1:length(location[,"row"]),1)
      row_random<-as.numeric(location[a,"row"])
      for (j in 1:nVar){
        if ((!is.na(A[row_random,j]))&&(A[row_random,j]==1)){
          if (row_random==j){
            A[row_random,j]=0
          }
          else{
            A[row_random,j]=0
            A[j,row_random]=1
          }
        }
      }
      A<-legal2(A,nVar)
      C<-A
    }
  }
}
#crossover withparticle best solution
Crossover1<-function(k,popc,particle_Best,c1,nVar){
  newpop1<-c()
  initial_pop1<-popc$position[[k]]
  initial_pop2<-particle_Best$position[[k]]
  initial_pop1<-t(initial_pop1)
  initial_pop2<-t(initial_pop2)
  initial_pop1<-as.vector(initial_pop1)
  initial_pop2<-as.vector(initial_pop2)
  for (i in 1:length(initial_pop1)) {
    if((runif(1,0,1)<c1)||(i==sample(1:(nVar*nVar),1))){
      newpop1[i]<-initial_pop1[i]
    }else{
      newpop1[i]<-initial_pop2[i]
    }
  }
  newpop1<-matrix(newpop1,nrow=nVar,ncol=nVar)
  newpop1<-t(newpop1)
  return(newpop1)
}
#crossover with global best solution
Crossover2<-function(m,popc,GlobalBest,c2,nVar){
  newpop2<-c()
  initial_pop3<-popc$position[[m]]
  initial_pop4<-GlobalBest$position
  initial_pop3<-t(initial_pop3)
  initial_pop4<-t(initial_pop4)
  initial_pop3<-as.vector(initial_pop3)
  initial_pop4<-as.vector(initial_pop4)
  for (i in 1:length(initial_pop3)) {
    if((runif(1,0,1)<=c2)||(i==sample(1:(nVar*nVar),1))){
      newpop2[i]<-initial_pop3[i]
    }else{
      newpop2[i]<-initial_pop4[i]
    }
  }
  newpop2<-matrix(newpop2,nrow=nVar,ncol=nVar)
  newpop2<-t(newpop2)
  return(newpop2)
}
#mutation operator
Mutation<-function(j,pop,w,nVar,nPop){
  ran_i<-c()
  F_scale<-0.5
  xigma<-0.001
  v<-c()
  ran_i<-sample(1:nPop,3,replace = FALSE)
  while((ran_i[1]==j)||(ran_i[1]==ran_i[2])||(ran_i[2]==ran_i[3])||(ran_i[1]==ran_i[3])||(ran_i[2]==j)||(ran_i[3]==j)){
    ran_i<-sample(1:nPop,3,replace = FALSE)
    if((ran_i[1]!=j)&&(ran_i[1]!=ran_i[2])&&(ran_i[2]!=ran_i[3])&&(ran_i[1]!=ran_i[3])&&(ran_i[2]!=j)&&(ran_i[3]!=j)){
      break
    }
  }
  Xbest<-pop$position[[ran_i[1]]]
  Xbest_score<-pop$score_dag[ran_i[1]]
  Xr1<-pop$position[[ran_i[2]]]
  Xr2<-pop$position[[ran_i[3]]]
  if((pop$score_dag[ran_i[2]]>Xbest_score)&&(pop$score_dag[ran_i[2]]>pop$score_dag[ran_i[3]])){
    Xbest<-pop$position[[ran_i[2]]]
    Xbest_score<-pop$score_dag[ran_i[2]]
    Xr1<-pop$position[[ran_i[1]]]
    Xr2<-pop$position[[ran_i[3]]]
  }else{
    if((pop$score_dag[ran_i[3]]>Xbest_score)&&(pop$score_dag[ran_i[3]]>pop$score_dag[ran_i[2]])){
      Xbest<-pop$position[[ran_i[3]]]
      Xbest_score<-pop$score_dag[ran_i[3]]
      Xr1<-pop$position[[ran_i[1]]]
      Xr2<-pop$position[[ran_i[2]]]
    }
  }
  Xr1<-t(Xr1)
  Xr2<-t(Xr2)
  Xbest<-t(Xbest)
  Xr1<-as.vector(Xr1)
  Xr2<-as.vector(Xr2)
  Xbest<-as.vector(Xbest)
  if(Xbest_score<pop$score_dag[j]){
    Ci<-pmin(c(rep(1,(nVar*nVar))),(F_scale*xor(Xr1,Xr2)+xigma))
  }else{
    Ci<-c(rep(xigma,(nVar*nVar)))
  }
  for (kb in 1:(nVar*nVar)) {
    if(Ci[kb]<w){
      v[kb]<-Xbest[kb]
    }else{
      v[kb]<-1-Xbest[kb]
    }
  }
  v<-matrix(v,nrow=nVar,ncol=nVar)
  v<-t(v)
  return(v)
}
#opposite_learning
#cc is the interation time
#eli_pop is the elite solution number
opposite_learning<-function(cc,eli_pop,nVar){ 
    M<-eli_pop$position[[cc]]
    M<-t(M)
    node_num<-c(1:nVar)
    position_node<-c()
    for (op in 1:nVar) {
      if((all(M[op,]==0))&&(all(M[,op]==0))){
        position_node<-append(position_node,op)
      }
    }
    if(length(position_node)!=0){
      node_number<-sample(position_node,1)#determine the node to add one edge
      nrow_number<-sample(node_num[-node_number],1)
      M[nrow_number,node_number]<-1
    }
    M<-legal2(M,nVar)
    return(M)
}

sparse_dag_in<-function(dag){
  arc_list <- arcs(dag)
  in_degree <- data.frame(table(arc_list[, 2]))
  # 提取网络中的边
  nodes_high_in_degree <- in_degree[in_degree$Freq>11,]
  nodes_high_in_degree<-as.vector(nodes_high_in_degree$Var1)
  while(length(nodes_high_in_degree)>0){
    for (node in nodes_high_in_degree) {
      # 找到该节点的父节点（即入边）
      parents_of_node <- parents(dag, node)
      if (length(parents_of_node) > 0) {
        parent_to_remove <- sample(parents_of_node, 1)
        dag <- drop.arc(dag, from = parent_to_remove, to = node)
      }
    }
    arc_list <- arcs(dag)
    in_degree <- data.frame(table(arc_list[, 2]))
    nodes_high_in_degree <- in_degree[in_degree$Freq>11,]
    nodes_high_in_degree<-as.vector(nodes_high_in_degree$Var1)
  }
  return(dag)
}
sparse_dag_out<-function(dag){
  arc_list <- arcs(dag)
  out_degree <- data.frame(table(arc_list[, 1]))
  # 提取网络中的边
  nodes_high_out_degree <- out_degree[out_degree$Freq>11,]
  nodes_high_out_degree<-as.vector(nodes_high_out_degree$Var1)
  while (length(nodes_high_out_degree)>0) {
    for (node2 in nodes_high_out_degree) {
      # 找到该节点的父节点（即入边）
      children_of_node <- children(dag, node2)
      if (length(children_of_node) > 0) {
        child_to_remove <- sample(children_of_node, 1)
        dag <- drop.arc(dag, from =node2, to = child_to_remove)
      }
    }
    arc_list <- arcs(dag)
    out_degree <- data.frame(table(arc_list[, 1]))
    nodes_high_out_degree <- out_degree[out_degree$Freq>11,]
    nodes_high_out_degree<-as.vector(nodes_high_out_degree$Var1)
  }
  return(dag)
}

#main program
nVar=27
MaxIt=100
nPop=50
eli_num<-25
nrow=1
ncol=50
w_start=0.9
w_end=0.35
c1_start=0.84
c1_end=0.52
c2_start=0.38
c2_end=0.81
succ<-0
Cx<-c()
popc<-list()
BestScore<-c()
order_score<-c()
sort_score<-c()

#initialization
row_name<-vector(mode="numeric")
col_name<-vector(mode="numeric")
#initial population
initial_Position<-list()
initial_Score<-matrix(nrow=nrow,ncol=ncol)
initial_pop<-list(position=initial_Position,score_dag=initial_Score)
#population
Position<-list()
Score<-matrix(nrow=nrow,ncol=ncol)
pop<-list(position=Position,score_dag=Score)
#elite_population
eli_Position<-list()
eli_Score<-matrix(nrow=nrow,ncol=eli_num)
eli_pop<-list(position=eli_Position,score_dag=eli_Score)
#opposite population
op_Position<-list()
op_Score<-matrix(nrow=nrow,ncol=eli_num)
op_pop<-list(position=op_Position,score_dag=op_Score)
#particle initialization
particle_Position<-list()
particle_Score<-matrix(nrow=nrow,ncol=ncol)
particle_Velocity<-list()
particle<-list(position=particle_Position,score_dag=particle_Score,velocity=particle_Velocity)
#particle best initialization
particle_Best_Position<-list()
particle_Best_Score<-matrix(nrow=nrow,ncol=ncol)
particle_Best<-list(position=particle_Best_Position,score_dag=particle_Best_Score)
#global best initialization
GlobalBest_Score<--100000000
GlobalBest_Position<-matrix(nrow=nVar,ncol=nVar)
GlobalBest<-list(GlobalBest_Score,GlobalBest_Position)
BestScore<-matrix(0,nrow=MaxIt,ncol = 1)

# load dataset
load('/Users/sunbaodan/Desktop/experiment network/insurance.rda')
new_cancer=rbn(bn,1000)
saveRDS(new_cancer,"/Users/sunbaodan/Desktop/sample_data/insurance_1000.rda")
new_loan<-readRDS("/Users/sunbaodan/Desktop/sample_data/insurance_pc_1000.rda")
#new_loan=data.frame(new_cancer)
#names(new_loan)<-c(1:nVar)
#structure priors(pc algorithm)
m=pc.stable(new_loan)
new_m=amat(m)
a=1
for (i in 1:nVar){
  for(j in 1:nVar){
    if ((new_m[i,j]==1)&&(new_m[j,i]==0)){
      row_name[a]=i
      col_name[a]=j
      a=a+1
    }
  }
}

#generate the first population
load('/Users/sunbaodan/Desktop/experiment network/insurance.rda')
dag_true<-amat(bn)
#saveRDS(new_cancer,"C:/Users/student/Documents/pc10/cancer_500.rda")
#new_loandata<-read.csv("C:/Users/student/Documents/cancer_1000.csv")
new_loandata<-readRDS("/Users/sunbaodan/Desktop/sample_data/insurance_1000.rda")
name_varaiable<-names(new_loandata)
e <-empty.graph(name_varaiable)
amat(e)<-dag_true
true_score<-score(e,new_loandata,type="bic")

#initial solution
for (k in 1:nPop){
  e<-random.graph(name_varaiable)
  dag<-amat(e)
  for (p in 1:length(row_name)){
    dag[row_name[p],col_name[p]]<-1
  }
  dag<-legal2(dag,nVar)
  Position[[k]]<-dag
  dg0<-empty.graph(name_varaiable)
  amat(dg0)<-dag
  Score[nrow,k]<-score(dg0,new_loandata,type="bic")#scores
}
initial_pop<-list(position=Position,score_dag=Score)
order_score<-order(initial_pop$score_dag,decreasing = TRUE)
sort_score<-sort(initial_pop$score_dag,decreasing = TRUE)
initial_pop$position<-initial_pop$position[order_score]

#obtain opposition-based learning population
for (i in 1:eli_num) {
  eli_Position[[i]]<-initial_pop$position[[i]]
  eli_Score[nrow,i]<-sort_score[i]
}
eli_pop<-list(position=eli_Position,score_dag=eli_Score)

#initial elite opposite solutions
for (elitn in 1:eli_num) {
  op_Position[[elitn]]<- opposite_learning(elitn,eli_pop,nVar)
  initial_pop$position[[nPop+elitn]]<-op_Position[[elitn]]#merge two lists position
  dg_eli<-empty.graph(name_varaiable)
  amat(dg_eli)<-op_Position[[elitn]]
  op_Score[nrow,elitn]<-score(dg_eli,new_loandata,type="bic")#scores
  initial_pop$score_dag[nPop+elitn]<-op_Score[nrow,elitn]#merge two lists scores
}
op_pop<-list(position=op_Position,score_dag=op_Score)
order_score<-order(initial_pop$score_dag,decreasing = TRUE)
sort_score<-sort(initial_pop$score_dag,decreasing=TRUE)
initial_pop$position<-initial_pop$position[order_score]
initial_pop<-list(position=initial_pop$position,score_dag=sort_score)
pop$position<-initial_pop$position[1:nPop]
pop$score_dag<-initial_pop$score_dag[1:nPop]


#initialize best solutions
for (k in 1:nPop){
  particle_Position[[k]]<-pop$position[[k]]
  particle_Score[nrow,k]<-pop$score_dag[k]
  #update the particle best solution
  particle_Best_Position[[k]]<-particle_Position[[k]]
  particle_Best_Score[nrow,k]<-particle_Score[nrow,k]
  if (particle_Best_Score[nrow,k]>=GlobalBest_Score){
    GlobalBest_Position<-particle_Best_Position[[k]]
    GlobalBest_Score<-particle_Best_Score[nrow,k]
  }
}
particle<-list(position=particle_Position,score_dag=particle_Score,velocity=particle_Velocity)
particle_Best<-list(position=particle_Best_Position,score_dag=particle_Best_Score)
GlobalBest<-list(position=GlobalBest_Position,score_dag=GlobalBest_Score)
BestSol<-list(position=pop$position[[1]],score_dag=sort_score[1])
GlobalBest<-BestSol
BestScore<-matrix(0,nrow=MaxIt,ncol = 1)

#iteration process
for (it in 1:MaxIt){
  Position2<-list()
  Score2<-matrix(nrow=nrow,ncol=ncol)
  popc<-list(position=Position2,score_dag=Score2)
  w=(w_start-w_end)*(1-it/MaxIt)+w_end*(succ/nPop)
  c1=c1_start-(c1_start-c1_end)*it/MaxIt
  c2=c2_start-(c2_start-c2_end)*it/MaxIt
  succ<-0.00
  #mutate
  for (j in 1:nPop){
    popc$position[[j]]<-Mutation(j,pop,w,nVar,nPop)
    particle$velocity[[j]]<-popc$position[[j]]
  }
  #crossover
  for (k2 in 1:nPop) {
    popc$position[[k2]]<-Crossover1(k2,popc,particle_Best,c1,nVar)
  }
  #crossover
  for (mc in 1:nPop) {
    popc$position[[mc]]<-Crossover2(mc,popc,GlobalBest,c2,nVar)
    E<-popc$position[[mc]]
    popc$position[[mc]]<-legal2(E,nVar)
    E<-popc$position[[mc]]
    dg3<-empty.graph(name_varaiable)
    amat(dg3)<-E
    dg3<-sparse_dag_in(dg3)
    dg3<-sparse_dag_out(dg3)
    popc$position[[mc]]<-amat(dg3)
    popc$score_dag[mc]<-score(dg3,new_loandata,type="bic")
    particle$position[[mc]]<-popc$position[[mc]]
  }
  #merge pop and popc
  Position2<-c(pop$position,popc$position)
  Score2<-c(pop$score_dag,popc$score_dag)
  pop<-list(position=Position2,score_dag=Score2)
  order_score<-order(pop$score_dag,decreasing = TRUE)
  sort_score<-sort(pop$score_dag,decreasing = TRUE)
  pop$position<-pop$position[order_score]
  pop<-list(position=pop$position,score_dag=sort_score)
  #repeat generate elite opoosite solutions
  #obtain elite population
  for (ii in 1:eli_num) {
    eli_Position[[ii]]<-pop$position[[ii]]
    eli_Score[nrow,ii]<-sort_score[ii]
  }
  eli_pop<-list(position=eli_Position,score_dag=eli_Score)
  #initial elite opposite solutions
  for (elitn in 1:eli_num) {
    op_Position[[elitn]]<- opposite_learning(elitn,eli_pop,nVar)
    dg_eli<-empty.graph(name_varaiable)
    amat(dg_eli)<-op_Position[[elitn]]
    dg_eli<-sparse_dag_in(dg_eli)
    dg_eli<-sparse_dag_out(dg_eli)
    op_Position[[elitn]]<-amat(dg_eli)
    pop$position[[nPop+elitn]]<-op_Position[[elitn]]#merge two lists position
    op_Score[nrow,elitn]<-score(dg_eli,new_loandata,type="bic")#scores
    pop$score_dag[nPop+elitn]<-op_Score[nrow,elitn]#merge two lists scores
  }
  op_pop<-list(position=op_Position,score_dag=op_Score)
  order_score<-order(pop$score_dag,decreasing = TRUE)
  sort_score<-sort(pop$score_dag,decreasing=TRUE)
  pop$position<-pop$position[order_score]
  pop<-list(position=pop$position,score_dag=sort_score)
  pop$position<-pop$position[1:nPop]
  pop$score_dag<-pop$score_dag[1:nPop]
  pop1<-list()
  for (r in 1:nPop){
    pop1$position[[r]]<-pop$position[[r]]
    pop1$score_dag[r]<-pop$score_dag[r]
  }
  pop<-pop1
  BestSol<-list(position=pop$position[[1]],score_dag=pop$score_dag[1])
  #update parameters of particles
  for (n in 1:nPop) {
    if (pop$score_dag[[n]]>particle$score_dag[n]){
      succ<-succ+1
      particle$position[[n]]<-pop$position[[n]]
      particle$score_dag[n]<-pop$score_dag[n]
    }
    Cx[n]=particle$score_dag[n]
  }
  r1<-which.max(particle$score_dag)
  BestScore[it]<-max(Cx)
  GlobalBest$score_dag<-particle$score_dag[r1]
  GlobalBest$position<-particle$position[[r1]]
  ##update parameters of particles
  for (num in 1:nPop){
    G<-particle$position[[num]]
    #G<-legal2(G,nVar)
    dg4=empty.graph(name_varaiable)
    amat(dg4)<-G
    particle$score_dag[num]=score(dg4,new_loandata,type="bic")
    #update the best solutions of particles
    if (particle$score_dag[num]>particle_Best$score_dag[num]){
      particle_Best$position[[num]]=particle$position[[num]]
      particle_Best$score_dag[num]=particle$score_dag[num]
      if (particle_Best$score_dag[num]>GlobalBest$score_dag){
        GlobalBest$position=particle_Best$position[[num]]
        GlobalBest$score_dag=particle_Best$score_dag[num]
      }
    }
  }
  for (num2 in 1:nPop) {
    if (particle$score_dag[num2]>=pop$score_dag[num2]){
      pop$position[[num2]]=particle$position[[num2]]
      pop$score_dag[num2]=particle$score_dag[num2]
    }
    Cx[num2]=pop$score_dag[num2]
  }
  r2<-which.max(pop$score_dag)
  BestScore[it]=max(Cx)
  GlobalBest$score_dag=pop$score_dag[r2]
  GlobalBest$position=pop$position[[r2]]
}