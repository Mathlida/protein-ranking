% read weighted network
data1=load('PPI_hint.txt');
inta=data1(:,1);
intb=data1(:,2);
vals=data1(:,3);
ints=[inta intb];
genelist=unique(ints);
edge=length(inta);
vertex=length(genelist);
geneorder=1:vertex;
genemap=containers.Map(genelist,geneorder);
degree=zeros(1,vertex);
laba=zeros(1,2*edge);
labb=zeros(1,2*edge);
val=zeros(1,2*edge);
clear data1;
clear ints;
disp('network is read!');
% degree
for k=1:edge
    degree(genelist==inta(k))=degree(genelist==inta(k))+vals(k);
    degree(genelist==intb(k))=degree(genelist==intb(k))+vals(k);
end
for k=1:edge
    laba(2*k-1)=genemap(inta(k));
    labb(2*k-1)=genemap(intb(k));
    val(2*k-1)=vals(k);
    laba(2*k)=labb(2*k-1);
    labb(2*k)=laba(2*k-1);
    val(2*k)=vals(k);
end
% adjacency matrix
clear inta;
clear intb;
clear vals;
A=sparse(laba,labb,val,vertex,vertex);
disp('A is built!');
clear laba;
clear labb;
clear val;
D=diag(degree);
I=eye(vertex);
L=I-inv(D)*A;
clear D;
clear A;
disp('L is calculated!');
arfa=0.5;
N=3;
M=I-(arfa/N)*L;
clear L;
clear I;
K=M;
for i=1:(N-1)
    disp(i);
    K=K*M;
end
%read all candidata genes
candigenes=genelist;
nallgene=length(candigenes);
p0=zeros(1,nallgene);
clear data;
%read disease genes
data_dis=load('KDPs.txt');
diseasegenes=data_dis(:,1);
clear data_dis;
[samedis,dis]=intersect(candigenes,diseasegenes);
[difdis,allindex]=setdiff(candigenes,diseasegenes);
clear difdis;
p0(dis)=1;
index=dis';
clear samedis;
p_all=p0*K;
p_dis=[];
for i=index
    p0(i)=0;
    p=p0*K;
    p_dis=[p_dis p(i)];
    p0(i)=1;
end
p_all(dis)=p_dis;
p_all=p_all';
%write rank to rank_global_dis.txt
save rank_GR_KDPs.txt p_all -ascii;
