% script qall.m to build global Q and compute mfpt

clear all
#load prob.dat
#mu=prob
mu=load("-ascii","cell.reverse.pi.dat")
#load allN
allN=load("-ascii","cell.reverse.Nij.dat")
#load allR
allR=load("-ascii","cell.reverse.Rij.dat")
#load allT
allT=load("-ascii","cell.reverse.dat")

%ncells = number of cells - the former and the last cell
ncells=length(mu);
nmiles=ncells-1;

% compute Nij from Nij^a

N=zeros(nmiles,nmiles);
for k=1:ncells;
N(k,k+1)=mu(k)*allN(k,1)/allT(k);
N(k+1,k)=mu(k)*allN(k,2)/allT(k);
end
N
% compute Ri from Ri^a

#R(1)=mu(1)*allR(1,1)/allT(1);
#R(nmiles)=mu(ncells)*allR(ncells,2)/allT(ncells);
#for j=2:nmiles-1;
#R(j)=mu(j-1)*allR(j-1,2)/allT(j-1)+mu(j)*allR(j,1)/allT(j);
#end

R=zeros(nmiles,1);
%j accumulate for right side cell edges
%  *the right hand edge for cell j is edge j
%k accumulate for left side cell edges
%  *the left hand edge for cell k is edge k-1
%Note j counts up from 1 to nmiles (number of edges)
% while k counts down from ncells to 2
for j=1:nmiles;
  k=ncells-j+1;
  R(j)=R(j)+mu(j)*allR(j,2)/allT(j);
  R(k-1)=R(k-1)+mu(k)*allR(k,1)/allT(k);
end
R

% compute global Q

Q=zeros(nmiles,nmiles);
for k=1:nmiles-1;
#disp(k)
#disp(N(k,k+1)/R(k))
#disp(N(k+1,k)/R(k+1))
Q(k,k+1)=N(k,k+1)/R(k);
Q(k+1,k)=N(k+1,k)/R(k+1);
#disp("----")
end
Q

a=Q;
aa=a;
na=size(a,1);

%% Reduce Q by eliminating zero lines and rows, when needed

ic = 0;
jkept=[];
for j=1:na;
 if (sum((a(j,:))==zeros(1,na))) == na
  aa (j-ic,:) = [];
  ic = ic+1;
 else
     jkept = [jkept j];   %% this array will keep track of the kept milestones in the old numeration
 end
end

ic = 0;
for j=1:na;
 if (sum((a(:,j))==zeros(na,1))) == na
  aa (:,j-ic) = [];
  ic = ic+1;
 end
end

b=aa;
na=size(a)
naa=size(aa)

%% Rebuild diagonal terms, and eliminate cemetery state

 for j=1:size(aa);
   aa(j,j) = -sum(aa(j,:));
 end

 b=aa;
 aa(:,naa) = [];
 aa(naa,:) = [];
 
 nr=size(aa,1)
    
%% Compute mfpts

 % take into account timestep & frames to have actual time
   % frq of writing frames
        nfrq = 500;
   % timestep 2 fs
     dts=.002 / 1000 #ts from ps to ns
     dt = nfrq*dts; %(2 time per snapshot in ns)
 
  % solving linear system equation (5) Maragliano et al. JCTC 2009
t1=aa\(-ones(nr,1));

tt1= sort(t1,'descend');
figure;
hold on; plot(tt1*dt,'b')
tt1*dt

% plot also MFPT progressive along the milestone index (i.e. along the
% reaction)
format long

t=tt1*dt;
t2=t(1,1)-t;
figure;
plot(t2)
dlmwrite('t2.reverse.dat',t2,'delimiter',' ','precision', '%.6e')