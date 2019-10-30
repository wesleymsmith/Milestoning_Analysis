
% script qall.m to build global Q and compute mfpt

clear all

msgbox("Load Equilibrium Window Distribution Data From PMF step (should end in .pi.dat)")
[piFileName,piFileDir]=uigetfile()
mu=load("-ascii",strcat(piFileDir,piFileName))

msgbox("Load Nij_alpha count data (should end in .Nij.dat)")
[nijFileName,nijFileDir]=uigetfile()
allN=load("-ascii",strcat(nijFileDir,nijFileName))

msgbox("Load Rij_alpha count data (should end in .Rij.dat)")
[rijFileName,rijFileDir]=uigetfile()
allR=load("-ascii",strcat(rijFileDir,rijFileName))

msgbox("Load window occupancy count data (should end in .T.dat)")
[tFileName,tFileDir]=uigetfile()
allT=load("-ascii",strcat(tFileDir,tFileName))

%ncells = number of cells - the former and the last cell
%this should be one less than the total number of entries
%in the equilibrium cell distribution data file
ncells=length(mu)-1;
nmiles=ncells+1;

% compute Nij from Nij^a

N=zeros(nmiles,nmiles);
for k=1:ncells;
N(k,k+1)=mu(k)*allN(k,1)/allT(k);
N(k+1,k)=mu(k)*allN(k,2)/allT(k);
end
N
% compute Ri from Ri^a

R(1)=mu(1)*allR(1,1)/allT(1);
R(nmiles)=mu(ncells)*allR(ncells,2)/allT(ncells);
for j=2:nmiles-1;
R(j)=mu(j-1)*allR(j-1,2)/allT(j-1)+mu(j)*allR(j,1)/allT(j);
end
R
% compute global Q

Q=zeros(nmiles,nmiles);
for k=1:nmiles-1;
Q(k,k+1)=N(k,k+1)/R(k);
Q(k+1,k)=N(k+1,k)/R(k+1);
end

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

Qij=aa
 
%% Compute mfpts 

 % take into account timestep & frames to have actual time
   % frq of writing frames
        nfrq = 100;
   % timestep 2 fs 
     dt = nfrq*2*1e-15; %(2 fs in s)
 
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
dlmwrite('t2.dat',t2,'delimiter',' ','precision', '%.6e')