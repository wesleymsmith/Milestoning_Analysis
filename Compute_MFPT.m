
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

%ncells = number of cells (milestone windows
%nmiles = number of edges between cells (ncells-1)
ncells=length(mu);
nmiles=ncells-1;

% compute Nij from Nij^a

%since this is 1D the milestone (edge) between
%cells i and j is equal to min(i,j)
%Nij tracks transitions from edge i to edge j
%these are found in the entry for cell min(i,j)
%in entry i of allN, column 1 is for edge i to i+1 (left edge to right edge)
% and column 2 is for edge i+1 to i (right edge to left edge)
N=zeros(nmiles,nmiles);
for k=1:nmiles;
N(k,k+1)=mu(k+1)*allN(k+1,1)/allT(k+1);
if k<nmiles;
  N(k+1,k)=mu(k+1)*allN(k+1,2)/allT(k+1);
endif
endfor
N
% compute Ri from Ri^a
% We use the same convention as for Nij above
% in allR, then entry for cell i contains the
% timecount for edge i in the first column and
% the timecount for edge i+1 in the second column
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
%Qij = Nij / Ri
Q=zeros(nmiles,nmiles);
for k=1:nmiles-1;
Q(k,k+1)=N(k,k+1)/R(k);
Q(k+1,k)=N(k+1,k)/R(k+1);
end

%a=Q;
%aa=a;
%na=size(a,1);

%% Reduce Q by eliminating zero lines and rows, when needed

%ic = 0;
%jkept=[];
%for j=1:na;
% if (sum((a(j,:))==zeros(1,na))) == na
%  aa (j-ic,:) = [];
%  ic = ic+1;
% else
%     jkept = [jkept j];   %% this array will keep track of the kept milestones in the old numeration
% end
%end

%ic = 0;
%for j=1:na;
% if (sum((a(:,j))==zeros(na,1))) == na
%  aa (:,j-ic) = [];
%  ic = ic+1;
% end
%end

%b=aa;
%na=size(a)
%naa=size(aa)

%% Rebuild diagonal terms, and eliminate cemetery state 
aa=Q

 for j=1:size(aa);
   aa(j,j) = -sum(aa(j,:)); 
 end

% b=aa;
% aa(:,naa) = [];
% aa(naa,:) = [];
 
 nr=size(aa,1)

Qij=aa
 
%% Compute mfpts 

 % take into account timestep & frames to have actual time
   % frq of writing frames
   %     nfrq = 100;
   nfrq=inputdlg("Enter number of timesteps between frames"){1}
   % timestep 2 fs 
   %  dt = nfrq*2*1e-15; %(2 fs in s)
   dts=str2double(inputdlg("Enter simulated timestep duration (in seconds)"){1})
   dt=nfrq*dts
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