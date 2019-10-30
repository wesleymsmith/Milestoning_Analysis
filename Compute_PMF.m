msgbox("Locate the rate matrix file")
[matFileName,matFileDir]=uigetfile("./","Load Escape Matrix File")
a=load("-ascii",strcat(matFileDir,matFileName));
spy(a)  ; %check nonzero entries

for j=1:size(a);
a(j,j) = 1-sum(a(j,:)); %fill the diagonal 
end

[u v ] = eig(a'); %compute eigenvl/vct and sort 
[vs is ] = sort(diag(v),'descend'); 
us = u(:,is);

% invariant distribution
mu = us(:,1);
% invariant distribution normalized
mu=mu/sum(mu);

% free energy
kbt=inputdlg("Enter kbT"){1} %0.596;
kbt
fe=-kbt*log(mu)

subplot(2,1,1)
imagesc(a)
colorbar()

subplot(2,1,2)
plot(fe-fe(1), 'b')

%msgbox("Save equilibirium distribution to file")
%[piFileName,piFilePath]=uiputfile()
%save("-ascii",strcat(piFilePath,piFileName),"mu")