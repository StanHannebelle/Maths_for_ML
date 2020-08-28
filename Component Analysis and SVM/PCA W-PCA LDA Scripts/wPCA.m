function U_reduct = wPCA(X,dim)

%Compared to the notation of lectures X corresponds to X^T here.
%So we transpose X to correspond to the notation of the lectures.
X = X';

% Mean computation
%mu = mean(X,'all').*ones(size(X,1),size(X,2));
mu = mean(X,1);
%Dot product matrix computation
DPM = (X-mu)'*(X-mu);
% we observe that norm(DPM - DPM') is not null ( close to 10^-8) whereas
% norm(0.5*(DPM+DPM') - (0.5*(DPM+DPM'))') is null
% so the following transformation turn DPM into a symetrical matrix
DPM = 0.5.*(DPM+DPM');


%Eigen Analysis
[V,Diag] = eig(DPM);

%Compute Eigen Vectors
Diag_minushalf = diag(1./sqrt(diag(Diag)));
%We add a multiplication by Diag_minushalf to have a whitening PCA
U = (X-mu)*V*Diag_minushalf*Diag_minushalf;

%Keeping specific number of first components
U = fliplr(U);
U_reduct = U(:,1:dim);
end