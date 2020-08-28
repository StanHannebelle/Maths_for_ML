function U_reduc = LDA(X,classes)
    
%Compared to the notation of lectures X corresponds to X^T here.
%So we transpose X to correspond to the notation of the lectures.
X = X';

n = size(X,2);
dim = n- (length(unique(classes))+1);

%Computation of M Block Diagonal Matrix
[bincounts,ind] = histc(classes,unique(classes));

e_matrix_array = cell(1,length(unique(classes)));
for i=1:length(unique(classes))
    e_matrix_array{i} = (1./bincounts(i)).*ones(bincounts(i));
end
M = blkdiag(e_matrix_array{i},e_matrix_array{i});
for i = 3:length(unique(classes))
    M = blkdiag(M,e_matrix_array{i});
end

%Computation of the Sw Matrix
Sw = (eye(size(M))-M)*(X')*X*(eye(size(M))-M);

% As we saw in PCA.m and wPCA.m we need to do the following transformation
% to ensure that Sw is a symetrical matrix.
Sw = 0.5.*(Sw+Sw');

%Then, we can perform an eigen analysis of Sw
[Vw,Diagw] = eig(Sw);


%Then we need to remove the null eigen values and corresponding vectors
%To do that we put the eigen vectors corresponding to the null eigen values at the bottom of Vw
[sorted_eigval,rearrangement_index] = sort(diag(Diagw),'descend');
Vw= Vw(:,rearrangement_index);


%Then, we can perform whitening on Sw
Diagw_minusone = real(diag(1./sorted_eigval));
Diagw_minusone = Diagw_minusone(:,1:dim);

U = X*(eye(size(M))-M)*Vw*Diagw_minusone ;


%Now we can project the data
Xb = X'*U;

%Then we perform PCA on Xb
Q = PCA(Xb, length(unique(classes)) - 1);

%Finally, we get the total transform
U_reduc = real(U*Q);
end
