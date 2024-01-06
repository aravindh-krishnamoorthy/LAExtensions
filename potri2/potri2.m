%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is a part of the package: MatrixAlgorithms
% Released under the MIT license, see LICENSE file for details.
% Copyright (C) 2023 Aravindh Krishnamoorthy and contributors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = potri2(X)

n = size(X,1) ;

% Cholesky decomposition
R = chol(X) ;

% Inversion
D = diag(1./diag(R)) ;
R(n,1:n) = conj(R(1:n,1:n)\D(1:n,n)) ;
for i=n-1:-1:1
    R(i,1:i) = conj(R(1:i,1:i)\(D(1:i,i) - R(1:i,i+1:n)*R(i+1:n,i))) ;
end
Y = tril(R) + tril(R,-1)' ;

end
