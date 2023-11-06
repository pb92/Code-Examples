function [outputMPO] = MPOCompressionLeft(inputMPO, dVector, dmax, eps)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm compresses the matrix product operator (MPO)
% inputMPO using the SVD compression algorithm. The resulting MPO outputMPO
% will be on a left-canonical form.
% Last updated: August 2021.

%% Description of input and output:
% Input:
% inputMPO is an MPO, corresponding to the structure described by dVector.
% dVector is a (N times 1) vector, where each entry details the degrees of
% freedom for each site in the N-site long chain.
% dmax and eps are parameters for the compression algorithm, with dmax
% being the largest allowed matrix dimension, and eps being the cutoff
% point for weighs in the Schmidt decomposition.

% Output:
% outputMPO is the resulting MPO from the compression algorithm. outputMPO
% is on a left-canonical form.


%% Initialization of environment
% Read out the number of sites:
N = size(dVector,1);

% A right-renormalization is made before the truncation in order to ensure 
% an MPO on right-canonical form.
outputMPO = MPORenormalization(inputMPO,dVector,'RCN',[]);

% The truncation may then be performed from the left, leaving us a
% left-canonical form MPO as the result. 


%% Compression algorithm
% We conduct a full-range compression sweep from left to right.
for cIndex = 1:1:N
    di = dVector(cIndex);

    % Load MPO matrices into a single matrix A; the matrices have to be
    % loaded in the same order as they are extracted at the end.
    A = [];

    for r2Index = 1:1:di
        for r1Index = 1:1:di
            
            A = [A; outputMPO{r1Index,cIndex,r2Index}];
            
        end
    end
    
    % Now perform an SVD on A to get U0, S0, and V0Dagger.
    [U0,S0,V0] = svd(A,'econ');
    V0Dagger = V0';
    
    % In order to do the compression, we need to calculate the Schmidt
    % weights:
    S2 = S0 * S0;
    if diag(S2) == 0 % Fix for cases where S2 is the zero matrix...
        Weights = ones(size(S2));
    else
        Weights = diag(S2)./S2(1,1);
    end
    
    % We now need to select the number of dimensions to keep; this is done
    % by considering dmax and eps:
    NewDimension = sum(Weights >= eps); % might want to do a cumulative sum and make the error a sum instead.
    NewDimension = min(NewDimension,dmax);
    
    % Resize the matrices to reflect the number of dimensions to keep.
    U0 = U0(:,1:1:NewDimension);
    S0 = S0(1:1:NewDimension,1:1:NewDimension);
    V0Dagger = V0Dagger(1:1:NewDimension,:);
    
    % Calculate c0, which is the matrix to be multiplied to the right.
    c0 = S0*V0Dagger;
    
    % Multiply c0 to the right OR, if we are at the final site, multiply c0
    % onto U0 (we preserve norm/phase when renormalizing MPOs):
    if cIndex == N
        U0 = U0*c0;
    else
        dip1 = dVector(cIndex+1);
        for r1Index = 1:1:dip1
            for r2Index = 1:1:dip1
                outputMPO(r1Index,cIndex+1,r2Index) = {c0*outputMPO{r1Index,cIndex+1,r2Index}};
            end
        end
    end
     
    % Finally, divide U0 into MPO matrices and return the resulting
    % matrixStorage:
    sizeMeasure = size(U0);
    intervalSizeMinor = sizeMeasure(1)/(di^2);
    intervalSizeMajor = sizeMeasure(1)/(di^1);    
    
    for r1Index = 1:1:di
        for r2Index = 1:1:di
            startIndex = (r1Index-1)*intervalSizeMinor + (r2Index-1)*intervalSizeMajor + 1;
            endIndex = r1Index*intervalSizeMinor + (r2Index-1)*intervalSizeMajor;
            outputMPO(r1Index,cIndex,r2Index) = {U0(startIndex:1:endIndex,:)};
        end
    end
    
end
    
end