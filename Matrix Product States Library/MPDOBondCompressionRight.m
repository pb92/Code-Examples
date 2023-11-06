function [outputMPDO] = MPDOBondCompressionRight(inputMPDO, chimax, eps)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm compresses the matrix product density 
% operator (MPDO) inputMPDO using the SVD compression algorithm. The 
% resulting MPDO outputMPDO will be on a right-canonical form.
% Last updated: November 2021.

%% Description of input and output:
% Input:
% inputMPDO is a matrix product density operator (MPDO).
% chimax and eps are parameters for the compression algorithm, with chimax
% being the largest allowed matrix dimension, and eps being the cutoff
% point for weighs in the Schmidt decomposition.

% Output:
% outputMPDO is the resulting MPDO from the compression algorithm. 
% outputMPDO is on a right-canonical form.


%% Initialization of environment
% Load dVector and kVector:
dVector = inputMPDO{1,end,1}(:,1);
kVector = inputMPDO{1,end,1}(:,2);

% Read out the number of sites:
N = size(dVector,1);

% A left-renormalization is made before the truncation in order to ensure 
% an MPDO on right-canonical form.
outputMPDO = MPDORenormalization(inputMPDO,'LCN',[]);

% The truncation may then be performed from the right, leaving us a
% left-canonical form MPDO as the result. 


%% Compression algorithm
% We conduct a full-range compression sweep from right to left.
for cIndex = N:-1:1
    di = dVector(cIndex);
    ki = kVector(cIndex);

    % Load MPDO matrices into a single matrix B; the matrices have to be
    % loaded in the same order as they are extracted at the end.
    B = [];
    
    for r2Index = 1:1:ki
        for r1Index = 1:1:di
            B = [B, outputMPDO{r1Index,cIndex,r2Index}];
        end
    end
    
    % Perform an SVD on the B matrix:
    [U0,S0,V0] = svd(B,'econ');
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
    % by considering chimax and eps:
    NewDimension = sum(Weights >= eps); % might want to do a cumulative sum and make the error a sum instead.
    NewDimension = min(NewDimension,chimax);
    
    % Resize the matrices to reflect the number of dimensions to keep.
    U0 = U0(:,1:1:NewDimension);
    S0 = S0(1:1:NewDimension,1:1:NewDimension);
    V0Dagger = V0Dagger(1:1:NewDimension,:);
    
    % Calculate c0, which is the matrix to be multiplied to the right.
    c0 = U0*S0;
    
    % Multiply c0 to the left OR, if we are at the final site, multiply c0
    % onto V0Dagger (we preserve norm/phase when renormalizing MPOs):
    if cIndex == 1
        V0Dagger = V0Dagger;
    else
        dim1 = dVector(cIndex-1);
        kim1 = kVector(cIndex-1);
        for r1Index = 1:1:dim1
            for r2Index = 1:1:kim1
                outputMPDO(r1Index,cIndex-1,r2Index) = {outputMPDO{r1Index,cIndex-1,r2Index}*c0};
            end
        end
    end
    
    % Finally, divide V0Dagger into MPDO matrices and return the resulting
    % matrixStorage:
    sizeMeasure = size(V0Dagger);
    intervalSizeMinor = sizeMeasure(2)/(ki*di);
    intervalSizeMajor = sizeMeasure(2)/(ki);    
    
    for r1Index = 1:1:di
        for r2Index = 1:1:ki
            startIndex = (r1Index-1)*intervalSizeMinor + (r2Index-1)*intervalSizeMajor + 1;
            endIndex = r1Index*intervalSizeMinor + (r2Index-1)*intervalSizeMajor;
            outputMPDO(r1Index,cIndex,r2Index) = {V0Dagger(:,startIndex:1:endIndex)};
        end
    end
    
end

end