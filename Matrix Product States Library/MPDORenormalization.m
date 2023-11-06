function [outputMPDO] = MPDORenormalization(inputMPDO,Mode,TargetColumn)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm renormalizes a matrix product density 
% operator (MPDO) according to a specified mode. Unlike the renormalization
% of an MPS, this renormalization algorithm preserves the input MPDO norm.
% Last updated: November 2021.


%% Description of input and output:
% Input:
% inputMPDO is a matrix product density operator.

% Mode = 'RCN', 'LCN', 'MixedL', 'MixedR'
% If 'MixedL','MixedR' chosen, TargetColumn must also be declared.
% 'MixedL' renormalizes the matrices 1 through TargetColumn
% 'MixedR' renormalizes the matrices TargetColumn through end.
% 'LCN' performs a left-canonical renormalization.
% 'RCN' performs a right-canonical renormalization.

% Output:
% outputMPDO is the renormalized version of inputMPDO (with the operator norm
% preserved).


%% Mode Selector
% The following logical statement allows the user to call a specific
% version of the renormalization algorithm.
dVector = inputMPDO{1,end,1}(:,1);


% Length of chain
N = size(dVector,1);

if strcmp(Mode,'RCN')
    outputMPDO = RightNormalization(inputMPDO,1);
elseif strcmp(Mode, 'LCN')
    outputMPDO = LeftNormalization(inputMPDO,N);
elseif strcmp(Mode, 'MixedL')
    outputMPDO = LeftNormalization(inputMPDO,TargetColumn);
elseif strcmp(Mode, 'MixedR')
    outputMPDO = RightNormalization(inputMPDO,TargetColumn);
end

end


%% Normalization algorithms
function matrixStorage = RightNormalization(inputMatrixStorage,TargetColumn)

%% Initialization of parameters
% Load dVector and kVector:
dVector = inputMatrixStorage{1,end,1}(:,1);
kVector = inputMatrixStorage{1,end,1}(:,2);

% Length of chain:
N = size(dVector,1);

% Initialize output structure identical to input structure:
matrixStorage = inputMatrixStorage;

%% Normalization sweep from the right to TargetColumn
for cIndex = N:-1:TargetColumn
    di = dVector(cIndex);
    ki = kVector(cIndex);
    
    % Load MPDO matrices into a single matrix B; the matrices have to be
    % loaded in the same order as they are extracted at the end.
    B = [];
    
    for r2Index = 1:1:ki
        for r1Index = 1:1:di
            B = [B, matrixStorage{r1Index,cIndex,r2Index}];
        end
    end
    
    % Perform an SVD on the B matrix:
    [U0,S0,V0] = svd(B,'econ');
    V0Dagger = V0';
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
                matrixStorage(r1Index,cIndex-1,r2Index) = {matrixStorage{r1Index,cIndex-1,r2Index}*c0};
            end
        end
    end
    
    
    % Finally, divide V0Dagger into MPDO matrices and return the resulting
    % matrixStorage:
    sizeMeasure = size(V0Dagger);
    intervalSizeMinor = sizeMeasure(2)/(di*ki);
    intervalSizeMajor = sizeMeasure(2)/(ki);    
    
    for r1Index = 1:1:di
        for r2Index = 1:1:ki
            startIndex = (r1Index-1)*intervalSizeMinor + (r2Index-1)*intervalSizeMajor + 1;
            endIndex = r1Index*intervalSizeMinor + (r2Index-1)*intervalSizeMajor;
            matrixStorage(r1Index,cIndex,r2Index) = {V0Dagger(:,startIndex:1:endIndex)};
        end
    end
    
end

end



function matrixStorage = LeftNormalization(inputMatrixStorage,TargetColumn)

%% Initialization of parameters
% Load dVector and kVector:
dVector = inputMatrixStorage{1,end,1}(:,1);
kVector = inputMatrixStorage{1,end,1}(:,2);

% Length of chain:
N = size(dVector,1);

% Initialize output structure identical to input structure:
matrixStorage = inputMatrixStorage;

%% Normalization sweep from the left to TargetColumn
for cIndex = 1:1:TargetColumn
    di = dVector(cIndex);
    ki = kVector(cIndex);

    % Load MPO matrices into a single matrix A; the matrices have to be
    % loaded in the same order as they are extracted at the end.
    A = [];

    for r2Index = 1:1:ki
        for r1Index = 1:1:di
            
            A = [A; matrixStorage{r1Index,cIndex,r2Index}];
            
        end
    end
    
    % Perform an SVD on the matrix A:
    [U0,S0,V0] = svd(A, 'econ');
    c0 = S0*(V0');

    % Multiply c0 to the right OR, if we are at the final site, multiply c0
    % onto U0 (we preserve norm/phase when renormalizing MPOs):
    if cIndex == N
        U0 = U0;
    else
        dip1 = dVector(cIndex+1);
        kip1 = kVector(cIndex+1);
        for r1Index = 1:1:dip1
            for r2Index = 1:1:kip1
                matrixStorage(r1Index,cIndex+1,r2Index) = {c0*matrixStorage{r1Index,cIndex+1,r2Index}};
            end
        end
    end
    
    % Finally, divide U0 into MPO matrices and return the resulting
    % matrixStorage:
    sizeMeasure = size(U0);
    intervalSizeMinor = sizeMeasure(1)/(ki*di);
    intervalSizeMajor = sizeMeasure(1)/(ki);    
    
    for r1Index = 1:1:di
        for r2Index = 1:1:ki
            startIndex = (r1Index-1)*intervalSizeMinor + (r2Index-1)*intervalSizeMajor + 1;
            endIndex = r1Index*intervalSizeMinor + (r2Index-1)*intervalSizeMajor;
            matrixStorage(r1Index,cIndex,r2Index) = {U0(startIndex:1:endIndex,:)};
        end
    end
    
end

end