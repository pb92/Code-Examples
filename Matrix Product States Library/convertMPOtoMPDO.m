function [outputMPDO] = convertMPOtoMPDO(MPO,dVector)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm converts a density matrix expressed as a 
% matrix product operator (MPO) into a matrix product density operator 
% (MPDO).
% Last updated: November 2021.

%% Description of input and output:
% Input:
% inputMPO is the matrix product operator array. The matrices of inputMPO
% are required to have dimensions that are squares of integers.
% dVector is a (N times 1)-vector containing the degrees of freedom for
% each site.

% Output:
% outputMPDO is the matrix product density operator array.


%% Initialize working environment
% Length of chain:
N = length(dVector);

% Create the kVector describing the dimensions of the Kraus indices. The
% vector entries are calculated in the decomposition algorithm.
kVector = ones(N,1);

% To create the outputMPDO object, we need to know the maximum dimension 
% for the Kraus indices:
maxKrausDimension = 1;
for siteIndex = 1:1:N
    maxKrausDimension = max(maxKrausDimension, dVector(siteIndex)*size(MPO{1,siteIndex,1},1)*size(MPO{1,siteIndex,1},2));
end

% Initialize the outputMPO array according to the system dimensions:
outputMPDO = cell(max(dVector),N+1,maxKrausDimension);


%% Convert MPO matrices into MPDO matrices
for siteIndex = 1:1:N
    
    di = dVector(siteIndex);
    
    % Read out the MPO matrix dimensions for site siteIndex:
    [xDimSquared, yDimSquared] = size(MPO{1,siteIndex,1});
    
    % Calculate the corresponding MPDO matrix dimensions:
    xDim = sqrt(xDimSquared);
    yDim = sqrt(yDimSquared);
    
    % To convert MPO matrices into MPDO matrices, we need to reorder
    % indices to from M^{sigma,sigma'}_{x x',y y'} to 
    % M_{sigma x y, sigma' x' y'}. This is done by defining a matrix C and
    % through four for-loops re-ordering the indices:
    C = NaN(di*xDim*yDim,di*xDim*yDim);
    
    for sigma = 1:1:di
        for sigmaprime = 1:1:di
            
            % Read out the (sigma,sigmaprime) matrix from the MPO:
            M = MPO{sigma,siteIndex,sigmaprime};
            
            for xIndex = 1:1:xDim
                for xprimeIndex = 1:1:xDim
                    for yIndex = 1:1:yDim
                        for yprimeIndex = 1:1:yDim
                            
                            % Calculate the indices in the old ordering:
                            xIndexIn = (xIndex-1)*xDim + xprimeIndex;
                            yIndexIn = (yIndex-1)*yDim + yprimeIndex;
                            
                            % Calculate the indices in the new ordering:
                            xIndexOut = (sigma-1)*xDim*yDim + (xIndex-1)*yDim + yIndex;
                            yIndexOut = (sigmaprime-1)*xDim*yDim + (xprimeIndex-1)*yDim + yprimeIndex;
                            
                            % Now assign the corresponding entries to the
                            % C-matrix:
                            C(xIndexOut,yIndexOut) = M(xIndexIn,yIndexIn);
                            
                        end
                    end
                end
            end
        end
    end

    % Do an SVD on the C-matrix:
    [U, S, V] = svd(C,'econ');
    
    % Define new matrices U and VDag by multiplying them with sqrtm(S). The
    % new U (VDag) will be used to create the MPDO matrices
    U = U*sqrtm(S);
    VDag = sqrtm(S)*(V');
    
    % Update the Kraus dimension in kVector:
    kVector(siteIndex) = length(diag(S));
    
    % Now split the U matrix into MPDO matrices:
    for sigma = 1:1:di
        for k = 1:1:kVector(siteIndex)
            
            A = NaN(xDim,yDim);
            
            for xIndex = 1:1:xDim
                for yIndex = 1:1:yDim
                    
                    xIndexIn = (sigma-1)*xDim*yDim + (xIndex-1)*yDim + yIndex;
                    A(xIndex,yIndex) = U(xIndexIn,k);
                    
                end
            end
            
            outputMPDO{sigma,siteIndex,k} = A;
            
        end
    end
end
                    
% Update the final entry of the MPDO with the dVector and kVector:
outputMPDO{1,end,1} = [dVector,kVector];

end