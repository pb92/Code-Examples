function [outputMPDO] = MPDOKrausCompression(inputMPDO, kappamax, eps)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm compresses the Kraus dimension of the 
% matrix product density operator (MPDO) inputMPDO using the SVD 
% compression algorithm.
% Last updated: November 2021.

%% Description of input and output:
% Input:
% inputMPDO is a matrix product density operator (MPDO).
% dmax and eps are parameters for the compression algorithm, with kappamax
% being the largest allowed Kraus dimension, and eps being the cutoff
% point for weighs in the Schmidt decomposition.

% Output:
% outputMPDO is the resulting MPDO from the compression algorithm. 


%% Initialization of environment
% Load dVector and kVector:
dVector = inputMPDO{1,end,1}(:,1);
kVector = inputMPDO{1,end,1}(:,2);

% Read out the number of sites:
N = size(dVector,1);

% Define a new kVector for the compressed MPDO:
kVectorCompressed = NaN(size(kVector));

% A right-renormalization is made before the truncation in order to ensure 
% an MPDO on left-canonical form.
outputMPDO = MPDORenormalization(inputMPDO,'RCN',[]);

% The truncation may then be performed from the left. The right-canonical
% form should not be necessary.


%% Compression algorithm
% We conduct a full-range compression sweep from left to right.
for cIndex = 1:1:N
    di = dVector(cIndex);
    ki = kVector(cIndex);
    
    % Read out the bond dimensions of the MPDO matrices at site cIndex:
    [xDim,yDim] = size(inputMPDO{1,cIndex,1});
    
    % Load MPDO matrices into a single matrix C; the matrices have to be
    % loaded in a particular order to ensure truncation of the Kraus
    % dimension.
    C = NaN(di*xDim*yDim,ki);
    
    for sigma = 1:1:di
        for k = 1:1:ki
            
            % Load the matrix corresponding to (sigma,k):
            A = outputMPDO{sigma,cIndex,k};
            
            for xIndex = 1:1:xDim
                for yIndex = 1:1:yDim
                    
                    % Outgoing indices:
                    xIndexOut = (sigma-1)*xDim*yDim + (xIndex-1)*yDim + yIndex;
                    yIndexOut = k;
                    
                    % Assign element to C matrix:
                    C(xIndexOut,yIndexOut) = A(xIndex,yIndex);
                    
                end
            end
        end
    end
    
    % Now perform an SVD on C to get S0 and V0Dagger (we can ignore 
    % V0Dagger due to gauge freedom):
    [U0,S0,V0] = svd(C,'econ');
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
    % by considering kappamax and eps:
    NewDimension = sum(Weights >= eps); % might want to do a cumulative sum and make the error a sum instead.
%     NewDimension = min(NewDimension,dVector(cIndex)*min(prod(dVector(1:cIndex-1)),prod(dVector(cIndex:end)))*min(prod(dVector(1:cIndex)),prod(dVector(cIndex+1:end))));
%     NewDimensionErrorKraus = sum(Weights(NewDimension+1:1:end))
    
    NewDimension = min(NewDimension,kappamax);
    
    % Resize the matrices to reflect the number of dimensions to keep.
    U0 = U0(:,1:1:NewDimension);
    S0 = S0(1:1:NewDimension,1:1:NewDimension);
    V0Dagger = V0Dagger(1:1:NewDimension,:);
    
    % Update the kVector for the compressed MPDO:
    kVectorCompressed(cIndex) = NewDimension;
    
    % Calculate c0, which is the matrix that is retained and will be sliced
    % into new matrices.
    C0 = U0*S0;
     
    % Finally, recover from C0 the MPDO matrices with a truncated Kraus
    % dimension:
    for sigma = 1:1:di
        for k = 1:1:NewDimension
            
            % Calculate the new MPDO matrix for (sigma,k):
            A = NaN(xDim,yDim);
            
            for xIndex = 1:1:xDim
                for yIndex = 1:1:yDim
                    
                    xIndexIn = (sigma-1)*xDim*yDim + (xIndex-1)*yDim + yIndex;
                    yIndexIn = k;
                    
                    A(xIndex,yIndex) = C0(xIndexIn,yIndexIn);
                    
                end
            end
            
            % Store the matrix in the outputMPDO array:            
            outputMPDO{sigma,cIndex,k} = A;
            
        end
    end
end

% Update the size of outputMPDO to reflect the compression:
outputMPDO = outputMPDO(:,:,1:1:max(kVectorCompressed));

% Update the final entry of the outputMPDO
outputMPDO{1,end,1} = [dVector, kVectorCompressed];
    
end