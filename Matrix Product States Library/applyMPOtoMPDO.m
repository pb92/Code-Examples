function [outputMPDO] = applyMPOtoMPDO(inputMPO,inputMPDO)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm applies a matrix product operator (MPO) to a 
% MPDO, which yields a new MPDO that is then returned.
% Last updated: December 2021.

%% Description of inputs and output:
% Input:
% inputMPO is a (dmax times N times dmax)-array, with each entry being a 
% matrix.
% inputMPDO is likewise a (dmax times N times kmax)-array.
% N is the number of sites in the 1D chain, while dmax is the largest
% number of degrees of freedom for a site in the 1D chain. kmax is the
% largest Kraus dimension in the MPDO.

% Output:
% outputMPDO outputs the resulting MPDO from applying inputMPO to 
% inputMPDO. The resulting array is again an MPDO.

%% Initialization of environment:
% Read out sizes:
dVector = inputMPDO{1,end,1}(:,1);
kVector = inputMPDO{1,end,1}(:,2);

% Define variables:
N = size(dVector, 1);

% Initialize output storage:
outputMPDO = cell(max(dVector),N+1,max(kVector));

%% Algorithm:

for siteIndex = 1:1:N
    
    dSite = dVector(siteIndex);
    kSite = kVector(siteIndex);
    
    for sigma = 1:1:dSite
        
        for a = 1:1:kSite
        
            Mtilde = 0;
            for sigmaprime = 1:1:dSite
                W = cell2mat(inputMPO(sigma,siteIndex,sigmaprime));
                M = cell2mat(inputMPDO(sigmaprime,siteIndex,a));
                Mtilde = Mtilde + kron(W,M);
            end
            
            outputMPDO(sigma,siteIndex,a) = {Mtilde};
            
        end
    end
end

outputMPDO{1,end,1} = inputMPDO{1,end,1};

end