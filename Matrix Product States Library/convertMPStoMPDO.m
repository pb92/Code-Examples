function [outputMPDO] = convertMPStoMPDO(inputMPS,dVector)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm turns a matrix product state (MPS) |Psi> into
% a matrix product density operator (MPDO) rho = |Psi><Psi| and returns the
% MPDO.
% Last updated: September 2021.

% Notes: DOES THE MPDO SATISFY APPROPRIATE PROPERTIES???

%% Description of inputs and output
% Input:
% inputMPS is a (dmax times N)-array, with each entry being a matrix. N is
% the number of sites in the 1D chain, while dmax is the largest number of
% degrees of freedom for a site in the 1D chain.
% dVector is a (N times 1)-vector containing the degrees of freedom for
% each site.

% Output:
% outputMPDO outputs the resulting MPDO.


%% Initialization of environment

% Define variables:
dmax = max(dVector);
N = size(dVector, 1);

% Initialize output storage:
outputMPDO = cell(dmax,N,dmax);


%% Create the output MPDO from the input MPS

for n = 1:1:N
    dn = dVector(n);
    
    for sigma = 1:1:dn
        
        for sigmaprime = 1:1:dn
            
            M = cell2mat(inputMPS(sigma,n));
            W = conj(double(cell2mat(inputMPS(sigmaprime,n))));
            
            outputMPDO(sigma,n,sigmaprime) = {kron(M,W)};            
        end
    end
end

end