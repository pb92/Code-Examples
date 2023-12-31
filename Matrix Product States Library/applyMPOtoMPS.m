function [outputMPS] = applyMPOtoMPS(inputMPO,inputMPS,dVector)

% Philip Daniel Blocher
% Quantum Optics Group
% Department of Physics and Astronomy, Aarhus University
% pblocher@phys.au.dk
% Algorithm for applying a matrix product operator (MPO) to a matrix
% product state (MPS). This yields a new MPS that is then returned to the
% caller.
% Last updated: November 2017.

%% Description of input and output:
% Input:
% inputMPO is a (dmax times L times dmax)-array, with each entry being one
% of the matrices making up the matrix product operator.
% inputMPS is a (dmax times L)-array where each entry is one matrix.
% L is the number of sites in the 1D chain.
% dmax is the largest number of degrees of freedom for a site in the 1D
% chain.
% dVector is a (L times 1)-vector with the number of degrees of freedom 
% for each site. 

% Output:
% outputMPS outputs the result of applying the MPO to the MPS. This
% resulting array is again an MPS.

%% Initialization of environment:

% Define variables:
dmax = max(dVector);
L = size(dVector, 1);

% Initialize output storage:
outputMPS = cell(dmax,L);

%% Algorithm:
for k = 1:1:L
    dk = dVector(k);
    for i = 1:1:dk % Sigma                    
        Mtilde = 0;
        for j = 1:1:dk % Sigma'
            N = cell2mat(inputMPO(i,k,j));
            M = cell2mat(inputMPS(j,k));
            Mtilde = Mtilde + kron(N,M);
        end
        outputMPS(i,k) = {Mtilde};
    end
end

end