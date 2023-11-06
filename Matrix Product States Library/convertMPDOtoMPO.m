function outputMPO = convertMPDOtoMPO(MPDO)

% Philip Daniel Blocher
% Center for Quantum Information and Control
% Department of Physics and Astronomy, University of New Mexico
% Email: blocher@unm.edu (current mail on Google Scholar)

% Description: This algorithm converts a matrix product density operator
% (MPDO) into a matrix product operator (MPO).
% Last updated: November 2021.

%% Description of input and output:
% Input:
% inputMPDO is the matrix product density operator array.

% Output:
% outputMPO is the matrix product operator array.


%% Initialize working environment
% Load dVector and kVector:
dVector = MPDO{1,end,1}(:,1);
kVector = MPDO{1,end,1}(:,2);

% Length of chain:
N = length(dVector);

% Initialize the outputMPO array according to the system dimensions:
outputMPO = cell(max(dVector),length(dVector),max(dVector));


%% Convert MPDO matrices into MPO matrices
for siteIndex = 1:1:N
    
    di = dVector(siteIndex);
    ki = kVector(siteIndex);
    
    for sigma = 1:1:di  
        for sigmaprime = 1:1:di
            
            % Read out the MPDO matrices corresponding to both the primed
            % and unprimed physical indices (sigma, sigmaprime), for the
            % first Kraus index k = 1:
            A_mat = MPDO{sigma,siteIndex,1};
            Aprime_mat = MPDO{sigmaprime,siteIndex,1};
            
            % The MPDO matrices are recovered as the Kronecker tensor
            % product between the MPDO matrices:            
            M_mat = kron(A_mat,conj(Aprime_mat));
            
            % and need to be summed over the different Kraus indices:
            for k = 2:1:ki
                
                A_mat = MPDO{sigma,siteIndex,k};
                Aprime_mat = MPDO{sigmaprime,siteIndex,k};
                
                M_mat = M_mat + kron(A_mat,conj(Aprime_mat));
                
            end
            
            % Store the MPO matrix in the outputMPO array:            
            outputMPO{sigma,siteIndex,sigmaprime} = M_mat;
            
        end
    end
            
end

end