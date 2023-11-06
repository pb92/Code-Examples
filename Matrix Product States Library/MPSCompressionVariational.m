function [outputMPS, stateNorm] = MPSCompressionVariational(inputMPS, targetMPS, dVector, nmax, eps)

% Philip Daniel Blocher
% Center for Quantum information and Control, University of New Mexico
% blocher@unm.edu
% Matrix Product State (MPS) compression algorithm for open boundary
% conditions (OBC) that utilizes the variational method rather than SVD.
% Last updated: July 2021.

% Output MPS should be left-canonical.
% Rewrite

%% Description of input and output:
% Input:
% inputMPS is a MPS, corresponding to the structure described by dVector.
% dVector is a (L times 1) vector, where each entry details the degrees of
% freedom for each site in the L-site long chain.
% dmax and eps are parameters for the compression algorithm, with dmax
% being the largest allowed matrix dimension, and eps being the cutoff
% point for weighs in the Schmidt decomposition.

% Output:
% outputMPS is the resulting MPS from the compression algorithm.
% stateNorm is the norm of the input MPS (in case this deviates from 1).

%% Initialization of environment

% A left-renormalization is made before the truncation in order to ensure 
% a left canonical normalized state.
[inputMPS,~] = MPSRenormalization(inputMPS,dVector,'LCN',[]);

% The following two parameters are relevant for the compression algorithm:
N = size(dVector,1);
stateNorm = 1;

% Finally, we initialize the output MPS to be the inputMPS (such that we
% can manipulate it).
outputMPS = inputMPS; 

%% Variational compression algorithm

% We run the algorithm until sufficiently converged, or until we hit a
% pre-determined abort. This will not be implemented for the first run.

% Replace with a while statement: while not converged sufficiently and not
% at some cutoff repetition number
for nRun = 1:1:4

    % For each site i in the N-site chain we optimize the matrices:
    for i = 1:1:N

        % Renormalize outputMPS to be mixed canonical centered on site i:
        if i == 1
            %[outputMPS,~] = MPSRenormalization(outputMPS,dVector,'LCN',[]);
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'MixedR',i+1);
        elseif i == N
            %[outputMPS,~] = MPSRenormalization(outputMPS,dVector,'RCN',[]);
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'MixedL',i-1);
        else
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'LCN',[]);
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'MixedR',i+1);
        end
        
        % Calculate the L and R tensors:
        if i == 1
            L = 1;
        else
            LMatStorage = cell(1,i);
            LMatStorage(1) = {1};
            
            for j = 1:1:(i-1)
                
                dj = dVector(j);
                D = 0;
                
                for k = 1:1:dj
                    A = cell2mat(outputMPS(k,j));
                    B = cell2mat(targetMPS(k,j));
                    C = cell2mat(LMatStorage(j));
                    
                    D = D + A'*C*B;
                end
                
                LMatStorage(j+1) = {D};
                
            end
            
            L = LMatStorage{i};
            
        end
        
        
        if i == N
            R = 1;
        else
            RMatStorage = cell(1,N-i+1);
            RMatStorage(N) = {1};
            
            for j = N:-1:(i+1)
                
                dj = dVector(j);
                D = 0;
                
                for k = 1:1:dj
                    A = cell2mat(outputMPS(k,j));
                    B = cell2mat(targetMPS(k,j));
                    C = cell2mat(RMatStorage(j));
                    
                    D = D + B*C*(A');
                end
                
                RMatStorage(j-1) = {D};
                
            end
            
            R = RMatStorage{i};
            
        end
                
        % Calculate new M matrices for site i:
        di = dVector(i);
        
        Convergence2Norm = 1;
        
        for k = 1:1:di
            outputMPS{k,i} = L*(targetMPS{k,i}*R);
            
            Convergence2Norm = Convergence2Norm - trace((outputMPS{k,i}')*outputMPS{k,i});
        end
        
        %Convergence2Norm
        
    end   
    
% For each site i in the N-site chain we optimize the matrices:
    for i = N:-1:1

        % Renormalize outputMPS to be mixed canonical centered on site i:
        if i == 1
            %[outputMPS,~] = MPSRenormalization(outputMPS,dVector,'LCN',[]);
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'MixedR',i+1);
        elseif i == N
            %[outputMPS,~] = MPSRenormalization(outputMPS,dVector,'RCN',[]);
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'MixedL',i-1);
        else
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'LCN',[]);
            [outputMPS,~] = MPSRenormalization(outputMPS,dVector,'MixedR',i+1);
        end
        
        % Calculate the L and R tensors:
        if i == 1
            L = 1;
        else
            LMatStorage = cell(1,i);
            LMatStorage(1) = {1};
            
            for j = 1:1:(i-1)
                
                dj = dVector(j);
                D = 0;
                
                for k = 1:1:dj
                    A = cell2mat(outputMPS(k,j));
                    B = cell2mat(targetMPS(k,j));
                    C = cell2mat(LMatStorage(j));
                    
                    D = D + A'*C*B;
                end
                
                LMatStorage(j+1) = {D};
                
            end
            
            L = LMatStorage{i};
            
        end
        
        
        if i == N
            R = 1;
        else
            RMatStorage = cell(1,N-i+1);
            RMatStorage(N) = {1};
            
            for j = N:-1:(i+1)
                
                dj = dVector(j);
                D = 0;
                
                for k = 1:1:dj
                    A = cell2mat(outputMPS(k,j));
                    B = cell2mat(targetMPS(k,j));
                    C = cell2mat(RMatStorage(j));
                    
                    D = D + B*C*(A');
                end
                
                RMatStorage(j-1) = {D};
                
            end
            
            R = RMatStorage{i};
            
        end
                
        % Calculate new M matrices for site i:
        di = dVector(i);
        
        Convergence2Norm = 1;
        
        for k = 1:1:di
            outputMPS{k,i} = L*(targetMPS{k,i}*R);
            
            Convergence2Norm = Convergence2Norm - trace((outputMPS{k,i}')*outputMPS{k,i});
        end
        
        %Convergence2Norm
        
    end       
    
    
    
    
end

[outputMPS,~] = MPSRenormalization(outputMPS,dVector,'LCN',[]);

end

