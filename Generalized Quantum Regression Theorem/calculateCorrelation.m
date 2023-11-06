% Philip Daniel Blocher
% Quantum Optics Group, Aarhus University
% pblocher@phys.au.dk - AU458007
% Last updated: May 2018.
clear all; close all; clc;
%% Variables:

% TLS parameters
Gamma = 1; % Decay in units of 1/s.
Omega = 2*Gamma;
delta = 3*Gamma;
dt = 1e-3;

% Radiation field
tDelay = 4;
m1t = sqrt(0.5);
m2t = sqrt(0.5);

m1r = sqrt(1 - m1t^2);
m2r = sqrt(1 - m2t^2);
beamsplitter = [m1t, m1r, m2t, m2r];

n = 0;

% EoM matrix:

M_coh = [0, 1i*Omega/2, -1i*Omega/2, 0;
     1i*Omega/2, -1i*delta, 0, -1i*Omega/2;
     -1i*Omega/2, 0, 1i*delta, 1i*Omega/2;
     0, -1i*Omega/2, 1i*Omega/2, 0];

M_incoh = [-Gamma*n, 0, 0, Gamma*(n+1);
           0, -(Gamma/2)*(2*n+1), 0, 0;
           0, 0, -(Gamma/2)*(2*n+1), 0;
           Gamma*n, 0, 0, -Gamma*(n+1);];
       
M = M_coh + M_incoh;
 
 
%% Determine the steady state of the system.

M_prime = [Gamma, 0, 0, Gamma; M(2:end,:)];

LHS_vec = [Gamma, 0, 0, 0]';

rho_steady = M_prime\LHS_vec;
rho_ss = reshape(rho_steady,2,2).'; 

rho_theo = Gamma*[delta^2 + Omega^2/4 + Gamma^2/4; 
                  delta*Omega/2+1i*Gamma*Omega/4;
                  delta*Omega/2-1i*Gamma*Omega/4;
                  Omega^2/4];
              
rho_theo = rho_theo*Gamma / (delta^2*Gamma^2+Gamma^4/4+Gamma^2*Omega^2/2);

%% Operators to be evaluated:

Id = eye(2);
sm = [0 1; 0 0];
sp = sm';

%% Noise
% Assign noise dependencies. Note that since we're ONLY interested in
% evaluating <sp sp sm sm> structures, I will only include noise terms
% needed for this evaluation.

N = zeros(16);

aad = Gamma*(n+1);
ada = Gamma*n;

% <P_g sp sm P_g> noise:
N(1,7)  =  aad; N(1,10)  =  ada;

% <P_g sp sm P_e> noise:
N(4,7)  = -aad; N(4,10)  = -ada;

% <sp sp sm sm> noise:
N(7,1)  =  ada; N(7,4)   = -ada; N(7,13)  = -ada; N(7,16)  =  ada;

% <sm sp sm sp> noise:
N(10,1) =  aad; N(10,4)  = -aad; N(10,13) = -aad; N(10,16) =  aad;

% <P_e sp sm P_g> noise:
N(13,7) = -aad; N(13,10) = -ada;

% <P_e sp sm P_e> noise:
N(16,7) = -aad; N(16,10) = -ada;

%% Evaluate four-operator correlation function
tDelayStart = 0;
tDelayEnd = 5;
dtDelay = 1/10;
tDelays = tDelayStart:dtDelay:tDelayEnd;

%tDelays = [0:1/10:1, 1.25:1/4:5];

tDelays = [0]


tau_start = 0;
tau_end = 5;
dtau = 1/50;
taus = tau_start:dtau:tau_end;

detectorARecord = NaN(size(taus,2), size(tDelays,2));
detectorARecordNoNoise = NaN(size(taus,2), size(tDelays,2));


tic
%parpool('local',30)
parfor index_tDelay = 1:1:size(tDelays,2)
    
    vA = NaN(size(taus,2),1);
    vAnoNoise = NaN(size(taus,2),1);
    
    for index_tau = 1:1:size(taus,2)

        tDelay = tDelays(index_tDelay);
        tau = abs(taus(index_tau));
        tau2 = taus(index_tau);

        valueA = correlationFunction(rho_ss, M, N, dt, sp, sp, sm, sm, [0, tau, tau, 0], tDelay, beamsplitter);
        valueAnoNoise = correlationFunction(rho_ss, M, zeros(size(N)), dt, sp, sp, sm, sm, [0, tau, tau, 0], tDelay, beamsplitter);


        vA(index_tau) = valueA;
        vAnoNoise(index_tau) = valueAnoNoise;
        
    end
    
    detectorARecord(:,index_tDelay) = vA;
    detectorARecordNoNoise(:,index_tDelay) = vAnoNoise;

    
end
toc

%% Export data

save('g2Test.mat');

%% Support functions

function returnValue = correlationFunction(rho_ss, M, N, dt, A, B, C, D, tVec, tDelay, beamsplitter)
  m1t = beamsplitter(1); 
  m1r = beamsplitter(2); 
  m2t = beamsplitter(3);
  m2r = beamsplitter(4);

  returnValue = 0;

% 'phase-augmented' beamsplitter setup
  returnValue = returnValue + (m1t*m2t)^4*            fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,0,0]);
  returnValue = returnValue + (m1t*m2t)^3*(m1r*m2r)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,0,1]);
  returnValue = returnValue + (m1t*m2t)^3*(m1r*m2r)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,1,0]);
  returnValue = returnValue + (m1t*m2t)^2*(m1r*m2r)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,1,1]);

  returnValue = returnValue + (m1t*m2t)^3*(m1r*m2r)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,0,0]);
  returnValue = returnValue + (m1t*m2t)^2*(m1r*m2r)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,0,1]);
  returnValue = returnValue + (m1t*m2t)^2*(m1r*m2r)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,1,0]);
  returnValue = returnValue + (m1t*m2t)*(m1r*m2r)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,1,1]);

  returnValue = returnValue + (m1t*m2t)^3*(m1r*m2r)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,0,0]);
  returnValue = returnValue + (m1t*m2t)^2*(m1r*m2r)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,0,1]);
  returnValue = returnValue + (m1t*m2t)^2*(m1r*m2r)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,1,0]);
  returnValue = returnValue + (m1t*m2t)*(m1r*m2r)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,1,1]);

  returnValue = returnValue + (m1t*m2t)^2*(m1r*m2r)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,0,0]);
  returnValue = returnValue + (m1t*m2t)*(m1r*m2r)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,0,1]);
  returnValue = returnValue + (m1t*m2t)*(m1r*m2r)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,1,0]);
  returnValue = returnValue + (m1r*m2r)^4*            fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,1,1]);
end

function returnValue = correlationFunctionSecondary(rho_ss, M, N, dt, A, B, C, D, tVec, tDelay, beamsplitter)
  m1t = beamsplitter(1);
  m1r = beamsplitter(2);
  m2t = beamsplitter(3);
  m2r = beamsplitter(4);

  returnValue = 0;

  % 'phase-augmented' beamsplitter setup
  returnValue = returnValue + (m1t*m2r)^4*            fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,0,0]);
  returnValue = returnValue - (m1t*m2r)^3*(m1r*m2t)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,0,1]);
  returnValue = returnValue - (m1t*m2r)^3*(m1r*m2t)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,1,0]);
  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,0,1,1]);

  returnValue = returnValue - (m1t*m2r)^3*(m1r*m2t)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,0,0]);
  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,0,1]);
  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,1,0]);
  returnValue = returnValue - (m1t*m2r)*(m1r*m2t)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,1,1]);

  returnValue = returnValue - (m1t*m2r)^3*(m1r*m2t)*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,0,0]);
  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,0,1]);
  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,1,0]);
  returnValue = returnValue - (m1t*m2r)*(m1r*m2t)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,0,1,1]);

  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,0,0]);
  returnValue = returnValue - (m1t*m2r)*(m1r*m2t)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,0,1]);
  returnValue = returnValue - (m1t*m2r)*(m1r*m2t)^3*  fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,1,0]);
  returnValue = returnValue + (m1r*m2t)^4*            fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[1,1,1,1]);
end

function returnValue = singleOTOC(rho_ss, M, N, dt, A, B, C, D, tVec, tDelay, beamsplitter)
  m1t = beamsplitter(1); 
  m1r = beamsplitter(2); 
  m2t = beamsplitter(3);
  m2r = beamsplitter(4);

  returnValue = 0;

  % 'phase-augmented' beamsplitter setup
  returnValue = returnValue + (m1t*m2r)^2*(m1r*m2t)^2*fourTimeCorrelationNoise(rho_ss, M, N, dt, A, B, C, D, tVec - tDelay*[0,1,0,1]);

end