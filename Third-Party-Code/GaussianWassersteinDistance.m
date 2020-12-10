function [gwd] = GaussianWassersteinDistance(m1,P1,m2,P2)

% Computes the Wasserstein Distance between to Gaussian distributions with
% mean m1 and m2 and covariance P1 and P2.
%
% Can be used as a performance measure for extended target tracking under
% assumed ellipse shapes. Input m1 and m2 as the kinematic vectors (i.e.
% position, velocity, etc.) and input P1 and P2 as the shape matrices.

% Code by Karl granstrom

%Modification by Angel Garcia to avoid warning with zero matrices

if(~any(P2,'all') || ~any(P1,'all'))
    %If P1 or P2 is zero
    gwd = norm(m1-m2,2) + trace(P1 + P2);
    
else
    sqrtP1 = sqrtm(P1);
    gwd = norm(m1-m2,2) + trace(P1 + P2 - 2*sqrtm(sqrtP1*P2*sqrtP1));
end