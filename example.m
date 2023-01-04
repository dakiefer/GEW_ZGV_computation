%% Compute ZGV points of guided elastic waves in an anisotropic plate
% Showcases three different computational techniques to locate ZGV points. 
% The matrices that describe the guided waves in the 1 mm thick plate are loaded
% from "matrices.mat" and are large enough (39x39) to account for frequencies up
% to 25 MHz.
%
% The methods have been presented in:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.

% % load matrices and parameters
load("matrices.mat") % matrices (L2, L1, L0, M), dispersion (dat), thickness (h), norm. param (np)
fmax = 25e6;         % maximum frequency of interest

% % here are some ZGV points with 2 digits accuracy (serve as initial guess for Newton's method): 
w0 = [0.29    0.52    0.69    0.72    0.75]*1e8; % circular frequencies
k0 = [0.34    0.35    0.71    0.25    0.07]*1e4; % wavenumbers

%% Newton-type iteration: 
% This method is super fast but requires initial guesses:
fprintf('\n\n++ Newton-type iteration: ++\nSearch close to initial guesses:\n')

% initialize:
kzgv = nan(length(k0),1);
wzgv = nan(length(w0),1); 
tic
for i=1:numel(w0) % for every initial guess (w0, k0)
    w0i = w0(i)*np.h0/np.fh0; k0i = k0(i)*np.h0; % normalize according to matrices
    [ki,wi] = ZGVNewtonBeta(L2, L1, L0, M, k0i, w0i);  % compute
    kzgv(i) = ki; wzgv(i) = wi; % store
end
toc
zgv.k = kzgv/np.h0; zgv.w = wzgv*np.fh0/np.h0; % save as structure in physical units

% print initial guesses and computed values 
disp('initial frequencies:')
disp(w0)
disp('converged frequencies:')
disp(zgv.w.')

% % plot
figure(1), clf, hold on, ylim([0, 25]); xlim([0, 25])
plot(dat.k*h, dat.w*h/2/pi/1e3);
plot(zgv.k(:)*h, zgv.w(:)*h/2/pi/1e3, 'r*');
xlabel('wavenumber-thickness kh in rad'), ylabel('frequency-thickness fh in MHz*mm')
title('Newton method: first five ZGV points'), drawnow;

%% Scanning the ZGV points
% This method does not need initial guesses and is likely to locate all ZGV
% points, but it is substantially slower than the Newton-type iteration. To speed
% up the computation, provide the paramter kMax (and optionally kStart), these
% define the wavenumber search interval [kStart, kMax].
fprintf('\n\n++ Scanning method: ++\n')

% method options:
clear opts, opts.kMax = 25e3*np.h0; % unit distance is np.h0

% use sparse matrices:
L0s = sparse(L0);
L1s = sparse(L1);
L2s = sparse(L2);
Ms = sparse(M);

tic; [kzgv, wzgv] = ZGV_MFRDScan(L2s, L1s, L0s, Ms, opts); toc % compute
zgv.k = kzgv/np.h0; zgv.w = wzgv*np.fh0/np.h0; % save as structure in physical units

disp('number of computed ZGV points:') % print number of located ZGV points 
disp(length(zgv.w(zgv.w < 2*pi*fmax))) 

% % plot
figure(2), clf, hold on, ylim([0, 25]); xlim([0, 25])
plot(dat.k*h, dat.w*h/2/pi/1e3); 
plot(zgv.k*h, zgv.w*h/2/pi/1e3, 'r*');
xlabel('wavenumber-thickness kh in rad'), ylabel('frequency-thickness fh in MHz*mm')
title('Scanning method'), drawnow;

%% Direct method
% This method does not need initial guesses and guarantees to find all ZGV
% points as long as the matrices defining the problem are not too big. This is
% rather slow and should not be used for matrices bigger than about 40x40.
fprintf('\n\n++ Direct method: ++\nThis will take a few minutes...\n')
if ~exist('threepar_delta', 'file') % check if function is on Matlab path
    error('To use the direct method, install MultiParEig from https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig');
end

tic, [kzgv, wzgv] = ZGVDirect(L2, L1, L0, M); toc % compute
zgv.k = kzgv/np.h0; zgv.w = wzgv*np.fh0/np.h0; % save as structure in physical units

disp('number of computed ZGV points:') % print number of located ZGV points 
disp(length(zgv.w(zgv.w < 2*pi*fmax))) 


% % plot
figure(3), clf, hold on, ylim([0, 25]); xlim([0, 25])
plot(dat.k*h, dat.w*h/2/pi/1e3);
plot(zgv.k*h, zgv.w*h/2/pi/1e3, 'r*');
xlabel('wavenumber-thickness kh in rad'), ylabel('frequency-thickness fh in MHz*mm')
title('Direct method'), drawnow;
