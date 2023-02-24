%% Compute ZGV points of guided elastic waves in an anisotropic plate
% Showcases three different computational techniques to locate ZGV points. 
% The matrices that describe the guided waves in the 1 mm thick plate are loaded
% from "matrices_small.mat" and are large enough (25x25) to account for frequencies up
% to 18 MHz. The second example, "matrices_big.mat" represents the same plate with a 
% finer discretization (39x39) and can represent frequencies up to 25 MHz.
%
% The methods have been presented in:
% D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing
% zero-group-velocity points in anisotropic elastic waveguides: Globally and
% locally convergent methods," The Journal of the Acoustical Society of America,
% vol. 153, no. 2, pp. 1386–1398, Feb. 2023, doi: 10.1121/10.0017252.


% % load matrices and parameters
% L2, L1, L0, M : spectral element matrices
% dat           : dispersion data for visualization
% h             : plate thickness
% np            : normalization of matrices 
load("matrices_small.mat") 
wmax = 2*pi*18e6;         % maximum frequency of interest

% % here are some approx. ZGV points (serve as initial guess for Newton's method): 
w0 = [0.25    0.5    0.7    0.7    0.8]*1e8; % circular frequencies
k0 = [0.3    0.35    0.7    0.2    0.1]*1e4; % wavenumbers

% % plot dispersion curves for reference: 
fig = figure; clf, hold on, 
plot(dat.k*h, dat.w*h/2/pi/1e3, 'Color', [.7, .7, .7], 'HandleVisibility', 'off');
ylim([0, wmax/2/pi/1e6]), xlim([0, wmax/2/pi/1e6]), legend('Location','southeast')
set(fig,'defaulttextinterpreter','latex'), set(fig, 'DefaultLegendInterpreter', 'latex')
xlabel('wavenumber-thickness $kh$ in rad'), ylabel('frequency-thickness $fh$ in MHz$\cdot$mm')
drawnow;

%% Newton-type iteration: 
% This method is super fast but requires initial guesses:
fprintf('\n\n++ Newton-type iteration: ++\nSearch close to initial guesses:\n')

% initialize:
kzgvN = nan(length(k0),1);
wzgvN = nan(length(w0),1); 
tic
for i=1:numel(w0) % for every initial guess (w0, k0)
    w0i = w0(i)*np.h0/np.fh0; k0i = k0(i)*np.h0; % normalize according to matrices
    [ki,wi] = ZGVNewtonBeta(L2, L1, L0, M, k0i, w0i);  % compute
    kzgvN(i) = ki; wzgvN(i) = wi; % store
end
toc
zgvN.k = kzgvN/np.h0; zgvN.w = wzgvN*np.fh0/np.h0; % save as structure in physical units

% print initial guesses and computed values 
disp('initial frequencies:'), disp(w0/2/pi)
disp('converged frequencies:'), disp(zgvN.w.'/2/pi)

% % plot
figure(fig);
plot(zgvN.k(:)*h, zgvN.w(:)*h/2/pi/1e3, 'rx', 'DisplayName', 'Newton method');
drawnow;

%% Scanning the ZGV points
% This method does not need initial guesses and is likely to locate all ZGV
% points, but it is substantially slower than the Newton-type iteration. To speed
% up the computation, provide the parameter kMax (and optionally kStart), these
% define the wavenumber search interval [kStart, kMax].
fprintf('\n\n++ Scanning method: ++\n')

% method options:
clear opts, opts.kMax = 25e3*np.h0; % unit distance is np.h0

% use sparse matrices:
L0s = sparse(L0);
L1s = sparse(L1);
L2s = sparse(L2);
Ms = sparse(M);

tic; [kzgvS, wzgvS] = ZGV_MFRDScan(L2s, L1s, L0s, Ms, opts); toc % compute
zgvS.k = kzgvS/np.h0; zgvS.w = wzgvS*np.fh0/np.h0; % save as structure in physical units

fprintf('Number of computed ZGV points: %d\n', length(zgvS.w(zgvS.w < wmax))) 

% % plot
figure(fig);
plot(zgvS.k(:)*h, zgvS.w(:)*h/2/pi/1e3, 'ksquare', 'DisplayName', 'Scanning method', 'MarkerSize', 8);
drawnow;

%% Direct method
% This method does not need initial guesses and guarantees to find all ZGV
% points as long as the matrices defining the problem are not too big. This is
% rather slow and should not be used for matrices bigger than about 40x40.
% NOTE: Requires MultiParEig from 
% https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig
fprintf('\n\n++ Direct method: ++\nThis will take several seconds...\n')

tic, [kzgvD, wzgvD] = ZGVDirect(L2, L1, L0, M); toc % compute
zgvD.k = kzgvD/np.h0; zgvD.w = wzgvD*np.fh0/np.h0; % save as structure in physical units

fprintf('Number of computed ZGV points: %d\n', length(zgvD.w(zgvD.w < wmax)))

% % plot
figure(fig);
plot(zgvD.k(:)*h, zgvD.w(:)*h/2/pi/1e3, 'k.', 'DisplayName', 'Direct method');
drawnow;
