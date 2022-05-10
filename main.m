load data/PANdata_P45.mat; % load data from Panasonic NMC cell
T = 45; % Test temperature

time    = DYNData.script1.time(:);   deltat = time(2)-time(1);
time    = time-time(1); 
current = DYNData.script1.current(:); 
voltage = DYNData.script1.voltage(:);
soc     = DYNData.script1.soc(:);

% Load cell-test data to be used for this batch experiment 
% (script 1 with

load data/PANmodel.mat; % load ESC model of Panasonic NMC cell
addpath helper_function\
% Reserve storage for computed results, for plotting
soc_hat = zeros(size(soc));
soc_bound = zeros(size(soc));

% Covariance values
SigmaX0 = diag([1e2 1e-2 1e-3]); 
SigmaV = 3e-1; %Voltage Sensor Noise
SigmaW = 4e0; %Current Sensor Noise

%% Data Initialization

spkfData.irInd = 1;
spkfData.hkInd = 2;
spkfData.zkInd = 3;

ir0   = 0;                           
hk0   = 0;                           
SOC0  = SOCfromOCVtemp(voltage(1),T,model); 
spkfData.xhat  = [ir0; hk0; SOC0]; % initial state

% Covariance values
spkfData.SigmaX = SigmaX0;
spkfData.SigmaV = SigmaV;
spkfData.SigmaW = SigmaW;
spkfData.Snoise = real(chol(diag([SigmaW; SigmaV]),'lower'));
spkfData.Qbump = 5;
 
% SPKF specific parameters
Nx = length(spkfData.xhat); spkfData.Nx = Nx; % state-vector length
Ny = 1; spkfData.Ny = Ny; % measurement-vector length
Nu = 1; spkfData.Nu = Nu; % input-vector length
Nw = size(SigmaW,1); spkfData.Nw = Nw; % process-noise-vector length
Nv = size(SigmaV,1); spkfData.Nv = Nv; % sensor-noise-vector length
Na = Nx+Nw+Nv; spkfData.Na = Na;     % augmented-state-vector length

h = sqrt(3);
spkfData.h = h; % SPKF/CDKF tuning factor  
Weight1 = (h*h-Na)/(h*h); % weighting factors when computing mean
Weight2 = 1/(2*h*h);      % and covariance
spkfData.Wm = [Weight1; Weight2*ones(2*Na,1)]; % mean
spkfData.Wc = spkfData.Wm;                     % covar

% previous value of current
spkfData.priorI = 0;
spkfData.signIk = 0;
  
% store model data structure too
spkfData.model = model;

%% Iteration
hwait = waitbar(0,'Processing');
for k = 1:length(voltage)
  %Measuring
  vk = voltage(k); 
  ik = current(k); 
  Tk = T;         

  % Load the cell model parameters
  Q  = getParamESC('QParam',Tk,model);
  G  = getParamESC('GParam',Tk,model);
  M  = getParamESC('MParam',Tk,model);
  M0 = getParamESC('M0Param',Tk,model);
  RC = exp(-deltat./abs(getParamESC('RCParam',Tk,model)))';
  R  = getParamESC('RParam',Tk,model)';
  R0 = getParamESC('R0Param',Tk,model);
  eta = getParamESC('etaParam',Tk,model);
  if ik<0, ik=ik*eta; end
  
  % Data from spkfData
  I = spkfData.priorI;
  SigmaX = spkfData.SigmaX;
  xhat = spkfData.xhat;
  Nx = spkfData.Nx;
  Nw = spkfData.Nw;
  Nv = spkfData.Nv;
  Na = spkfData.Na;
  Snoise = spkfData.Snoise;
  Wc = spkfData.Wc;
  irInd = spkfData.irInd;
  hkInd = spkfData.hkInd;
  zkInd = spkfData.zkInd;
  if abs(ik)>Q/100, spkfData.signIk = sign(ik); end
  signIk = spkfData.signIk;

  % KF Preparation (Sigma Points)
  [sigmaXa,p] = chol(SigmaX,'lower'); 
  if p>0
    fprintf('Cholesky error.  Recovering...\n');
    theAbsDiag = abs(diag(SigmaX));
    sigmaXa = diag(max(sqrt(max(0,theAbsDiag)),sqrt(max(0,spkfData.SigmaW))));
  end
  sigmaXa=[real(sigmaXa) zeros([Nx Nw+Nv]); zeros([Nw+Nv Nx]) Snoise];
  xhata = [xhat; zeros([Nw+Nv 1])];
  Xa = xhata(:,ones([1 2*Na+1]))+spkfData.h*[zeros([Na 1]), sigmaXa, -sigmaXa];

  % KF Prediction
  % State Equation Calculation
  xold = Xa(1:Nx,:);
  xnoise = Xa(Nx+1:Nx+Nw,:);
  current_n = I + xnoise; % noise adds to current
  Xx = 0*xold;
  Xx(irInd,:) = RC*xold(irInd,:) + (1-diag(RC))*current_n;
  Ah = exp(-abs(current_n*G*deltat/(3600*Q)));  % hysteresis factor
  Xx(hkInd,:) = Ah.*xold(hkInd,:) + (Ah-1).*sign(current_n);
  Xx(zkInd,:) = xold(zkInd,:) - current_n/3600/Q;
  Xx(hkInd,:) = min(1,max(-1,Xx(hkInd,:)));
  Xx(zkInd,:) = min(1.05,max(-0.05,Xx(zkInd,:)));

  xhat = Xx*spkfData.Wm;
  Xs = Xx - xhat(:,ones([1 2*Na+1]));
  SigmaX = Xs*diag(Wc)*Xs';
  I = ik; yk = vk;
  ynoise = Xa(Nx+Nw+1:end,:);
  Y = OCVfromSOCtemp(Xx(zkInd,:),Tk,model);
  Y = Y + M*Xx(hkInd,:) + M0*signIk;
  Y = Y - R*Xx(irInd,:) - R0*I + ynoise(1,:);
  yhat = Y*spkfData.Wm;
  
  %KF Correction
  Ys = Y - yhat(:,ones([1 2*Na+1]));
  SigmaXY = Xs*diag(Wc)*Ys';
  SigmaY = Ys*diag(Wc)*Ys';
  L = SigmaXY/SigmaY; 

  r = yk - yhat;
  if r^2 > 100*SigmaY, L(:,1)=0.0; end 
  xhat = xhat + L*r; 
  xhat(zkInd)=min(1.05,max(-0.05,xhat(zkInd)));
  xhat(hkInd) = min(1,max(-1,xhat(hkInd)));

  SigmaX = SigmaX - L*SigmaY*L';
  [~,S,V] = svd(SigmaX);
  HH = V*S*V';
  SigmaX = (SigmaX + SigmaX' + HH + HH')/4; 
  
  % Q-bump code
  if r^2>4*SigmaY % bad voltage estimate by 2-SigmaX, bump Q 
    fprintf('Bumping sigmax\n');
    SigmaX(zkInd,zkInd) = SigmaX(zkInd,zkInd)*spkfData.Qbump;
  end  

  spkfData.priorI = ik;
  spkfData.SigmaX = SigmaX;
  spkfData.xhat = xhat;
  soc_hat(k) = xhat(zkInd);
  soc_bound(k) = 3*sqrt(SigmaX(zkInd,zkInd));

  % update waitbar periodically, but not too often (slow procedure)
  if mod(k,1000)==0
    waitbar(k/length(voltage),hwait);
  end
end
close(hwait)
%% Plot
figure(1)
hold on
plot(time/60,100*soc_hat,time/60,100*soc);
plot([time/60; NaN; time/60],[100*(soc_hat+soc_bound); NaN; 100*(soc_hat-soc_bound)]);
title('SOC estimation using SPKF'); 
grid on
xlabel('Time (min)'); ylabel('SOC (%)'); legend('Estimate','Truth','Bounds');
hold off

figure(2)
hold on
plot(time/60,100*(soc-soc_hat));
plot([time/60; NaN; time/60],[100*soc_bound; NaN; -100*soc_bound]);
title('SOC estimation errors using SPKF');
xlabel('Time (min)'); ylabel('SOC error (%)'); ylim([-4 4]); 
legend('Estimation error','Bounds'); 
grid on
hold off
