%  
%  
% [alpha_eeg,ai,varargout] = get_sim_alpha(Nsec,sr,aband,pcoef,varargin) 
%  
%  
% DESCRIPTION 
% ------------------------------------------------------------------------| 
% function to simulate an alpha dominant EEG time series
%
%  
% INPUTS 
% ------------------------------------------------------------------------| 
% Nsec: scaler length of time sample to simulate. Does not have to be an 
%       integer. Defaults to 10.
% sr: scaler sampling rate. Defaults to 1000.
% aband: two-element vector specifying the beginning and end of the
%        dominate alpha band. Defaults to [8 13].
% pcoef: scaler specifying the pink noise exponent by which 1/f noise falls
%        off (i.e. 1/f.^pcoef). Defaults to 1.
% varargin{1} As: scaler alpha scaling paramter. Determines how big the
%                 alpha peak is relative to the pink noise. Defaults to
%                 0.7.
% varargin{2} tw: tukey window width paramater. Scaler which specifies how
%                 a tukey window shapes power spectrum in the dominate
%                 alpha band. High values for this parameter result in
%                 smaller contribution to the waveform from the frequencies
%                 towards the edge of the specified band. A low value
%                 results in larger contributions from these freuqencies.
% varargin{any>2} 'plots': string tag to indicate return plots of the time
%                          series and power spectrum.
%  
%
% OUTPUTS 
% ------------------------------------------------------------------------| 
% alpha_sim: vector of simulated alpha data (1 x Nsec*sr)
% ai: alpha information structure containing fields detailing the various
%     parameter choices and indices.
% varargout{1} Pfx: power spectrum of simulated alpha.
% varargout{2} Pf: power spectrum on which the simulated alpha was modelled.
%
%  
% NOTES 
% ------------------------------------------------------------------------| 
% Examples of usage:
% [alpha_eeg,ai,Px] = get_sim_alpha(5,500,[8 12],0.9);
% [alpha_eeg,ai,Px] = get_sim_alpha(5,500,[8 12],0.9,1);
% [alpha_eeg,ai,Px] = get_sim_alpha(5,500,[8 12],0.9,1,0.5);
% [alpha_eeg,ai,Px] = get_sim_alpha(5,500,[8 12],0.9,1,0.5,'plots');
% [alpha_eeg,ai,Px] = get_sim_alpha([],[],[],[],[],[],'plots');
%
%  
% Written 04/22/2014 
% By Sam Thorpe 


function [alpha_eeg,ai,varargout] = get_sim_alpha(Nsec,sr,aband,pcoef,As,tw,varargin) 


% % input parameters
% ----------------------------------------------|
if ~nargin || isempty(Nsec), Nsec = 10; end
if nargin<2 || isempty(sr), sr = 1000; end
if nargin<3 || isempty(aband), aband = [8 13]; end
if nargin<4 || isempty(pcoef), pcoef = 1; end
if nargin<5 || isempty(As), As = 0.7; end
if nargin<6 || isempty(tw), tw = 0.75; end


% % build simulated alpha power spectrum
% ----------------------------------------------|
Nsec0 = Nsec + 2;
t0 = 0:1/sr:Nsec0;
Df0 = 1/Nsec0;
f0 = 0:Df0:sr;
L = @(t,A,K,B,M)(A + (K - A)./(1 + exp(-B*(t - M))));
Pf = 1./f0.^(pcoef);
ab = find(f0>=aband(1) & f0<=aband(2));
db = find(f0>=0 & f0<=1);
Pf(ab) = Pf(ab) + As*tukeywin(length(ab),tw).';
A =  Pf(end); 
K = Pf(db(end));
B = 10;
M = 0.4;
Pf(db) = L(f0(db),A,K,B,M);
if iseven(length(Pf))
    Pf(length(Pf)/2+1:end) = fliplr(Pf(1:length(Pf)/2));
else
    Pf((length(Pf)-1)/2+1:end) = fliplr(Pf(1:(length(Pf)-1)/2+1));
end


% % get simulated alpha time series data
% ----------------------------------------------|
tmp = rand(1,length(t0));
fcoef = sqrt(Pf).*(cos(tmp) + 1i*sin(tmp));
alpha_eeg = ifft(fcoef,'symmetric');
alpha_eeg = alpha_eeg(sr:sr+Nsec*sr);
alpha_eeg = 50*(alpha_eeg./max(abs(alpha_eeg)));
t = 0:1/sr:Nsec;
Df = 1/Nsec;
f = 0:Df:sr;
Pfx = abs(fft(alpha_eeg)).^2;


% % organize outputs
% ----------------------------------------------|
alpha_scale = As;
tukey_width = tw;
ai = structure(Nsec,sr,aband,pcoef,alpha_scale,tukey_width,t);
if nargout>2
    varargout{1} = Pfx; 
    ai.f = f;
    ai.Df = Df;
end;
if nargout>3,
    varargout{2} = Pf; 
    ai.f0 = f0;
    ai.Df0 = Df0;
end;


%                                 PLOTS 
% % ----------------------------------------------------------------------| 
%  


if any(strcmpi(varargin,'plots'))
    
    figure;
    plot(t,alpha_eeg,'linewidth',2);
    set(gca,'xlim',[0 Nsec],'fontweight','bold','fontsize',14);
    xlabel('Time (sec)','fontweight','bold','fontsize',15);
    ylabel('Amplitude','fontweight','bold','fontsize',15);
    title('Simulated alpha time series','fontweight','bold','fontsize',15);
    
    figure;
    plot(f,Pfx,'linewidth',2);
    set(gca,'xlim',[0 50],'fontweight','bold','fontsize',14);
    xlabel('Frequency (Hz)','fontweight','bold','fontsize',15);
    ylabel('Power','fontweight','bold','fontsize',15);
    title('Simulated alpha power spectrum','fontweight','bold','fontsize',15);
    
    if nargout>3,
        figure;
        plot(f0,Pf,'linewidth',2);
        set(gca,'xlim',[0 50],'fontweight','bold','fontsize',14);
        xlabel('Frequency (Hz)','fontweight','bold','fontsize',15);
        ylabel('Power','fontweight','bold','fontsize',15);
        title('Input Power Spectrum','fontweight','bold','fontsize',15);
        grid('on');
    end
end;


%                              SUBFUNCTIONS 
% % ----------------------------------------------------------------------| 
%  


%                                END ALL 
% % ----------------------------------------------------------------------| 
%  


