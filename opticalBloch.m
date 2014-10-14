% OPTICALBLOCH  Numerically integrate optical Bloch equations.
%
% A. Hunter, Wed Oct 16 12:55:01 MDT 2013
% --------------------------------------------------------------------


clear('all');
close('all');


%  Set up conditions for the simulation  -----------------------------
%  Output file
outFile = 'floppingphasepulses10T2.dat';

%  Time range over which to integrate
x = -1000:10:5000;  % fs

%  Physical constants
hbar = 660;  % meV fs

%  Constants specific to this simulation
MU    = 1;
OMEGA = 1548/hbar; % fs^-1
T1    = 5e4;       % population decay time in fs
T2    = 2000;      % dephasing time in fs

%  Initialize variables
h = x(2) - x(1);    % time step, in fs
t = x(1);           % fs
w = -1;             % population inversion (see Boyd Sec. 6.4)
p = 0;              % polarization
K = 2*MU/hbar;      % coupling constant


%  Optical Bloch equations  ------------------------------------------
%
%  Assumptions:
%
%  1. Dipole interaction.
%
%  2. Ad hoc inclusion of T1 and T2, the population decay time and
%     dephasing time.
%
%  3. Rotating wave approximation.  This is valid as long as the
%     driving light field is near resonance (for now, it is on
%     resonance), and weakly coupled, so that the Rabi frequency is
%     much smaller than the transition frequency.
%  
%  -------------------------------------------------------------------

p_t = inline('-1*p/T2-1i*hbar/4*K^2*E*w', ...
             't', 'p', 'w', 'T2', 'hbar', 'K', 'E');

w_t = inline('-1*(w+1)/T1-4/hbar*imag(E*conj(p))', ...
             't', 'p', 'w', 'T1', 'hbar', 'E');


%  With local-field correction  --------------------------------------
% $$$ p_t = inline('-1*p/T2-1i*hbar/4*K^2*(E + 20*p)*w', ...
% $$$              't', 'p', 'w', 'T2', 'hbar', 'K', 'E');
% $$$ 
% $$$ w_t = inline('-1*(w+1)/T1-4/hbar*imag((E - 20*1i*p)*conj(p))', ...
% $$$              't', 'p', 'w', 'T1', 'hbar', 'E');


%  Run simulation  ---------------------------------------------------
%  Initialize value(s) to plot
n = 1;
pulse = zeros(1, length(x));
p11 = zeros(1, length(x));
p10 = zeros(1, length(x));
p11(1) = (w+1)/2;
p10(1) = p/MU;

%  Progress bar
fprintf(1, '[                                                  ]\n ');
step = floor(length(x)/50);

%  RK4 iteration
while t < x(end)

    %  Calculate RK4 k values:
    kPOL1 = h*p_t(t,p,w,T2,hbar,K,E(t,OMEGA));
    kPOP1 = h*w_t(t,p,w,T1,hbar,E(t,OMEGA));
    
    kPOL2 = h*p_t(t+h/2,p+kPOL1/2,w+kPOP1/2,T2,hbar,K,E(t,OMEGA));
    kPOP2 = h*w_t(t+h/2,p+kPOL1/2,w+kPOP1/2,T1,hbar,E(t,OMEGA));
    
    kPOL3 = h*p_t(t+h/2,p+kPOL2/2,w+kPOP2/2,T2,hbar,K,E(t,OMEGA));
    kPOP3 = h*w_t(t+h/2,p+kPOL2/2,w+kPOP2/2,T1,hbar,E(t,OMEGA));
    
    kPOL4 = h*p_t(t+h,p+kPOL3,w+kPOP3,T2,hbar,K,E(t,OMEGA));
    kPOP4 = h*w_t(t+h,p+kPOL3,w+kPOP3,T1,hbar,E(t,OMEGA));

    %  Update t, p, and w
    t = t + h;
    p = p + 1/6*(kPOL1+2*kPOL2+2*kPOL3+kPOL4);
    w = w + 1/6*(kPOP1+2*kPOP2+2*kPOP3+kPOP4);
    
    %  Save whatever value(s) I want to plot
    pulse(n) = E(t,OMEGA)*conj(E(t,OMEGA));
    n = n+1;
    p11(n) = (w+1)/2;
    p10(n) = p/MU;
    
    %  Update progress bar
    if mod(n,step) == 0
        fprintf(1, '.');
    end

end


%  Plot results and save data  ---------------------------------------
figure(1);
plotyy(x, p11, x, pulse);
figure(2);
plot(x, real(p10), x, imag(p10));

% Ask user if results should be saved
prompt = sprintf('\nSave data to %s? Y/N [N]: ', outFile);
reply = input(prompt, 's');
if isempty(reply)
    reply = 'N';
end

% If yes, save results
if reply == 'Y' || reply == 'y'

    %  Write to file
    %  Columns 1,2,3,4,5 --> t,pulse,Re[pol],Im[pol],pop
    fid = fopen(outFile, 'w');

    fprintf(1, '[                                                  ]\n ');
    for j = 1:length(x)

        fprintf(fid, '%e %e %e %e %e\n', ...
                x(j), pulse, real(p10(j)), imag(p10(j)), p11(j));

        if mod(j,step) == 0
            fprintf(1, '.');
        end
        
    end
    
    fclose(fid);
    fprintf(1, '\nWrote %s.\n', outFile);

else
    fprintf(1, '\nDid not write data to file.\n');
end
