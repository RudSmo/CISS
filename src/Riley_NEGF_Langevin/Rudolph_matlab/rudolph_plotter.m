%% 1s MODEL CONFIGURATION TESTER
% can probably be easily modified to only test one of these. We loop over
% coordinates for each DoF to test quantities in each configuration.

%% NOTES
% We seem to have issues running simulations with less than 3 nuclei
% because the condition number on Dmat becomes so large that sqrtm(Dmat)
% cannot be computed. Not sure why yet.

%% CODE PREPARATION
clear
close all 
clc

%% LOOP OVER COORDINATES

Lvec = linspace(-6, 6, 15);
qvec = linspace(-6, 6, 15);
% Tvec = linspace(-pi/4, pi/4, 7);

% Lvec = [10];
% Tvec = [0];
% qvec = [0.0001];

% Define what coordinates we are looping over (for pre-allocation and
% indexing). Element 1 is outer loop, element 2 is inner loop.
coords = {Lvec, qvec};
strings = {'L', 'q'};


%% PRE-ALLOCATE SPACE

F_field = zeros(length(coords{1}), length(coords{2}), 3);
JL0_upstore = zeros(length(coords{1}), length(coords{2}));
JL0_downstore = zeros(length(coords{1}), length(coords{2}));
D1 = zeros(length(coords{1}), length(coords{2}));
D2 = zeros(length(coords{1}), length(coords{2}));
D3 = zeros(length(coords{1}), length(coords{2}));
Xisym1 = zeros(length(coords{1}), length(coords{2}));
Xisym2 = zeros(length(coords{1}), length(coords{2}));
Xisym3 = zeros(length(coords{1}), length(coords{2}));
Xiasym1 = zeros(length(coords{1}), length(coords{2}));
Xiasym2 = zeros(length(coords{1}), length(coords{2}));
Xiasym3 = zeros(length(coords{1}), length(coords{2}));


##########################################################
# Hamiltonian Parameters
num_laps = 1;
num_ions = 3;
N_rud = num_laps * num_ions;
ws  = 0.002*ones(N_rud);
T0  = 0.04/ws(1) ; # in eV
t0_rud  = -T0;
eu  = -6*T0;
ed  = -6*T0;
E_rud   = [eu,ed];
gPh = zeros(2);
a_rud   = 1.0;
c_rud   = 1.0;
tso = T0/40;
tPh = 0.1*T0*ones(2);
tsoPh = 0.1*tso;



## Leads
P = 0.5;
tl = T0/4;
tL = tl * (ones(2) .+ P * [1,-1])  ;
tR = tl * ones(2);
##########################################################


%% INITIATE LOOPS

count1 = 1;

for L = Lvec

    count1
    
    count2 = 1;
    
    for q = qvec
        
        %% SYSTEM PARAMETERS
        % System Space
        k_length = 0; %0.008
        k_bend = 0; %0.0001
        k_twist = 0;

        % Lead Space
        kTL = 0.00095 / ws(1); %left lead temperature
        kTR = 0.00095 / ws(1); %right lead temperature
        muL = 7.0 * T0 %1*1.2*T0; %left lead fermi level
        muR = -muL; %right lead fermi level
        c = 50*T0 %0.21; %leads energy band width

        %% INITIAL NUCLEAR GEOMETRY
        % number of sites
        N = 3;



        % Define the equilibrium geometry for our DoF (which may be different to
        % the initial conditions)
        L_eq = 0;
        q_eq = 0;
        X_eq = [L_eq; q_eq];
        
        X = [L, q];

        %% CLASSICAL DEGREES OF FREEDOM
        % Specify the number of classical DoF
        Ndof = 2;

        %% SPECIFY GAMMA MATRICES
        % Constant leads couplings, no off-diagonals
        GammaL = zeros(2 * N);
        GammaR = zeros(2 * N);
        
        gam0 = T0/4;
        
        GammaL(1, 1) = gam0 * (1 + 0.5); %first nucleus coupled to left up/down
        GammaL(2, 2) = gam0 * (1 - 0.5);
        GammaR(end - 1, end - 1) = gam0 * (1); %last nucleus coupled to right up/down
        GammaR(end, end) = gam0 * (1);
        


        %% ENERGY GRID ALLOCATION
        % For our numerical calculations, we require two separate energy grids
        % which we will use to calculate different quantities. A 'wide' grid for
        % calculating the adiabatic quantities, and a 'narrow' grid for calculating
        % higher order quantities.

        % The wide grid has a lower band-width c, the narrow grid quantities should
        % get cut off at the chemical potentials anyway
        % Our omega dimension goes into 3D
        nw = 2000; %500
        omega_grid = linspace((muL < muR) * muL + (muR < muL) * muR - 20, ...
            (muR > muL) * muR + (muL > muR) * muL +20, nw); %small grid
        omega(1,1,:) = omega_grid;

        nw1 = 5000; %5000
        omega1_grid = linspace(- c, (muL > muR) * muL + (muR > muL) * muR + 0.02, nw1); %large grid
        omega1(1,1,:) = omega1_grid;

        % Prepare rotated arrays used for integration
        omegaint(1,:,1) = omega(1,1,:);
        omega1int(1,:,1) = omega1(1,1,:); 

        %% PRELIMINARY CALCULATIONS 
        % Calculate the Fermi-Dirac distributions for the leads. 
        fL = 1 ./ (exp((omega - muL) / kTL) + 1);
        fR = 1 ./ (exp((omega - muR) / kTR) + 1);
        fL_w = 1 ./ (exp((omega1 - muL) / kTL) + 1);
        fR_w = 1 ./ (exp((omega1 - muR) / kTR) + 1);

        % Calculate the first energy derivative of the Fermi-Dirac distributions.
        fLdw = - (1 / kTL) * fL .* (1 - fL); 
        fRdw = - (1 / kTR) * fR .* (1 - fR); 
        fLdw_w = - (1 / kTL) * fL_w .* (1 - fL_w); 
        fRdw_w = - (1 / kTR) * fR_w .* (1 - fR_w); 
        

        %% CALCULATE ELECTRONIC HAMILTONIAN
        % Calculate electronic Hamiltonian and its derivatives.

        % Our default electronic Hamiltonian
        h = FranssonHamiltonian(N_rud,E_rud,t0_rud,tso,tPh,tsoPh,a_rud,c_rud,num_laps,ws,X);

        %% SELF-ENERGIES
        % Specify self-energy quantities
        ER_L = - (1i / 2) * GammaL;
        ER_R = - (1i / 2) * GammaR;
        ER = ER_L + ER_R;
        EA_L = (1i / 2) * GammaL;
        EA_R = (1i / 2) * GammaR;
        EA =  EA_L + EA_R;
        EL_L = 1i * fL .* GammaL;
        EL_R = 1i * fR .* GammaR;
        EL = EL_L + EL_R;
        EG_L = - 1i * (1 - fL) .* GammaL;                                      
        EG_R = - 1i * (1 - fR) .* GammaR;
        EG = EG_L + EG_R;

        EL_w_L = 1i * fL_w .* GammaL;
        EL_w_R = 1i * fR_w .* GammaR;
        EL_w = EL_w_L + EL_w_R;

        % First Derivatives
        ELdw_L = 1i * (fLdw .* GammaL);
        ELdw_R = 1i * (fRdw .* GammaR);
        ELdw = ELdw_L + ELdw_R;

        EGdw_L = ELdw_L;
        EGdw_R = ELdw_R;
        EGdw = EGdw_L + EGdw_R;

        %% GREEN'S FUNCTIONS
        % Specify adiabatic Green's functions. 
        GA = zeros(2*N,2*N,nw);
        GR = zeros(2*N,2*N,nw);
        GA_w = zeros(2*N,2*N,nw1);
        GR_w = zeros(2*N,2*N,nw1);

        for u = 1:nw
            GinvR = omega(u) * eye(2*N) - h - ER;
            GR(:,:,u) = inv(GinvR);
        %             cond(GinvR)

            GinvA = omega(u) * eye(2*N) - h - EA;
            GA(:,:,u) = inv(GinvA);
        end

        % We utilise the external mtimesx function for multiplying 3D matrices
        GL = pagemtimes(pagemtimes(GR, EL), GA);
        GG = pagemtimes(pagemtimes(GR, EG), GA);

        for u = 1:nw1
            GinvR = omega1(u) * eye(2*N) - h - ER;
            GR_w(:,:,u) = inv(GinvR);

            GinvA = omega1(u) * eye(2*N) - h - EA;
            GA_w(:,:,u) = inv(GinvA);
        end

        GL_w = pagemtimes(pagemtimes(GR_w, EL_w), GA_w);


        %% GREEN'S FUNCTION OMEGA DERIVATIVES
        % Define the first derivatives of adiabatic Green's functions.

        GRdw = - pagemtimes(GR, GR);
        GAdw = - pagemtimes(GA, GA);
        GLdw = - pagemtimes(GR, GL) - pagemtimes(GL, GA) + pagemtimes(pagemtimes(GR, ELdw), GA);
        GGdw = - pagemtimes(GR, GG) - pagemtimes(GG, GA) + pagemtimes(pagemtimes(GR, EGdw), GA);

        %% ADIABATIC CURRENT
        % Compute the adiabatic current characteristics
        JL0 = (1 / pi) * (pagemtimes(GL, EA_L) + pagemtimes(GR, EL_L));

        % Take a trace of each layer of the 3D matrix
        JL0uptr = 0;
        JL0downtr = 0;

        % Take trace over the up/down components individually
        for b = 1:size(JL0,2)/2

            % Up current
            JL0uptr = JL0uptr + JL0(2*b - 1, 2*b - 1, :);

            % Down current
            JL0downtr = JL0downtr + JL0(2*b, 2*b, :);

        end

        JL0uptr = real(JL0uptr);
        JL0downtr = real(JL0downtr);

        JL0upint(1,:,1) = JL0uptr(1,1,:); %rotate so we can integrate
        JL0downint(1,:,1) = JL0downtr(1,1,:);

%         plot(omegaint, JL0downint)
%         pause(0.01) 

        JL0_up = trapz(omegaint, JL0upint);
        JL0_down = trapz(omegaint, JL0downint);
        
        JL0_upstore(count1, count2) = JL0_up;
        JL0_downstore(count1, count2) = JL0_down;

        %% CALCULATE THE LANGEVIN FORCES
        % We calculate the adiabatic force on each DoF, while we calculate the
        % diffusion and viscosity for each pair of DoF.

        % Define cell array for derivative components wrt each DoF
        hdq = cell(1,Ndof);
        GL1_q = cell(1,Ndof);

        % Loop over DoF, where in a given iteration we calculate the
        % required derivatives WRT to that DoF
        for u = 1:Ndof



            %% HAMILTONIAN DERIVATIVES WRT DoF    
            % Shift the current DoF by a small value so we can numerically
            % calculate the derivative of Hsoc using forward difference
            shift = 0.0001;
            X2 = X;
            X2(u) = X2(u) + shift;
        
            % Calculate Hsoc corresponding to the shifted X2
            h2 = FranssonHamiltonian(N_rud,E_rud,t0_rud,tso,tPh,tsoPh,a_rud,c_rud,num_laps,ws,X2);

            % Use forward difference to find the derivative
            h_dq = (h - h2) / shift;
            
            % Store the total h derivative
            hdq{u} = h_dq;

            %% DEFINE GREEN'S FUNCTION DERIVATIVES
            GAdq = pagemtimes(pagemtimes(GA, hdq{u}), GA);
            GRdq = pagemtimes(pagemtimes(GR, hdq{u}), GR);
            GLdq = pagemtimes(pagemtimes(GRdq, EL), GA) + pagemtimes(pagemtimes(GR, EL), GAdq);

            GAdq_w = pagemtimes(pagemtimes(GA_w, hdq{u}), GA_w);
            GRdq_w = pagemtimes(pagemtimes(GR_w, hdq{u}), GR_w);
            GLdq_w = pagemtimes(pagemtimes(GRdq_w, EL_w), GA_w) + pagemtimes(pagemtimes(GR_w, EL_w), GAdq_w);

            GA1_q = (1 / 2i) * (pagemtimes(pagemtimes(GA, GA), pagemtimes(hdq{u}, GA))...
                - pagemtimes(pagemtimes(GA, hdq{u}), pagemtimes(GA, GA))); 
            GR1_q = (1 / 2i) * (pagemtimes(pagemtimes(GR, GR), pagemtimes(hdq{u}, GR))...
                - pagemtimes(pagemtimes(GR, hdq{u}), pagemtimes(GR, GR)));

            % We store GL1 because we need it to calculate the viscosity
            GL1_q{u} = pagemtimes(pagemtimes(GR, EL), GA1_q) ...
                - (1 / 2i) * pagemtimes(pagemtimes(GR, ELdw), GAdq) ...
                + (1 / 2i) * pagemtimes(pagemtimes(GR, hdq{u}), GLdw) ...
                + (1 / 2i) * pagemtimes(GR, GLdq);

            %% ADIABATIC FORCE
            % Calculate the x,y,z components of the adiabatic force for the
            % current nucleus

            Fq = (1i / (2 * pi)) * pagemtimes(hdq{u}, GL_w);
            Fq = real(Fq);

            % Take trace over 3D matrix
            Fqtr = 0;
            for b = 1:size(Fq,2)
                Fqtr = Fqtr + Fq(b,b,:);
            end  

            % Integrate each component
            Fint(1,:,1) = Fqtr(1,1,:); %rotate so we can integrate

%                 plot(omega1int,Fint)
%                 pause(0.01)

            % Compute the adiabatic force. 
            Fq = trapz(omega1int, Fint);     
%             Fq = 0;
            Fqvec(u) = Fq;

        end

        %% LOOP OVER DIFFUSION AND VISCOSITY TENSOR ELEMENTS
        % Now loop over each pair of DoF to calculate each element of the
        % viscosity and diffusion matrices

        for u = 1:Ndof

            for w = 1:Ndof

                %% DIFFUSION COEFFICIENT

                D = (1 / (2 * pi)) * pagemtimes(pagemtimes(hdq{u}, GG), pagemtimes(hdq{w}, GL));

                D = real(D);

                % Take trace over 3D matrix
                Dtr = 0;
                for b = 1:size(D,2)
                    Dtr = Dtr + D(b,b,:);
                end  


                Dint(1,:,1) = Dtr(1,1,:); %rotate so we can integrate

%                     plot(omegaint,Dint)
%                     pause(0.01)

                % Compute the adiabatic force. 
                D = trapz(omegaint, Dint);      
                Dmat(u,w) = D;


                %% VISCOSITY COEFFICIENT

                Xi = - (1i / (2 * pi)) * pagemtimes(hdq{u}, GL1_q{w});
                Xi = real(Xi);

                % Take trace over 3D matrix
                Xitr = 0;
                for b = 1:size(Xi,2)
                    Xitr = Xitr + Xi(b,b,:);
                end      


                Xiint(1,:,1) = Xitr(1,1,:); %rotate so we can integrate

                plot(omegaint,Xiint)
                pause(0.01)

                % Compute the adiabatic force. 
                Xi = trapz(omegaint, Xiint);      
                Ximat(u,w) = Xi;

            end

        end

        %% CALCULATE SYMMETRIC AND ANTISYMMETRIC VISCOSITY COMPONENTS
        Xisym = 0.5 * (Ximat + Ximat');
        Xiasym = 0.5 * (Ximat - Ximat');
        
        % Store symmetric viscosity components
        Xisym1(count1, count2) = Xisym(1,1);
        Xisym2(count1, count2) = Xisym(2,2);
        
        % Store asymmetric component (take +'ve value)
        Xiasym1(count1, count2) = (Xiasym(1,2));
        
        D1(count1, count2) = Dmat(1,1);
        D2(count1, count2) = Dmat(2,2);

        %% UPDATE VELOCITY HALF-STEP

        % Put adiabatic forces in the appropriate form
        Fvec = Fqvec';

        % Define a function which calculates the force on a DoF as a function
        % of the distance away from the 'equilibrium' point of that DoF
        pot_force_length = @(q) -k_length * q;
        pot_force_bend = @(q) -k_bend * q;
        pot_force_twist = @(q) -k_twist * q;

        % Calculate separation of current coordinates with equilibrium
        % coordinates
        L_sep = X(1) - L_eq;
        q_sep = X(2) - q_eq;

        % Calculate classical force due to this separation
        L_force = pot_force_length(L_sep);
        q_force = pot_force_bend(q_sep);
        pot_force = [L_force; q_force];

        % Force vector is given by potential and adiabatic force
        Ftot = pot_force + Fvec;
        
        % Store the 2D vector field for L,q,T in layers (1,2,3)
        #F_field(count1, count2) = Ftot(1:2);
        
        % Update count
        count2 = count2 + 1;
        
    end
    
    count1 = count1 + 1;
    
end
%% CALCULATIONS
Ximat

%% FIGURES

%% CONFIGURATION
    

%% FORCE VECTOR FIELD

xmat = repmat(coords{2}, length(coords{1}), 1);
ymat = repmat(transpose(coords{1}), 1, length(coords{2}));

% y-axis corresponds to coords{1} (outer loop)
% x-axis corresponds to coords{2} (inner loop)

##figure %x-axis|y-axis|x-component | y-component
##quiver(xmat, ymat, F_field(:,:,2), F_field(:,:,3))
##title('Force field')
##xlabel(strings{2})
##ylabel(strings{1})


%% EFFECTIVE BENDING POTENTIAL
% (for a given length)

figure
plot(Lvec, F_field(:,1,2))
% 
% U = -cumtrapz(qvec, F_field(1,:,2));
% figure
% plot(qvec, U)

%% ADIABATIC FORCE CURL

%% DIAGONAL DIFFUSIONS

figure
contourf(coords{2},coords{1},D1)
title('L diffusion')
colorbar

figure
contourf(coords{2},coords{1},D2)
title('q diffusion')
colorbar


%% SYMMETRIC VISCOSITIES

figure
contourf(coords{2},coords{1},Xisym1)
title('L viscosity')
colorbar

figure
contourf(coords{2},coords{1},Xisym2)
title('q viscosity')
colorbar


%% ASYMMETRIC VISCOSITIES

figure
contourf(coords{2},coords{1},Xiasym1)
title('L-q')
colorbar


%% EFFECTIVE TEMPERATURES
figure
contourf(coords{2},coords{1}, D1 ./ (2*Xisym1));
title('L eff temp')
colorbar

figure
contourf(coords{2},coords{1}, D2 ./ (2*Xisym2));
title('q eff temp')
colorbar



%% ADIABATIC CURRENT

figure
contourf(coords{2},coords{1},JL0_upstore)
title('adiabatic JL up')
colorbar

figure
contourf(coords{2},coords{1},JL0_downstore)
title('adiabatic JL down')
colorbar

