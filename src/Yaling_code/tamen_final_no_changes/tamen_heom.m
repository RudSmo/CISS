%define initial launch directory
launchdir=pwd
%set a diary file for the output
diary('logfile.txt');
% Setup TT-Toobox (not necessary if globally performed when starting matlab)
%cd /home/lello/MATLAB/TT-Toolbox
cd /home/fr/fr_fr/fr_rs1121/PhD/CISS/tamen/TT-Toolbox
setup
% Setup tamen (not necessary if globally performed when starting matlab)
%cd /home/lello/MATLAB/tamen
cd /home/fr/fr_fr/fr_rs1121/PhD/CISS/tamen/source/
addpath(pwd)
% back to the initial directory
cd(launchdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameters for integration and scaling
Nsteps=80000;
Nsave=2;
nt = 6; %numer of time points for the algorithm
tol = 1e-3; %convergence of the wavefunction (1e-2 is pretty low, but not too bad)
autofs = 2.418884326505E-2; % fs to atomic unit time conversion
Tfinfs=0.05; % in fs
Tfin=Tfinfs/autofs; % Tfin in atomic units
initProp = 1;
propStp = 0;
ts = [propStp*Nsteps*Tfinfs:Tfinfs:Nsteps*Tfinfs*(propStp+1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsites = 3;
nel = 2*nsites;
nedof = 2*nel; % number of electronic states for each electronic site

% Rudi: shouldn't this be 2*number of electronic states (total, all sites) because of twin-space?
% otherwise, the following will simply break.
num_vib = nsites-1;
nvdof = 2*num_vib; % 2*number of vibrations
mvib = 50;
NB = 80; %number of bath dofs

%total number of degrees of freedom
dof= NB + nvdof + nedof;
% Change to HEOM folder (just to read the data correctly)
% b is the Initial state of all ones
[A,b] = read_ttheom(dof,'XOP.mat','RHO.mat');


disp(b)

b = round(b,1e-8);
% rescale A
A = Tfin*A; 
% basis set on each dof
n = A.n;
% Read current operator
CurrOp00 = read_ttop(dof,'CurrOp00.mat')
CurrOp01 = read_ttop(dof,'CurrOp01.mat')
CurrOp10 = read_ttop(dof,'CurrOp10.mat')
CurrOp11 = read_ttop(dof,'CurrOp11.mat')
%% X and b have the following order: bath - electronic states - vibrations

% Rudi: in these Holstein model calculations, we have 10 pade poles, 2 leads, and forward and backward
% propagation = 10 * 2 * 2 = 40
% check that
if (NB ~= dof-nedof-nvdof)
	error("Please check your that the number of degrees of freedom and bath modes are consistent.")
end
%dof = length(n); %redefine total number of dof: 2*elec+vib+tilde vib+bath SD

%% Auxiliary vectors for expectation values
v1 = NB+1; v2 = dof; v3 = 0;
% Create unit-operator |1> = \sum_n | n \tilde n >
% |1> = (\sum_m1 |m1 m1>) (\sum_m2 |m2 m2>) ....
% Remeber: in this model the electronic part is described using one single operator! (and its tilde)

for i=v1:2:v2
    nl = zeros(1,n(i));
    Id = tkron(tt_tensor(nl),tt_tensor(nl)); % Init to null vector
    for k=1:n(i)
        nl(k) = 1.0+0.0j; 
        I2 = tt_tensor(nl);
        Id = Id + tkron(I2,I2); 
        nl(k) = 0.0;
    end
    Id = round(Id,1e-14); % remove null vector
    if (i == v1) 
       I1ev = Id;
    else 
       I1ev = tkron(I1ev,Id);
    end
end

% Create the vector |0_1 0_2 ...0_N> where N = number of bath spectral densities (or Matsubara expansion terms)
for i=1:NB
    nl = zeros(1,n(v3+i));
    nl(1) = 1.0+0.0j;
    v = tt_tensor(nl);
    if (i == 1) 
       v0 = v;
    else 
       v0 = tkron(v0,v);
    end
end

% This is used to enforce norm preservation
Ix = mtkron(v0,I1ev)

% Prepare population operators
Pop = cell(nel,1);
Nvib = cell(num_vib,1);
%for k=1:nsites
    Op1 = zeros(2,2);
    Op2 = zeros(2,2);
    Op1(1, 1) = 1.0;
    Op2(2, 2) = 1.0;
    Op1 = tt_matrix(Op1, 1e-14);
    Op2 = tt_matrix(Op2, 1e-14);
    Op11 = mtkron(Op1,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2));
 %   disp(Op11) 
 %   disp(tt_eye(n(NB+3:dof)))
    Op11 = tkron(Op11,tt_eye(n(NB+nedof+1:dof)));
 %   disp(Op11)
 % Rudi: Not sure if this should be n(NB+2:dof) or n(NB+1:dof). What is the most general way to write this? probably like this
 % 1:NB = BATH, NB+1:NB+1+sites = ELECTRONIC DOF, NB+2+sites:end = VIBRATION DOF
 % Could one simplify this to tkron(Op11,tt_eye(nvdof))? Probably not because one also needs the dim of each dof following on the right 
 % to the electronic dof
 % OK, so apparently, n is taking into account twin space (no surprise), and thats why we need to start NB+2*nelec+1 when contracting Op11 with the vib operators
    Op11 = tkron(tt_eye(n(1:NB)),Op11);
 %   disp(Op11)
    Op22 = mtkron(Op2,tt_eye(2),Op1,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2));   	
    Op22 = tkron(Op22,tt_eye(n(NB+nedof+1:dof)));
    Op22 = tkron(tt_eye(n(1:NB)),Op22);
    
    Op33 = mtkron(Op2,tt_eye(2),Op2,tt_eye(2),Op1,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2));
    Op33 = tkron(Op33,tt_eye(n(NB+nedof+1:dof)));
    Op33 = tkron(tt_eye(n(1:NB)),Op33);

    Op44 = mtkron(Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op1,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2));
    Op44 = tkron(Op44,tt_eye(n(NB+nedof+1:dof)));
    Op44 = tkron(tt_eye(n(1:NB)),Op44);

    Op55 = mtkron(Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op1,tt_eye(2),Op2,tt_eye(2));
    Op55 = tkron(Op55,tt_eye(n(NB+nedof+1:dof)));
    Op55 = tkron(tt_eye(n(1:NB)),Op55);

    Op66 = mtkron(Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op2,tt_eye(2),Op1,tt_eye(2));
    Op66 = tkron(Op66,tt_eye(n(NB+nedof+1:dof)));
    Op66 = tkron(tt_eye(n(1:NB)),Op66);

 %   Op21 = mtkron(Op2,tt_eye(2),Op1,tt_eye(2));
 %   Op21 = tkron(Op21,tt_eye(n(NB+5:dof)));
 %   Op21 = tkron(tt_eye(n(1:NB)),Op21);
 %   Op12 = mtkron(Op1,tt_eye(2),Op2,tt_eye(2));
 %   Op12 = tkron(Op12,tt_eye(n(NB+5:dof)));
 %   Op12 = tkron(tt_eye(n(1:NB)),Op12);
 %   Op = tkron(Op,tt_eye(n(NB+2*k:dof)));
 %   Op = tkron(tt_eye(n(1:NB+2*(k-1))),Op);
    Pop{1} = Op11*Ix;
    Pop{2} = Op22*Ix;
    Pop{3} = Op33*Ix;
    Pop{4} = Op44*Ix;
    Pop{5} = Op55*Ix;
    Pop{6} = Op66*Ix;
 %   Pop{3} = Op21*Ix;
 %   Pop{4} = Op12*Ix;
%    Pop{k} = Op*Ix;
%end
Opv = zeros(mvib,mvib);
for k=1:mvib
	Opv(k,k) = k-1
end
Opv = tt_matrix(Opv,1e-14)
Opv1 = mtkron(Opv,tt_eye(mvib),tt_eye(mvib),tt_eye(mvib));
Opv2 = mtkron(tt_eye(mvib),tt_eye(mvib),Opv,tt_eye(mvib));
Opv1 = tkron(tt_eye(n(1:NB+nedof)),Opv1)
Opv2 = tkron(tt_eye(n(1:NB+nedof)),Opv2)
Nvib{1} = Opv1*Ix;
Nvib{2} = Opv2*Ix;
%%
opts = struct;
opts.time_error_damp = 10.0; % The default value is 100. Setting 1 removes almost completely the time adaptivity.
opts.time_scheme = 'cheby';
opts.resid_damp = 10; % was 1e2 for the dimer
opts.nswp = 20;
opts.obs{1} = Ix; % norm conservation law
opts.kickrank = 4
%opts.max_full_size = 100
%opts.rmax = 200
fprintf('Running adaptive tamen...\n\n');
tic;
%pop=zeros(nsys,Nsteps);
curr00=zeros(1,Nsteps);
curr01=zeros(1,Nsteps);
pop=zeros(nel,Nsteps);
curr10=zeros(1,Nsteps);
curr11=zeros(1,Nsteps);
norm=zeros(1,Nsteps);
rankx=zeros(1,Nsteps);
rankm=zeros(1,Nsteps);
nvib=zeros(num_vib,Nsteps);
%for k=1:nsys
%pop(k,1) = dot(Pop{k},b);
%end
norm(1) = 1.0;
if (initProp == 1)
	U=tkron(b,tt_ones(nt));

	rankx(1) = max(U.r);
	rankm(1) = mean(U.r);

	% curSL (S=spin, L=lead)
	curr00(1) = real(dot(CurrOp00,b));
	curr01(1) = real(dot(CurrOp01,b));
	curr10(1) = real(dot(CurrOp10,b));
	curr11(1) = real(dot(CurrOp11,b));
	for i=1:nel
		pop(i,1) = real(dot(Pop{i},b));
	end
	for i=1:num_vib 
		nvib(i,1) = real(dot(Nvib{i},b));
	end
	% Write data to file
	dlmwrite('res.txt',[ts(1), norm(1), rankm(1), rankx(1)],'delimiter',' ','precision',16)
	dlmwrite('curr.txt',[ts(1), real(curr00(1)), real(curr10(1)), real(curr01(1)), real(curr11(1))],'delimiter',' ','precision',16)
	dlmwrite('population.txt',[ts(1), real(pop(1)), real(pop(2)), real(pop(3)), real(pop(4)), real(pop(5)), real(pop(6))],'delimiter',' ','precision',16)
	dlmwrite('nvib.txt',[ts(1), real(nvib(1)), real(nvib(2))],'delimiter',' ','precision',16)
	j = 2;
elseif (initProp == 0)
	state = load('state.xyz','-mat')
	U = state.U
	j = 1;
end
for i=1:Nsteps-1
    fprintf('Step %g\n', i);
    [U,t] = tamen(U, A, tol, opts);
    u = extract_snapshot(U,t,1);
    dlmwrite('rank.txt',u.r','-append')
    [rankx(j),pos] = max(u.r);
    rankm(j) = mean(u.r);
    curr00(j) = real(dot(CurrOp00,u));
    curr01(j) = real(dot(CurrOp01,u));
    curr10(j) = real(dot(CurrOp10,u));
    curr11(j) = real(dot(CurrOp11,u));
    for k=1:nel
    	pop(k,j) = real(dot(Pop{k},u));
    end
    for k=1:num_vib    
	nvib(k,j) = real(dot(Nvib{k},u));
    end
    nrm = real(dot(Ix,u))
    norm(j) = nrm;
    %save -mat DensityMatrix.mat pop norm rankm rankx
    dlmwrite('res.txt',[ts(j), norm(j), rankm(j), rankx(j)],'delimiter',' ','-append','precision',16)
    dlmwrite('curr.txt',[ts(j), real(curr00(j)), real(curr10(j)), real(curr01(j)), real(curr11(j))],'delimiter',' ','-append','precision',16)
    dlmwrite('population.txt',[ts(j), real(pop(1,j)), real(pop(2,j)), real(pop(3,j)), real(pop(4,j)), real(pop(5,j)), real(pop(6,j))],'delimiter',' ','-append','precision',16)
    dlmwrite('nvib.txt',[ts(j), real(nvib(1,j)), real(nvib(2,j))],'delimiter',' ','-append','precision',16)
    if i==Nsave
        save('state.xyz','U','-mat')
    end
    j = j+1; 
end
t_adapt = toc;
fprintf('CPU time in the adaptive regime is %g\n', t_adapt);
fprintf('Produced %d time intervals\n\n', numel(t));

