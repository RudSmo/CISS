%define initial launch directory
launchdir=pwd
%set a diary file for the output
diary('logfile.txt');
% Setup TT-Toobox (not necessary if globally performed when starting matlab)
cd /home/lello/MATLAB/TT-Toolbox
setup
% Setup tamen (not necessary if globally performed when starting matlab)
cd /home/lello/MATLAB/tamen
addpath(pwd)
% back to the initial directory
cd(launchdir)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameters for integration and scaling
Nsteps=10000;
nt = 6; %numer of time points for the algorithm
tol = 1e-2; %convergence of the wavefunction (1e-2 is pretty low, but not too bad)
autofs = 2.418884326505E-2; % fs to atomic unit time conversion
Tfinfs=0.10; % in fs
Tfin=Tfinfs/autofs; % Tfin in atomic units
ts = [0:Tfinfs:Nsteps*Tfinfs];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%total number of degrees of freedom
dof=104;
% Change to HEOM folder (just to read the data correctly)
% b is the Initial state of all ones
[A,b] = read_ttheom(dof,'XOP.mat','RHO.mat');
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
nedof = 16; % number of electronic states for each electronic site
nvdof = 8; % 2*number of vibrations
NB = 80; %number of bath dofs
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
Pop = cell(1,4);
%for k=1:nsites
    Op1 = zeros(2,2);
    Op2 = zeros(2,2);
    Op1(1, 1) = 1.0;
    Op2(2, 2) = 1.0;
    Op1 = tt_matrix(Op1, 1e-14);
    Op2 = tt_matrix(Op2, 1e-14);
    Op11 = mtkron(Op1,tt_eye(2),Op1,tt_eye(2));
    Op11 = tkron(Op11,tt_eye(n(NB+5:dof)));
    Op11 = tkron(tt_eye(n(1:NB)),Op11);
    Op22 = mtkron(Op2,tt_eye(2),Op2,tt_eye(2));
    Op22 = tkron(Op22,tt_eye(n(NB+5:dof)));
    Op22 = tkron(tt_eye(n(1:NB)),Op22);
    Op21 = mtkron(Op2,tt_eye(2),Op1,tt_eye(2));
    Op21 = tkron(Op21,tt_eye(n(NB+5:dof)));
    Op21 = tkron(tt_eye(n(1:NB)),Op21);
    Op12 = mtkron(Op1,tt_eye(2),Op2,tt_eye(2));
    Op12 = tkron(Op12,tt_eye(n(NB+5:dof)));
    Op12 = tkron(tt_eye(n(1:NB)),Op12);
 %   Op = tkron(Op,tt_eye(n(NB+2*k:dof)));
 %   Op = tkron(tt_eye(n(1:NB+2*(k-1))),Op);
    Pop{1} = Op11*Ix;
    Pop{2} = Op22*Ix;
    Pop{3} = Op21*Ix;
    Pop{4} = Op12*Ix;
%    Pop{k} = Op*Ix;
%end

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
curr10=zeros(1,Nsteps);
curr11=zeros(1,Nsteps);
norm=zeros(1,Nsteps);
rankx=zeros(1,Nsteps);
rankm=zeros(1,Nsteps);
%for k=1:nsys
%pop(k,1) = dot(Pop{k},b);
%end
norm(1) = 1.0;
U=tkron(b,tt_ones(nt));
rankx(1) = max(U.r);
rankm(1) = mean(U.r);
curr00(1) = real(dot(CurrOp00,b));
curr01(1) = real(dot(CurrOp01,b));
curr10(1) = real(dot(CurrOp10,b));
curr11(1) = real(dot(CurrOp11,b));
% Write data to file
dlmwrite('res.txt',[ts(1), norm(1), rankm(1), rankx(1)],'delimiter',' ','precision',16)
dlmwrite('curr00.txt',[ts(1), real(curr00(1)), imag(curr00(1))],'delimiter',' ','precision',16)
dlmwrite('curr01.txt',[ts(1), real(curr01(1)), imag(curr01(1))],'delimiter',' ','precision',16)
dlmwrite('curr10.txt',[ts(1), real(curr10(1)), imag(curr10(1))],'delimiter',' ','precision',16)
dlmwrite('curr11.txt',[ts(1), real(curr11(1)), imag(curr11(1))],'delimiter',' ','precision',16)
%dlmwrite('res.txt',[ts(1), pop(1:nsys,1)', norm(1), rankm(1), rankx(1)],'delimiter',' ','precision',16)
%dlmwrite('rho.txt',[ts(1), real(coh(1:ncoh,1)'), imag(coh(1:ncoh,1)')],'delimiter',' ','-append','precision',16)
for i=1:Nsteps-1
    fprintf('Step %g\n', i);
    [U,t] = tamen(U, A, tol, opts);
    u = extract_snapshot(U,t,1);
    dlmwrite('rank.txt',u.r','-append')
    [rankx(i+1),pos] = max(u.r);
    rankm(i+1) = mean(u.r);
%    p = pos
%   Compute current
    curr00(i+1) = real(dot(CurrOp00,u));
    curr01(i+1) = real(dot(CurrOp01,u));
    curr10(i+1) = real(dot(CurrOp10,u));
    curr11(i+1) = real(dot(CurrOp11,u));
    %for k=1:nsys
    %  pop(k,i+1) = real(dot(Pop{k},u));
    %end
    nrm = real(dot(Ix,u))
    norm(i+1) = nrm;
    %save -mat DensityMatrix.mat pop norm rankm rankx
    dlmwrite('res.txt',[ts(i+1), norm(i+1), rankm(i+1), rankx(i+1)],'delimiter',' ','-append','precision',16)
    dlmwrite('curr00.txt',[ts(i+1), real(curr00(i+1)), imag(curr00(i+1))],'delimiter',' ','-append','precision',16)
    dlmwrite('curr01.txt',[ts(i+1), real(curr01(i+1)), imag(curr01(i+1))],'delimiter',' ','-append','precision',16)
    dlmwrite('curr10.txt',[ts(i+1), real(curr10(i+1)), imag(curr10(i+1))],'delimiter',' ','-append','precision',16)
    dlmwrite('curr11.txt',[ts(i+1), real(curr11(i+1)), imag(curr11(i+1))],'delimiter',' ','-append','precision',16)
    %dlmwrite('res.txt',[ts(i+1), pop(1:nsys,i+1)', norm(i+1), rankm(i+1), rankx(i+1)],'delimiter',' ','-append','precision',16)
    %dlmwrite('rho.txt',[ts(i+1), real(coh(1:ncoh,i+1)'), imag(coh(1:ncoh,i+1)')],'delimiter',' ','-append','precision',16)
    %subplot(2,2,1); plot(ts(1:i+1),pop(:,1:i+1)'); title('Populations');
    %subplot(2,2,3); plot(ts(1:i+1),norm(1:i+1)-1); title('Norm-1');
    %subplot(2,2,4); plot(ts(1:i+1),rankx(1:i+1), ts(1:i+1),rankm(1:i+1)-1); title('Max and Mean Rank');
    %drawnow;
end
t_adapt = toc;
fprintf('CPU time in the adaptive regime is %g\n', t_adapt);
fprintf('Produced %d time intervals\n\n', numel(t));

