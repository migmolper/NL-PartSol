
clear 
format long

% CODED FOR: ELASTIC: 3PB
% EXPLICIT
% FRACTURE
% MESHFREE - LME+OTM 
% PLAIN STRAIN
% CONTACT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global INIT_file FB_BAR MODEL DOF FRAC

INIT_file=0;    % Start from a previous stage
PLOT_ini=1;     % Geometry at the beggining

D1=0;           % Flag for 1D
SHAPE=3;        % 1-FE lineal, 2- FE quad, 3- LME
FB_BAR=0;       % 0 - No one, 1 - BBar, 2 - FBar 3 - FBar W

%DISCRETIZATION
DISC=5;  %1- quad (1IP), 2- quad (4IP) 3- 2P1P0, 4- reverse (3), 5- 4P1P0

%PROBLEM
MODEL=1;        % 0-St. Venant;   1-Neo-Hookean
FRAC=1;         % 0-No frac;  1- Eigenerosion
DOF=1;          % Number of phases

%FORCE or DISPL
FRC=0;          % Applied ext. forces = 1;  0 if imposed displacements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nodes sp elements NNE df x_0 elem
AMP=1000;          % Amplification of the geometry
[xg]=geo_mesh(SHAPE,D1,DISC,PLOT_ini,AMP);
        

[nodes,sp]=size(x_0);
[elements,NNE]=size(elem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Area to SEARCH for the fracture
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global W1 M1 H1 H2 s1 s2 ws wi
W1 = 50;        % width of the crack area (mm)
M1 = 210;       % Middle point  (mm)
H1 = 100;       % Region Top (mm)
H2 = 30;        % Region Button (mm)
s1 = 60;        % Middle support 1 (mm)
s2 = 360;       % Middle support 2 (mm)
ws = 20;        % Wide of the support (mm)
wi = 6;        % Wide of the support (mm)

imp = 138.1;   % Different material of the hammer (mm)


W1=W1*AMP;
M1=M1*AMP;
H1=H1*AMP;
H2=H2*AMP;
s1=s1*AMP;
s2=s2*AMP;
ws=ws*AMP;
wi=wi*AMP;

imp =imp*AMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add a physical problem: FRACTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global MAT_el MAT boundary vad dad Material Mat_nds thickness tf ext_forces_s 

%--------------------------------------------------------------------------
%LOCALIZATION of the bodies
%--------------------------------------------------------------------------
Material= zeros(elements,1);
Mat_nds = zeros(nodes,1);
for nodo=1:nodes
    if x_0(nodo,2)<=H1*1.0009
        Mat_nds(nodo)=1;
    else
        Mat_nds(nodo)=2;
    end
end

for e=1:elements
    if xg(e,2)<=H1*1.0009
        Material(e)=1;
    else
        Material(e)=2;
    end
end

%--------------------------------------------------------------------------
%LOCALIZATION of the properties
%--------------------------------------------------------------------------
MAT_el= zeros(elements,1);
for e=1:elements
    if xg(e,2)<=H1*1.001
        MAT_el(e)=1;
    elseif xg(e,2)<=imp
        MAT_el(e)=2;
    else
        MAT_el(e)=3;
    end
end


%--------------------------------------------------------------------------
%PARAMETERS
%--------------------------------------------------------------------------

%%%% MAT 1 %%%%
Mat=zeros(16,1);
%Elastic
Mat(1)=31e-3;                            % E (N/um2)
Mat(2)=0.2;                             % nu
Mat(3)=2368e-24;                        % rho (Tn^3/um3)
Mat(4)=Mat(1)/2/(1+Mat(2));             % G (Pa)
Mat(5)=2*Mat(4)*Mat(2)/(1-2*Mat(2));    % Lambda (Pa)
M=Mat(1)*(1-Mat(2))/((1+Mat(2))*(1-2*Mat(2)));
Mat(6)=sqrt(M/Mat(3));                  % cp (m/s)
%Fracture
Mat(7)=0.1475e-3;                       % Gc  (N/um)
Mat(8)=1.5;                             % C_epsilon

%save
MAT(:,1)=Mat;

%%%% MAT 2 %%%%
Mat=zeros(16,1);
%Elastic
Mat(1)=211e-3;                            % E (MPa)
Mat(2)=0.3;                             % nu
Mat(3)=7800e-24;                        % rho (Tn/mm3)
Mat(4)=Mat(1)/2/(1+Mat(2));             % G (Pa)
Mat(5)=2*Mat(4)*Mat(2)/(1-2*Mat(2));    % Lambda (Pa)
M=Mat(1)*(1-Mat(2))/((1+Mat(2))*(1-2*Mat(2)));
Mat(6)=sqrt(M/Mat(3));                  % cp (m/s)
%Fracture
%Mat(7)=0.1475;                          % Gc  (N/mm)
%Mat(8)=1.5;                             % C_epsilon

%save
MAT(:,2)=Mat;

%%%% MAT 3 %%%%
Mat=zeros(16,1);
%Elastic
Mat(1)=211e-3;                            % E (MPa)
Mat(2)=0.3;                             % nu
Mat(3)=214000e-24;                        % rho (Tn/mm3)
Mat(4)=Mat(1)/2/(1+Mat(2));             % G (Pa)
Mat(5)=2*Mat(4)*Mat(2)/(1-2*Mat(2));    % Lambda (Pa)
M=Mat(1)*(1-Mat(2))/((1+Mat(2))*(1-2*Mat(2)));
Mat(6)=sqrt(M/Mat(3));                  % cp (m/s)
%Fracture
%Mat(7)=0.1475;                          % Gc  (N/mm)
%Mat(8)=1.5;                             % C_epsilon

%save
MAT(:,3)=Mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BOUNDARY CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GEOM
df=DOF*sp;
L=max(x_0(:,1));
H=max(x_0(:,2));
thickness=100e3;         % If plain strain

%--------------------------------------------------------------------------
%%FORCES && BOUNDARY
%--------------------------------------------------------------------------
boundary = zeros(nodes*df,1);
vad      = zeros(nodes*df,1);
dad      = zeros(nodes*df,1);
tf       = 0;         % Time of application of the load, s
VMAX     = -2640e3;          % mm/s

for nodo=1:nodes
    
    if x_0(nodo,2)>=H1*1.001
        vad(nodo*df)=VMAX;
    end
   
    %Inforced displacement BC
%     if  x_0(nodo,2)==min(x_0(:,2))
%          if (x_0(nodo,1)>=(s1-ws/2) && x_0(nodo,1)<=(s1+ws/2)) ||...
%                  (x_0(nodo,1)>=(s2-ws/2) && x_0(nodo,1)<=(s2+ws/2))
%              boundary(nodo*df)=1;
%              dad(nodo*df)=0;
%          end
%     end
end

if FRC
    Hf=1;                    % Hf=1, Contour Horizontal,  Hf=0, Vertical
    i=0;
    nod_f(1,1)=0;
    for nodo=1:nodes
        if Hf==1;
            if x_0(nodo,2)==H && x_0(nodo,1)>=5.0
                i=i+1;
                nod_f(i,1)=nodo;
                Hf=1;
            end
        else
            if x_0(nodo,1)==L
                i=i+1;
                nod_f(i,1)=nodo;
            end
        end
    end
    %Distributed force
    f_d = -1e6; %N
    Dir_f = 2;  %1 - x, 2 - y     Direction of the force

    [ext_forces_s,~]=ext_forces(nod_f,x_0,f_d,Dir_f,0,0,Hf);
else
    ext_forces_s      = zeros(nodes*sp,1);
end

%--------------------------------------------------------------------------
% Contact conditions
%--------------------------------------------------------------------------
% CONTOUR
 global MASTER M_LIST nCC mu
% Friction
mu=0.0;
% Initial node
x0=[196 114.1];
% Circle contours: center
xc=[210 114.1];

x0=x0*AMP;
xc=xc*AMP;
% Tolerance
tol = 0.01;
% Flag 1: 1 = vertical 2 = horizontal 3 = slope 4 = circle
% Flag 2: Slope (0 if vert or horiz)
% If circle:  2) radius  3) Angle
% ANTIHOUR!!!!!!!!
CC=[4 14*AMP pi;
    1 0 0;
    2 0 0;
    1 0 0];
 
 MASTER=2;      % Master material
 if MASTER>0
    [M_LIST,nCC]=M_location(x_0,x0,xc,CC,tol);
 end


%--------------------------------------------------------------------------
% Potential Surface
%--------------------------------------------------------------------------
% CONTOUR
global K mu_2 X_PS PS_LIST
K=8;
K=K*MAT(1,1);  % Falta *h
mu_2=0.0;
% All nodes, ordered, % ANTIHOUR!!!!!!!!
x_ps1=[50.0 -1000;
      50.0 0.0;
      70.0 0.0;
      70.0 -1000]*AMP;
x_ps2=[350.0 -1000;
      350.0 0.0;
      370.0 0.0;
      370.0 -1000]*AMP; 
[nPS1,~]=size(x_ps1);
[nPS2,~]=size(x_ps2);
nPS1=nPS1-1;
nPS2=nPS2-1;
    for i=1:nPS1
        list(1)=i;
        list(2)=i+1;
        PS_LIST1{i}=list;
        clear list
    end
    for i=1:nPS2
        list(1)=i;
        list(2)=i+1;
        PS_LIST2{i}=list;
        clear list
    end
    
 X_PS{1}=x_ps1;
 X_PS{2}=x_ps2;
 PS_LIST{1}=PS_LIST1;
 PS_LIST{2}=PS_LIST2;  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SOLVER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% Solver variables: Time Integration Scheme
%--------------------------------------------------------------------------

global Time_final time_step SAVE_I SAVE_F IMPLICIT

Time_final= 6.0e-4;   %seconds
time_step = 1.5e-7;   %seconds

IMPLICIT=0;
NR = 10;      % Update K in NR every xx iterations;

     %IMPLICIT
%1-Newmark 1st
%2-Newmark 2nd
%3-Generalized-alpha
%4-HHT
%5-Wilson
%6-WBZ
%7-Collocation method
     %EXPLICIT
%1-Newmark central difference
TIS=1;
DAMP=0;

SAVE_I=50;     % Save info each XX steps
SAVE_F=10;     % Save the file each XX steps

TIME_INT_var(IMPLICIT,TIS,DAMP);

if IMPLICIT==0
	EXPL_frac_2(xg);
else
	Implicit_def(xg,NR);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




