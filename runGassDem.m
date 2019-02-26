% ************************************************************************* 
% Script for running GassDem
% ************************************************************************* 
%
% GassDem (Gassmann Differential Effective Medium) designed for
%    an anisotropic elastic solid A with variable amounts of liquid B 
%    with a shape specified by an ellipsoid of arbitary orientation with
%    semi-axes of variable length parallel to X,Y,Z described in
%    (Mainprice, 1997) "Modelling the anisotropic seismic properties of
%    partially molten rocks found at Mid-Ocean Ridges", 
%    Tectonophysics, 279, 161-179.
% 
% GassDem program and analytical schemes are described in the paper of 
%    Kim, E., Kim, Y., and Mainprice, D. (2019), "GassDem: A MATLAB program
%    for modeling the anisotropic seismic properties of porous medium 
%    using differential effective medium theory and Gassmann?s poroelastic 
%    relationship", Computers and Geosciences.
% 
% MATLAB version LAST REVISED: February 26, 2019
%    EUNYOUNG KIM
%    Seoul National University, Seoul, Republic of Korea
%    e-mail: brilliant@snu.ac.kr
%
% ************************************************************************* 
%% Setup

% Output file name and directory
inputGassDem = [];
inputGassDem.title = 'GassDem_ModelSC3_Melt_AspectRatio5';
inputGassDem.outputDir = '/Users/ekim/Documents/GassDem/';
inputGassDem.inputDir = '/Users/ekim/Documents/GassDem/';

% Add path (input directory)
if exist(inputGassDem.inputDir,'dir') == 0
    mkdir(inputGassDem.inputDir);
    addpath(genpath(inputGassDem.inputDir));
else
    addpath(genpath(inputGassDem.inputDir));
end

% Add path (output directory)  
if exist(inputGassDem.outputDir,'dir') == 0
    mkdir(inputGassDem.outputDir);
    addpath(genpath(inputGassDem.outputDir));
else
    addpath(genpath(inputGassDem.outputDir));
end


%**************************************************************************
%                        INPUT SECTION (user defined)                     *
%**************************************************************************

%
% Fluid compressibility Bf
% Bf = 1/Kf (1/GPa), Kf: Bulk modulus
% e.g. Water  = 0.4878 1/GPa
%      Basalt = 0.0546 1/GPa at 1233C
%
% Fluid viscosity (Pa.s) 
% e.g. Water  = 1.0E-03 Pa.s 
%      Basalt = 10 Pa.s 
%
% Fluid Vp (km/s) 
% e.g. Water  = 1.480 km/s at 20C 
%      Basalt = 2.610 km/s at 1233C 
%
% Fluid density (g/cm3)
% e.g. Water  = 1.000 g/cm3
%      Basalt = 2.680 g/cm3 at 1233C
%

% errmax: Max. %error for Greens Tensor *(e.g. 0.1)*
inputGassDem.errmax = 0.1; 

% Fluid properties
inputGassDem.bf = 1/14.91; % Fluid compressibility Bf = 1/Kf (1/GPa), Kf: Bulk modulus
inputGassDem.visc = 10; % Fluid viscosity (Pa.s) 
inputGassDem.vpf = 2.35; % Fluid Vp (km/s)
inputGassDem.rhof = 2.70; % Fluid density (g/cm3)

% Preallocation
inputGassDem.mineral = cell(2,1);
inputGassDem.efile = cell(2,1);
inputGassDem.rho = zeros(2,1);
inputGassDem.ioptel = zeros(2,1);
inputGassDem.ioptoe = zeros(2,1);
inputGassDem.axis = zeros(2,3);


% Name of phase
inputGassDem.mineral(1,1) = {'ModelSC3'};
inputGassDem.mineral(2,1) = {'Basalt'};

% Elastic Constant Data (GPa) for phases 1 and 2
inputGassDem.efile(1,1) = {'Morales2018_ModelSC3_GPa.txt'};
inputGassDem.efile(2,1) = {'Basalt1200C_GPa.txt'};

% Density (g/cm3)
inputGassDem.rho(1,1) = 2.7292;
inputGassDem.rho(2,1) = 2.70;


% Ellipsoid inclusion semi-axes A1,A2,A3
% Same for all grains (or fluid) of this phase
% 1 = Sphere
% 2 = Ellipsoid
inputGassDem.ioptel(1,1) = 1; % host: sphere
inputGassDem.ioptel(2,1) = 2; % inclusion: ellipsoid

% Orientation of shape ellipsoid semi-axes');
% -I\--I/- ');
% -\I--\I/ ');
% 1 = parallel to elastic axes of each crystal');
%     ---- or  \\\\\  or IIIII or  //// etc   ');
%     ---- or  \\\\\  or IIIII or  //// etc   ');
% 2 = parallel to user-defined orientation (SPO)');
%    (type 2 for fluid inclusion)');
inputGassDem.ioptoe(1,1) = 1; % host: parallel to elastic axes of an inclusion phase
inputGassDem.ioptoe(2,1) = 2; % inclusion (fluid)

% ioptoe(MIN,1) = 1;
% Ellipsoid axes A // Elastic tensor axes X
% A1 // X1 (eg olivine = [100]: North)
% A2 // X2 (eg olivine = [010]: East)
% A3 // X3 (eg olivine = [001]: Up)
% X1 = 100 olivine, X1 = a in calcite etc.
%                                           
%                      North=0              
%                         .                 
%                        /I\                
%                       X1//A1                   
%                       *****               
%                    *    *     *                
%                 *       *        *             
%               *         *          *      
%              *          *           *     
%         X2//*           *            *    
% West=270 A2-* * * * * X3//A3 * * * * *-- East=90
%             *           *           *     
%              *          *          *     
%               *         *         *       
%                 *       *       *         
%                    *    *    *            
%                       *****               
%                         I                 
%                         I                            
%                      South=180            
%                                           
% Ellipsoid Semi-axes lengths: [A1,A2,A3]           
inputGassDem.axis(1,:) = [1,1,1]; % host: sphere 
inputGassDem.axis(2,1) = 5; % Inclusion A1
inputGassDem.axis(2,2) = 5; % Inclusion A2
inputGassDem.axis(2,3) = 1; % Inclusion A3

% User defined shape ellipsoid orientation
% * same for all grains/inclusions of this phase
% (mineral or fluid) in sample coordinates *
% A1 should be 90 degrees to A3
%                              
% al = Azimuth (degree) of ellipsoid semi-axis A1            
% xl = Inclination (degree) of ellipsoid semi-axis A1
% af = Azimuth (degree) of ellipsoid semi-axis A3
% xf = Inclination (degree) of ellipsoid semi-axis A3
inputGassDem.al = 90;
inputGassDem.xl = 0;
inputGassDem.af = 0;
inputGassDem.xf = 0;
        
% Which phase is the inclusion (1 or 2)? 
% IDEM = 2;
inputGassDem.idem = 2;

% vliq = Max. volume Fraction of fluid (0-1)
inputGassDem.vliq = 0.2;

%**************************************************************************
%                   END OF INPUT SECTION (user defined)                   *
%**************************************************************************

% Preallocation                                                      
eaxis  = zeros(1,3); 
RBACKI = zeros(3,3); 
% C      = zeros(6,6);
EC     = zeros(6,6);
Cdem1  = zeros(6,6);
Cdem2  = zeros(6,6); 
CG     = zeros(6,6); 
Cpore  = zeros(6,6); 
Csolid = zeros(6,6); 
Cdry   = zeros(6,6); 
Csat   = zeros(6,6); 
Chigh  = zeros(6,6); 
Clow   = zeros(6,6); 

NMIN    = 2; % Two phases (including host and inclusion)
mineral = cell(NMIN,1);
rho     = zeros(NMIN,1); 
axis    = zeros(NMIN,3);
ioptel  = zeros(NMIN,1);
ioptoe  = zeros(NMIN,1);
ECM     = zeros(NMIN,6,6); 
REMIN   = zeros(NMIN,3,3); 

ijkl = reshape([1,6,5,6,2,4,5,4,3],[3,3]);
l1   = reshape([1,2,3,2,3,1],[1,6]);
l2   = reshape([1,2,3,3,1,2],[1,6]);

column  = cell(1,22);
column2 = cell(1,23);
column3 = cell(1,23);

column(1,:) = {'Vb',...
               'C11','C12','C13','C14','C15','C16',...
               'C22','C23','C24','C25','C26',...
               'C33','C34','C35','C36',...
               'C44','C45','C46',... 
               'C55','C56',... 
               'C66'}; 
column2(1,:) = {'Vb','Density',...
                'VpXL','Vs1XL','Vs2XL',...
                'VpYL','Vs1YL','Vs2YL',...
                'VpZL','Vs1ZL','Vs2ZL',...
                'Q VpX','Q Vs1X','Q Vs2X',...
                'Q VpY','Q Vs1Y','Q Vs2Y',...
                'Q VpZ','Q Vs1Z','Q Vs2Z',...
                'fmax-X','fmax-Y','fmax-Z'};
column3(1,:) = {'Vb','Density',...
                'VpXH','Vs1XH','Vs2XH',...
                'VpYH','Vs1YH','Vs2YH',...
                'VpZH','Vs1ZH','Vs2ZH',...
                'Q VpX','Q Vs1X','Q Vs2X',...
                'Q VpY','Q Vs1Y','Q Vs2Y',...
                'Q VpZ','Q Vs1Z','Q Vs2Z',...
                'fmax-X','fmax-Y','fmax-Z'};


%**************************************************************************
%                           INPUT SECTION                                 *
%**************************************************************************

% Output file name
title = inputGassDem.title;

outputDir = inputGassDem.outputDir;

% Max. %error for Greens Tensor *(e.g. 0.1)*
errmax = inputGassDem.errmax; 

% Fluid compressibility Bf
% Bf = 1/Kf (1/GPa), Kf: Bulk modulus
% e.g. Water  = 0.4878 1/GPa
%      Basalt = 0.0546 1/GPa at 1233C 
bf = inputGassDem.bf;

% Fluid viscosity (Pa.s) 
% e.g. Water  = 1.0E-03 Pa.s 
%      Basalt = 10 Pa.s 
visc = inputGassDem.visc;

% Fluid Vp (km/s) 
% e.g. Water  = 1.480 km/s at 20C 
%      Basalt = 2.610 km/s at 1233C 
vpf = inputGassDem.vpf;

% Fluid density (g/cm3)
% e.g. Water  = 1.000 g/cm3
%      Basalt = 2.680 g/cm3 at 1233C
rhof = inputGassDem.rhof;

% Name of phase
mineral(:,1) = inputGassDem.mineral;

% Elastic Constant Data (GPa) for phases 1 and 2
efile(:,1) = inputGassDem.efile;

% Density (g/cm3)f
rho(:,1) = inputGassDem.rho;

% Ellipsoid inclusion semi-axes A1,A2,A3
% Same for all grains (or fluid) of this phase
% 1 = Sphere
% 2 = Ellipsoid
ioptel(:,1) = inputGassDem.ioptel;

% Orientation of shape ellipsoid semi-axes');
% -I\--I/- ');
% -\I--\I/ ');
% 1 = parallel to elastic axes of each crystal');
%     ---- or  \\\\\  or IIIII or  //// etc   ');
%     ---- or  \\\\\  or IIIII or  //// etc   ');
% 2 = parallel to user-defined orientation (SPO)');
%    (type 2 for fluid inclusion)');
ioptoe(:,1) = inputGassDem.ioptoe; 

% ioptoe(MIN,1)= 1 
% Ellipsoid axes A // Elastic tensor axes X
% A1 // X1 (eg olivine = [100]: North)
% A2 // X2 (eg olivine = [010]: East)
% A3 // X3 (eg olivine = [001]: Up)
% X1 = 100 olivine, X1 = a in calcite etc.
%                                           
%                      North=0              
%                         .                 
%                        /I\                
%                       X1//A1                   
%                       *****               
%                    *    *     *                
%                 *       *        *             
%               *         *          *      
%              *          *           *     
%         X2//*           *            *    
% West=270 A2-* * * * * X3//A3 * * * * *-- East=90
%             *           *           *     
%              *          *          *     
%               *         *         *       
%                 *       *       *         
%                    *    *    *            
%                       *****               
%                         I                 
%                         I                            
%                      South=180            
%                                           
% Ellipsoid Semi-axes lengths: [A1,A2,A3]           
axis(:,:) = inputGassDem.axis; 

% User defined shape ellipsoid orientation
% * same for all grains/inclusions of this phase
% (mineral or fluid) in sample coordinates *
% A1 should be 90 degrees to A3
%                              
% Azimuth (degree) of ellipsoid semi-axis A1            
al = inputGassDem.al;

% Inclination (degree) of ellipsoid semi-axis A1
xl = inputGassDem.xl;

% Azimuth (degree) of ellipsoid semi-axis A3
af = inputGassDem.af;

% Inclination (degree) of ellipsoid semi-axis A3
xf = inputGassDem.xf;
        
% Which phase is the inclusion (1 or 2)? 
% idem = 2;
idem = inputGassDem.idem;

% Max. volume Fraction of fluid (0-1)
vliq = inputGassDem.vliq;

%**************************************************************************
%                          END OF INPUT SECTION                           *
%**************************************************************************


% Create text files to save results
pfile  = [outputDir,deblank(title),'_GDem_Print.txt'];
fid_12 = fopen(pfile,'w+');

% Create text file to save velocities for low frequency 
cfile1 = [outputDir,deblank(title),'_GDem_VpL.txt'];
fid_16 = fopen(cfile1,'w+');
for i = 1:length(column2)
    fprintf(fid_16,'%8s \t',cell2mat(column2(i))); 
end
fprintf(fid_16,' \n');

% Create text file to save velocities for high frequency 
cfile2 = [outputDir,deblank(title),'_GDem_VpH.txt'];
fid_17 = fopen(cfile2,'w+');
for i = 1:length(column3)
    fprintf(fid_17,'%8s \t',cell2mat(column3(i))); 
end
fprintf(fid_17,' \n');

% Create text file to save Cij for low frequency 
cfile3 = [outputDir,deblank(title),'_GDem_CijL.txt'];
fid_18 = fopen(cfile3,'w+');
for i = 1:length(column)
    fprintf(fid_18,'%8s \t',cell2mat(column(i))); 
end
fprintf(fid_18,' \n');

% Create text file to save Cij for high frequency 
cfile4 = [outputDir,deblank(title),'_GDem_CijH.txt'];
fid_19 = fopen(cfile4,'w+');
for i = 1:length(column)
    fprintf(fid_19,'%8s \t',cell2mat(column(i))); 
end
fprintf(fid_19,' \n');

% Hemisphere Convention:
%    +1 = dip +VE upper hemisphere
%    -1 = dip +VE lower hemisphere
ihemi = 1;

% Hemisphere Convention For Contouring:
%    +1 = dip +VE hemisphere
%    -1 = dip -VE hemisphere
icont = 1;

fprintf(fid_12,['\n  Seismic Anisotropy Analysis', ...
                repmat(' ',1,8),'%s \n\n'],title);

%**************************************************************************
%                        MINERAL LOOP                                     *                    
%**************************************************************************

for MIN = 1:NMIN
        
    % Elastic Constant Data (GPa) for phase 1 or 2
    fid_10 = fopen(cell2mat(efile(MIN,1)),'r+');    
    
    % Read elastic constants from fid_10
    xref = cell(2,1);
    xref(1,1) = cellstr(fgetl(fid_10));
    xref(2,1) = cellstr(fgetl(fid_10));
    a_temp  = fgetl(fid_10);
    a_temp2 = textscan(a_temp,'%8.4f');
    a = a_temp2{1,1}(1:3);
    e = a_temp2{1,1}(4:6);
    indexx = a_temp2{1,1}(end);   
    
    EC = zeros(6,6);
    
    for i = 1:6        
        ec_temp  = fgetl(fid_10);
        ec_temp2 = textscan(ec_temp,'%8.4f');
        EC(i,:)  = ec_temp2{1,1}(:); % GPa         
    end
    
    ECM(MIN,:,:) = EC(:,:);
    
    fclose(fid_10);

    % Ellipsoid inclusion semi-axes A1,A2,A3
    % Same for all grains (or fluid) of this phase
    % 1 = Sphere
    % 2 = Ellipsoid
    %
    % Spherical case:
    % axes all equal one and rotation matrix equal unitary matrix
    if ioptel(MIN,1) == 1 
        axis(MIN,1) = 1;
        axis(MIN,2) = 1;
        axis(MIN,3) = 1;        
        ioptoe(MIN,1) = 1;     
        for i = 1:3
           REMIN(MIN,i,i) = 1; 
        end
                
    % Ellipsoid Semi-axes lengths: [A1,A2,A3]    
    elseif ioptel(MIN,1) == 2        
        switch ioptoe(MIN,1)
            case 1 
                for i = 1:3
                    REMIN(MIN,i,i) = 1; 
                end
            case 2
                % Orientation
                [re,az1,xinc1,az2,xinc2,az3,xinc3] = ...
                    orientation(al,xl,af,xf,ihemi);
                REMIN(MIN,:,:) = re(:,:);
        end
    end 
        
    %**********************************************************************
    %                      PRINT ELASTIC CONSTANTS                        *
    %**********************************************************************

    fprintf(fid_12,['\n  Single Crystal Elastic Stiffness Matrix \n',...
                    '    %s \n    %s \n\n'],...
                    cell2mat(xref(1,1)),cell2mat(xref(2,1)));

    for i = 1:6
        fprintf(fid_12,[repmat('%9.2f',1,6),' \n\n'],EC(i,1:6)); 
    end

    fprintf(fid_12,'  Density = %7.4f  g/cm3 \n\n',(rho(MIN,1)));
    fprintf(fid_12,['  No. of Medium: %3i   %s \n',...
                    '  Ellipsoid semi-axes lengths A1:A2:A3 = %4i:%4i:%4i \n'], ...
                    MIN,cell2mat(mineral(MIN)),axis(MIN,1),axis(MIN,2),axis(MIN,3));

    if ioptoe(MIN) == 2 
        fprintf(fid_12,['  Fixed orientation of shape ellipsoid axes: \n',...
                        '  A1 : AZ = %6.2f INC = %6.2f \n',...
                        '  A2 : AZ = %6.2f INC = %6.2f \n',...
                        '  A3 : AZ = %6.2f INC = %6.2f \n\n\n'],...
                        az1,xinc1,az2,xinc2,az3,xinc3);
    end
end

%**************************************************************************
%                       END OF MINERAL LOOP                               *
%**************************************************************************



% *************************************************************************
%  DEM LOOP
%  with TWO DEM determinations
%  a) The high frequency elatic constants 
%    (classical DEM solid-liquid)
%  b) The low frequency elatic constants
%    (classical DEM solid-empty pore plus Gassmann poro-elastic analysis)
% *************************************************************************
% 
% The self-consistent average is often close to the VRH (or Hill) average. 
% It is considered a reliable estimate of the 'true' elastic constants
%    of anisotropic rock comprising of several mineral phases
%    of contrasting elastic constants (i.e. elastically heterogeneous).
% In this program the SC polycrystalline polyphase average elastic
%    constants are used for 1OO percent solid (aka background medium)
%    at pole A of the Differential Effective Medium (DEM). 
%    Pole B is the 1OO percent fluid (aka inclusion).
% The self-consistent scheme is an interative method which requires
%    an initial value, the Voigt (upper bound) average is used as the
%    initial value of Cij and the convergence of the SC method always
%    results in a final SC value that is less than the Voigt average.
%   
%   SC = I * Csolid
%       I   *
%       I  +   *
%       I    +    *
%       I       +      *
%       I           +        *
%       I                +           *  high frequency
%       I                     +
%       I                           +  low frequency
%       I
%       I
%       ________________________________________
%       A 100% solid (background)              B 100% fluid (inclusion)
%     POLE                                   POLE
%
% High frequency (unrelaxed) elastic constants (Cij)
% (classical DEM solid-liquid)
% ____________________________________________
%           ------------------------------
%           I                  III       I
%           I   III            III       I
%           I   III                      I
%           I                            I
%           I            III             I
%           I            III             I
%           I                 III        I
%           I       III       III        I
%           I       III                  I
%           I-----------------------------
%
% The high frequency unrelaxed constants correspond to UNCONNECTED
%    fluid inclusions, hence fluid pressure is considered to be variable.
% 
% Calculated using the classical DEM equation
%
% where Cij(high frequency) = f(Csolid,Cinclusion,Vinclusion)
%
%   Csolid = Csc
%   Cinclusion = Cfluid
%   Vinclusion = volume faction of inclusion (liquid)
%
%   or 
%
%   dCdem = dVinclusion/(1-Vinclusion) (Cinclusion-Cdem) ei
%   where ei is the strain inside the inclusion
%
% Low frequency (relaxed) elastic constants (Cij)
% (classical DEM solid-empty pore plus Gassmann poro-elastic analysis)
% ____________________________________________
%           ------------------------------
%           I                  III       I
%           I   III            III       I
%           I   III\          //         I
%           I      \\        //          I
%           I       \\  /III//           I
%           I        \\//III/            I
%           I         //    \\III        I
%           I       III      \III        I
%           I       III                  I
%           I-----------------------------
%
% The high frequency unrelaxed constants correspond to CONNECTED
%    fluid inclusions, hence fluid pressure is considered to be uniform.
%
% Step 1: Calculating the Cij of Dry (empty) porous medium using
%         the classical DEM equation
%
% where Cij(dry) = f(Csolid,Cinclusion,Vinclusion)
%
% Csolid = Csc
% Cinclusion = Cpore (void)
% Vinclusion = volume faction of (empty) pores
% 
% Step 2: Calculate the Cij of Saturated (filled) porous medium using
%         the anisotropic generalization of Brown and Korringa (1975) of
%         Gassmann's (1951) poroelastic theory.
%
% where Cij (sat) = f(Csolid,Cdry,Bfluid,Vpore)
%
% Csat   = C low frequency
% Csolid = Csc
% Bfluid = fluid compressibility = 1/bulk modulus
% Vpore  = volume fraction of (filled) pores = volume fraction of inclusion
% 
% In the program 
%    the arrays correspond to the following elastic constant (Cij)
%
% Csolid(3,3) = Cij (Self-consistent average of 100 percent polycrystalline solid)
% Cpore(3,3)  = Cij (empty pore or void)
% Cdem1(3,3)  = Cij (high frequency solid/liquid mixture)
% Cdem2(3,3)  = Cij (low frequency dry porous medium)
% Csat(3,3)   = Cij (low frequency saturated porous medium)


% disp(' DEM analysis:');
% disp('              ');
% disp(' DEM: Two-component system of A and B');
% disp('  I *');
% disp('  I   *');
% disp('  I     *');
% disp('  I        *');
% disp('  I             *');   
% disp('  I                  *');   
% disp('  I___________________');
% disp('  A                   B');                 
% disp(' A: Mineral of the host (background medium)');
% disp(' B: Fluid (inclusion)');
% disp(' Inclusion is the *FLUID* phase');

% idem2 - Background medium
if idem == 1
    idem2 = 2;
elseif idem == 2
    idem2 = 1;
end

fprintf(' Phase %d is the inclusion: %s \n\n',idem,cell2mat(mineral(idem)));

fprintf(' Max. volume Fraction of fluid: %g \n\n',vliq);

% Background medium
% idem2 = 1
Cbg = squeeze(ECM(idem2,:,:)); % GPa
Csolid(:,:) = Cbg(:,:);

%**************************************************************************
%     INCLUSION (FLUID) VOLUME FRACTION LOOP
%**************************************************************************

% Stating values for Cdem1, Cpore and Cdem2
% Starting value of Cdem = Csc

% Mineral A the host (background)
Cdem1(:,:) = Cbg(:,:);
        
% Porous Mineral A the host (background)
Cdem2(:,:) = Cbg(:,:);
 
% Elastic constants for pore
Cpore(:,:) = 0;

% Starting position 99.9:0.1
v1 = 0.001;
v2 = 0.999;

% Mineral 'idem' chosen as (fluid) inclusion
MIN = idem;

% No. of integration steps
nstep = (v2-v1)*1000;

% Step size
dVb = (v2-v1)/(nstep);

% Starting value v1 will correspond to Vb in loop 13
Vb = v1;
fprintf(' dVb = %d \n',dVb);

% vliq: Max. vol. Fraction liq. (0-1), dvb: step size
nstep2 = fix(vliq/dVb); 
fprintf(' nstep2 = %d \n',nstep2);
disp(' ');

% Fluid inclusion parameters for DEM loop
% Ellipsoidal axes for (fluid) inclusion
eaxis(1,:) = axis(MIN,:); 

% Elastic constants for (fluid) inclusion
EC(:,:) = ECM(MIN,:,:);

% Rotation matrix RBACK for (fluid) Inclusion
% Background - shape parallel specimen axes
RBACK(:,:) = REMIN(MIN,:,:);

% Inverse rotation matrix for (background) matrix RBACKI
RBACKI(:,:) = RBACK';

% No rotate of (fluid) inclusion elastic constants
% as fluid is isotropic ECij = CGij
CG(:,:) = EC(:,:);

%**************************************************************************
%                   INCLUSION (FLUID) VOLUME FRACTION LOOP                *
%**************************************************************************

for  k = 1:nstep2
    if k > 1
        % Computational values of volume fraction Vb=fluid, Va=solid
        Vb = Vb + dVb;
    end
    Va = 1 - Vb;
    
    % Error analysis for Green's Tensor
    % Mineral form in background medium Cdem2 (dry) porous
    [nb,x,w,ntcheb] = error4G(Cdem2,errmax,eaxis);    
    
    % Rotate DEM (background) matrix C* 
    % (high frequency Cij of solid/liquid mixture)
    [Cdemr1] = rot(l1,l2,ijkl,RBACK,Cdem1);
    
    % Rotate DEM (background) matrix C* 
    % (dry porous medium)
    [Cdemr2] = rot(l1,l2,ijkl,RBACK,Cdem2);
    
    % Solve Eshelby inclusion problem for Green's tensor Gij (Gt1 and Gt2)
    [GtR1] = green3(Cdemr1,eaxis,x,w,ntcheb);
    [GtR2] = green3(Cdemr2,eaxis,x,w,ntcheb);
    
    % Rotate Gij to (background) matrix/specimen axes
    [Gt1] = rot(l1,l2,ijkl,RBACKI,GtR1);
    [Gt2] = rot(l1,l2,ijkl,RBACKI,GtR2);
    
    % (dC/dVb) = 1/(1-Vb)*(Ci-Cdem(Vb))*inv(I+G(Ci-Cdem(Vb))
    % high frequency solid/liquid mixture
    [Cdem1] = difdem2(Cdem1,Vb,CG,Gt1,dVb);
    
    % Low frequency dry porous solid
    [Cdem2] = difdem2(Cdem2,Vb,Cpore,Gt2,dVb);
    
    % Poro-elastic formulation
    % dry porous solid
    Cdry(:,:) = Cdem2(:,:);
    
    % Gassmann's (1951) poro-elastic formula for saturated porous medium
    vpore = Vb;    
    [Csat] = gasman(Cdry,Csolid,bf,vpore,Csat);

    % Density = Fluid(inclusion) + Rock(background) density
    drock = Vb*rho(idem) + Va*rho(idem2);
    
    % Porous dry (unfilled pores) rock density
    rhodry = Va*rho(idem2);
       
    % Cij [GPa]
    Clow(:,:)  = Csat(:,:);
    Chigh(:,:) = Cdem1(:,:);
    Cdry(:,:)  = Cdry(:,:);
    
    
    % Low frequency phase velocites in X,Y,Z directions
    % bulkl [GPa]; shearl [GPa]; vpxl [km/s]; drock [g/cm3] 
    [bulkl,shearl,vpxl,vs1xl,vs2xl,vpyl,vs1yl,vs2yl,vpzl,vs1zl,vs2zl] = ...
        xyzv(Clow,drock,ijkl);
        
    % High frequency phase velocites in X,Y,Z directions
    % bulkh [GPa]; shearh [GPa]; vpxh [km/s]; drock [g/cm3]  
    [bulkh,shearh,vpxh,vs1xh,vs2xh,vpyh,vs1yh,vs2yh,vpzh,vs1zh,vs2zh] = ...
        xyzv(Chigh,drock,ijkl);
    
    % Estimate Fmax    
    % Calclulate moduli associated with P and S waves in X,Y,Z directions
    % bulk [GPa]; shear [GPa]; cvpx [GPa]; rhodry [g/cm3]  
    [bulk,shear,cvpx,cvs1x,cvs2x,cvpy,cvs1y,cvs2y,cvpz,cvs1z,cvs2z] = ...
        xyzv3(Cdry,rhodry,ijkl);

    % Calclulate fmax associated with P waves in X,Y,Z directions
    % cvpx [GPa]; visc [Pa.s]; vpf [km/s]; rhof [g/cm3]; fvpx [1/s] 
    [fvpx,fvpy,fvpz] = fmax3(cvpx,cvpy,cvpz,visc,vpf,rhof);
        
    % 1/Qs and 1/Qp max
    % drock [g/cm3]; vpxl [km/s]; vpxh [km/s]; 
    [qpx,qs1x,qs2x] = qinv(drock,vpxl,vs1xl,vs2xl,vpxh,vs1xh,vs2xh);
    [qpy,qs1y,qs2y] = qinv(drock,vpyl,vs1yl,vs2yl,vpyh,vs1yh,vs2yh);
    [qpz,qs1z,qs2z] = qinv(drock,vpzl,vs1zl,vs2zl,vpzh,vs1zh,vs2zh);
    
    
    %**********************************************************************
    %                             PRINT SECTION                           *
    %**********************************************************************
        
    % Write to screen and file
    vbp = 100*Vb;
    fprintf(fid_16,[repmat('%12.4f \t',1,23),'\n'],...
                    vbp,drock,vpxl,vs1xl,vs2xl,...
                              vpyl,vs1yl,vs2yl,...
                              vpzl,vs1zl,vs2zl,...
                              qpx,qs1x,qs2x,...
                              qpy,qs1y,qs2y,...
                              qpz,qs1z,qs2z,...
                              fvpx,fvpy,fvpz);
    fprintf(fid_17,[repmat('%12.4f \t',1,23),'\n'],...
                    vbp,drock,vpxh,vs1xh,vs2xh,...
                              vpyh,vs1yh,vs2yh,...
                              vpzh,vs1zh,vs2zh,...
                              qpx,qs1x,qs2x,...
                              qpy,qs1y,qs2y,...
                              qpz,qs1z,qs2z,...
                              fvpx,fvpy,fvpz);
    fprintf(fid_18,[repmat('%12.4f \t',1,22),'\n'],...
                    Vb,Clow(1,1:6),Clow(2,2:6),Clow(3,3:6), ...
                       Clow(4,4:6),Clow(5,5:6),Clow(6,6));
    fprintf(fid_19,[repmat('%12.4f \t',1,22),'\n'],...
                    Vb,Chigh(1,1:6),Chigh(2,2:6),Chigh(3,3:6), ...
                       Chigh(4,4:6),Chigh(5,5:6),Chigh(6,6));
    
    % Write elastic constants and velocities Clow to fid_20 and fid_30   
    iprint = 10; % (e.g. iprint = nstep/100)

    if (fix(k/iprint)*iprint-k) == 0
    
        ivf = fix(vbp);
        
        % Convert numerical value percent volume fraction 'VbP' to text 'textVF'
        textvf = num2str(ivf);
        cijfile = [outputDir,deblank(title),'-Clow',deblank(textvf),'.txt'];
        %fprintf(1,' WRITING :  %s \n',cijfile);
        fid_20 = fopen(cijfile,'w+');
        
        % Write elastic constants Clow to fid_20
        fprintf(fid_20,'%s \n',[deblank(title),'-Clow',deblank(textvf),'.txt']);
        fprintf(fid_20, ...
            ['Cij low in GPa : Volume Fraction of liquid = ','%6.2f ',...
            'density = ','%7.4f ','g/cm3\n'],Vb,drock);

        for j = 1:3
            a(j) = 1;
            e(j) = 90;
        end

        indexx = 1;
        fprintf(fid_20,...
               [repmat('%8.4f',1,6),'%1i','\n ', ...
                repmat('%8.2f',1,6),'\n ',repmat('%8.2f',1,6),'\n ', ...
                repmat('%8.2f',1,6),'\n ',repmat('%8.2f',1,6),'\n ', ...
                repmat('%8.2f',1,6),'\n ',repmat('%8.2f',1,6) ' \n'],...
                a,e,indexx,Clow);
        
        fclose(fid_20);

        % Write seismic velocities Clow to fid_30        
        [~,~,~,~,~] = vpfile(Clow,drock,ihemi,icont,outputDir,title,textvf,ijkl);

        % 1/Qmax
        dvp = 200*(vpxh-vpxl)/(vpxh+vpxl);

        fprintf(1,['%d  Vb = %0.15g \n',...
                   '    VpXL = %0.15g   VpXH = %0.15g   dVp  = %0.15g \n',...
                   '    QpX  = %0.15g   Qs1X = %0.15g   Qs2X = %0.15g \n',...
                   '    QpY  = %0.15g   Qs1Y = %0.15g   Qs2Y = %0.15g \n',...
                   '    QpZ  = %0.15g   Qs1Z = %0.15g   Qs2Z = %0.15g \n\n'],...
                   k,Vb,vpxl,vpxh,dvp,...
                   qpx,qs1x,qs2x,...
                   qpy,qs1y,qs2y,...
                   qpz,qs1z,qs2z);
    end
    %**********************************************************************
    %                      END OF PRINT SECTION                           *
    %**********************************************************************   
end

%**************************************************************************
%            END OF INCLUSION (FLUID) VOLUME FRACTION LOOP                *
%**************************************************************************

fclose(fid_16);
fclose(fid_17);
fclose(fid_18);
fclose(fid_19);

%**************************************************************************
%                             END OF DEM LOOP                             *
%**************************************************************************






function [re,az1,xinc1,az2,xinc2,az3,xinc3] = orientation(al,xl,af,xf,ihemi)
% User defined shape ellipsoid orientation
% * same for all grains/inclusions of this phase
%   (mineral or fluid) in sample coordinates *
% A1 should be 90 degrees to A3
% 
% Inputs
%    al    - Azimuth (degree) of ellipsoid semi-axis A1
%    xl    - Inclination (degree) of ellipsoid semi-axis A1
%    af    - Azimuth (degree) of ellipsoid semi-axis A3
%    xf    - Inclination (degree) of ellipsoid semi-axis A3
%
%    ihemi - Hemisphere Convention: +1 = dip +VE upper hemisphere
%                                   -1 = dip +VE lower hemisphere
% 
% Outputs
%    re     - Rotation matrix considering shape ellipsoid and 
%                     user defined orientation (SPO) 

% X=L Z=FN  Y= Z x X
[re] = azr(al,xl,af,xf);
            
% az/inc of 3 axes of ellipsoid (rows of REij)
% convert from dir.cos to AZ/INC
disp('Shape ellipsoid rotation matrix');
disp(re);

% [az1,xinc1] = azi(re(1,1),re(2,1),re(3,1),ihemi);
% [az2,xinc2] = azi(re(1,2),re(2,2),re(3,2),ihemi);
% [az3,xinc3] = azi(re(1,3),re(2,3),re(3,3),ihemi);
[az1,xinc1] = azi(re(1,1),re(1,2),re(1,3),ihemi);
[az2,xinc2] = azi(re(2,1),re(2,2),re(2,3),ihemi);
[az3,xinc3] = azi(re(3,1),re(3,2),re(3,3),ihemi);


fprintf([' Fixed orientation of shape ellipsoid axes: \n',...
         ' A1 : AZ = %6.2f INC = %6.2f \n',...
         ' A2 : AZ = %6.2f INC = %6.2f \n',...
         ' A3 : AZ = %6.2f INC = %6.2f \n\n'],az1,xinc1,az2,xinc2,az3,xinc3);
     
end



function [r] = azr(al,xl,af,xf)
% AZ/INC Lineation & Foliation -> Rij 
%
% ellipsoid semi-axis A1
%    al: Azimuth 
%    xl: Inclination 
% ellipsoid semi-axis A3
%    af: Azimuth 
%    xf: Inclination 

r   = zeros(3,3); 
deg = 180/pi;

[x1,y1,z1] = xyz(al,xl);
[x3,y3,z3] = xyz(af,xf);

[x1,y1,z1] = dircos(x1,y1,z1);
[x3,y3,z3] = dircos(x3,y3,z3);

dot = x1*x3 + y1*y3 + z1*z3;

if dot > 1
    dot = 1; 
elseif dot < -1
    dot = -1; 
end

value = acos(dot)*deg;

fprintf(' Check 1 :- Angle A1 to A3 = %0.15g \n',value);

% repeat to round-off errors
[x2,y2,z2] = cross(x3,y3,z3,x1,y1,z1);
[x2,y2,z2] = dircos(x2,y2,z2);
[x3,y3,z3] = cross(x1,y1,z1,x2,y2,z2);
[x1,y1,z1] = cross(x2,y2,z2,x3,y3,z3);
[x2,y2,z2] = cross(x3,y3,z3,x1,y1,z1);

dot = x1*x3 + y1*y3 + z1*z3;

if dot > 1
    dot = 1; 
elseif dot < -1
    dot = -1; 
end

value = acos(dot)*deg;

fprintf(' Check 2 :- Angle A1 to A3 = %0.15g \n',value);

% stored by rows
r(1,1) = x1;
r(1,2) = y1;
r(1,3) = z1;
r(2,1) = x2;
r(2,2) = y2;
r(2,3) = z2;
r(3,1) = x3;
r(3,2) = y3;
r(3,3) = z3;

rdet = det(r);

fprintf([' Check 3 :- Det|R| should be +1 \n',...
         ' Det|R| = %0.15g \n'],rdet);

end




function [az,xinc] = azi(x,y,z,ihemi)
% Convert cartesian to geographic
%
% x,y,z --> az,xinc
% x = North y = West z = Vertical
% x,y,z dir cosines - anticlockwise +ve right handed
% xinc = 90 vertical, xinc = 0 horizontal
% convert to +ve HEMIsphere

if z < 0
    xt = -x;
    yt = -y;
    zt = -z;
else
    xt = x;
    yt = y;
    zt = z;
end

% convert to dir.cos
[xt,yt,zt] = dircos(xt,yt,zt);

% az anticlockwise +ve
[az] = arctan(xt,yt);

% convert to clockwise +ve (geographic) - left hand
az = 360 - az;

% polar angle (vertical to pole)
r = sqrt(xt^2 + yt^2 + zt^2);

if zt == 0
    polar = 0;
else
    polar = zt/r;
end

[~,xinc] = arccos(polar);

% convert to inclination (horizontal to pole)
xinc = 90 - xinc;

% if ihemi = -1 then make a distinction between +/-ve hemispheres
if ihemi == -1 && z < 0
    az = az + 180;
    xinc = -xinc;
end
if az > 360
    az = az - 360; 
end

end




function [angle] = arctan(x,y)
% arctan(y/x) = angle 0-360 degrees anticlockwise
% avoids underflow problems with Microsoft F77 compiler

small = 0.00001;
deg = 180/pi;
ax = abs(x);
ay = abs(y);

% Exceptions 0,90,180,270
if ay < small && ax < small 
    angle = 0; 
    return;
elseif ay <= small && x > 0
    angle = 0;
    return;
elseif ax <= small && y > 0
    angle = 90; 
    return;
elseif ay <= small && x < 0
    angle = 180; 
    return;
elseif ax <= small && y < 0
    angle = 270; 
    return;
end

% General case
ratio = y/x;
if ratio > 350
    ratio = 350; 
elseif ratio < -350
    ratio = -350; 
end

angle = deg*atan(ratio);

if x > 0 && y > 0 
    return;
elseif x < 0 && y > 0
    angle = 180 + angle; 
    return;
elseif x < 0 && y < 0
    angle = 180 + angle; 
    return;
elseif x > 0 && y < 0
    angle = 360 + angle;
end

end




function [x,angle] = arccos(x)

if x > 1
    x = 1; 
elseif x < -1
    x = -1; 
end

angle = acos(x)*180/pi;

end






function [x,y,z] = xyz(az,xinc) 
% az,xinc --> x,y,z (direction cosines)
%
% Right handed system 
%       x = North y = West z = Vertical

rad  = pi/180;
raz  = az*rad;
rinc = xinc*rad;

x = cos(raz)*cos(rinc);
y = -sin(raz)*cos(rinc);
z = sin(rinc);

end



function [x,y,z] = dircos(x,y,z)
% Convert to direction cosines

small = 0.0001;
xm = sqrt(x^2 + y^2 + z^2);

if abs(x) > small
    x = x/xm; 
elseif abs(x) <= small
    x = 0; 
end

if abs(y) > small
    y = y/xm; 
elseif abs(y) <= small
    y = 0; 
end

if abs(z) > small
    z = z/xm; 
elseif abs(z) <= small
    z = 0; 
end

end



function [x3,y3,z3] = cross(x1,y1,z1,x2,y2,z2)
% Vector cross-product 1 x 2 = 3 right-handed

x3 = (y1*z2) - (z1*y2);
y3 = (z1*x2) - (x1*z2);
z3 = (x1*y2) - (y1*x2);

end



function [nb,x,w,ntcheb] = error4G(cm,errmax,eaxis)
% Error analysis of Green's Tensor for medium Cij
%
% Cm(6,6)  = Background medium
% Errmax   = Max. Error % (0.1 is a typical value) on diagonal terms
% eaxis(3) = Ellipsoid semi-axes
% nb       = Number of iterations necessary to obtain
%      error < max. error
% ** revised June 1996 **

gtold = zeros(6,6); 

% loop over possible number of integration steps
for inb = 1:100
    nb = inb*4;
    
    [ntcheb,x,w] = tcheb(nb);
    
    % Calculate Green's tensor
    [gt] = green3(cm,eaxis,x,w,ntcheb);
    
    diag = 0; 
    
    % upper diagonal of symmetric matrix
    for i = 1:6
        for j = i:6            
            % sum diagonal terms only to avoid rounding errors
            % on off-diagonal terms in isotropic case
            if j == i
                diag = diag + abs(gt(j,i) - gtold(j,i))/abs(gt(j,i)); 
            end
            gtold(j,i) = gt(j,i);
        end
    end
    
    % average difference of diagonal terms as percentage
    diagav = 100*diag/6;
    
    % quit loop if convergence
    if diagav < errmax
        [ntcheb,x,w] = tcheb(nb); 
        return;
    end
end

end



function [ntcheb,x,w] = tcheb(nb)

x  = zeros(1,nb/2);
w  = zeros(1,nb/2);

ntcheb = nb/2;
x1 = 0;
x2 = 1;

[xp,wp] = gauleg(x1,x2,nb);

for i = 1:nb/2
    x(i) = 0.5 - xp(i);
    w(i) = wp(i);
end

end




function [x,w] = gauleg(x1,x2,n)
% Given the lower and upper limits of integration X1 and X2, and given N, 
%   this routine returns 
% arrays X and W of length N, containing the abscissas and weights of the 
%   Gauss-Legendre 
% N-point quadrature formula

m = n/2;
x = zeros(1,n);
w = zeros(1,n);

xm = 0.5*(x1+x2);
xl = 0.5*(x2-x1);
xn = n;

for i = 1:m
    xi = i;
    z = cos(pi*(xi-0.25)/(xn+0.5));
    p1 = 1;
    p2 = 0;
  
    for j = 1:n
        xj = j;
        p3 = p2;
        p2 = p1;
        p1 = ((2*j-1)*z*p2-(xj-1)*p3)/xj;
    end
    pp = n*(z*p1-p2)/(z*z-1);
    z1 = z;
    z = z1-p1/pp;
  
    while abs(z-z1) > eps        
        p1 = 1;
        p2 = 0;
        for j = 1:n
            xj = j;
            p3 = p2;
            p2 = p1;
            p1 = ((2*j-1)*z*p2-(xj-1)*p3)/xj;
        end
        pp = n*(z*p1-p2)/(z*z-1);
        z1 = z;
        z = z1-p1/pp;
    end
     
    x(i) = xm - xl*z;
    x(n+1-i) = xm + xl*z;
    w(i) = 2*xl/((1-z*z)*pp*pp);
    w(n+1-i) = w(i);
  
end

end


  

function [G] = green3(C2,axis,x,w,ntcheb)
% Green's tensor for Eshelby inclusion problem
%   using extended analysis of Kinoshita and Mura (1971)
%   phys. stat. sol. (a) Vol.5 pp.759-768
%   for anisotropic inclusion and in anisotropic medium
%   following Appendix of Lebensohn & Tome (1993)
%   Acta metall. mater. Vol.41 pp.2611-2624.
%
% N.B. Green's Tensor for inclusion with axes parallel to
%   the reference frame of the elastic medium (X1,X2,X3).
% Based on routine by Ricardo Lebensohn & Gilles Canova
%   modified for single site and Voigt tensors
%
% C2(6,6) = Matrix (background) stiffness Voigt tensor
% G(6,6)  = Green's Voigt tensor
% axes(3) = Ellipsoidal inclusion semi-axes

ijkl = reshape([1,6,5,6,2,4,5,4,3],[3,3]);
ijv  = reshape([1,2,3,2,3,1,1,2,3,3,1,2],[6,2]);

G  = zeros(6,6);
g2 = zeros(6,6);
ck = zeros(3,3);
cp = zeros(3,3);
xk = zeros(1,3);

% POINTS & WEIGHTS FOR GAUSS INTEGRATION
nint = 2*ntcheb;

xp = zeros(1,nint);
wp = zeros(1,nint);

for i = 1:ntcheb
    xp(i)        = 0.5 + x(i);
    wp(i)        = w(i);
    xp(i+ntcheb) = 0.5 - x(i);
    wp(i+ntcheb) = w(i);
end

% INTEGRATION INTERVALS: [0,PI/2][PI/2,PI] BOTH FOR THETA AND PHI
pas = pi/2;

stheta = zeros(nint,2);
ctheta = zeros(nint,2);
sphi   = zeros(nint,2);
cphi   = zeros(nint,2);

for in = 1:2
    xx1 = pas*(in-1);
    for ig = 1:nint
        xx = xx1 + pas*xp(ig);        
        
        stheta(ig,in) = sin(xx);
        ctheta(ig,in) = cos(xx);
        sphi(ig,in)   = sin(xx);
        cphi(ig,in)   = cos(xx);
        
    end
end

% INITIALIZE SYMMETRIC Greens Tensor
G(:,:) = 0;

% BIG LOOP (OVER INTERVALS AND INTEGRATION POINTS)
for it = 1:2
    for ip = 1:2
        for jt = 1:nint
            for jp = 1:nint
                
                sth = stheta(jt,it);
                cth = ctheta(jt,it);
                sph = sphi(jp,ip);
                cph = cphi(jp,ip);
                
                % ALFA: VERSOR in K DIRECTION
                xk(1) = sth*cph;
                xk(2) = sth*sph;
                xk(3) = cth;

                % (k**2 tch_1 Tfourier(G))=f(alfa)
                % Christoffel matrix
                for i = 1:3
                    for k = 1:i
                        ckik = 0;
                        for j = 1:3
                            for l = 1:3                                
                                mm = ijkl(j,i);
                                nn = ijkl(l,k);
                                ckik = ckik+C2(nn,mm)*xk(j)*xk(l);                                
                            end
                        end
                        ck(k,i) = ckik;
                        ck(i,k) = ckik;
                    end
                end
                
                cp(:,:) = ck(:,:);

                [ck] = inv(cp);
                
                % Symmetric alfa*alfa*ck
                % /2 instead of /4, because of integration
                % in phi=[0,pi] instead of [0,2pi]
                for i = 1:6
                    i1 = ijv(i,1);
                    i2 = ijv(i,2);
                    for j = 1:6
                        j1 = ijv(j,1);
                        j2 = ijv(j,2);
                        g2(j,i) = ck(i1,j1)*xk(i2)*xk(j2) + ...
                                  ck(i2,j1)*xk(i1)*xk(j2) + ...
                                  ck(i1,j2)*xk(i2)*xk(j1) + ...
                                  ck(i2,j2)*xk(i1)*xk(j1);
                        g2(j,i) = 0.5*g2(j,i);
                    end
                end

                % F(th,ph)/ro**6 CALCULATION
                ro = sqrt((axis(1)*xk(1))^2 + ...
                          (axis(2)*xk(2))^2 + ...
                          (axis(3)*xk(3))^2);
                alfa = 2*ro;
                sgalfa = abs(alfa)/alfa;
                t = (-alfa^3 + 6*ro*alfa^2 - 6*alfa*ro^2)*sgalfa;
                f1 = 2*t/ro^6;

                % sum
                for i = 1:6
                    for j = 1:6
                        G(j,i) = G(j,i) + g2(j,i)*sth*f1*wp(jt)*wp(jp);
                    end
                end
                
                % END OF BIG LOOP
            end
        end
    end
end

fact = axis(1)*axis(2)*axis(3)/32/pi;

G(:,:) = fact*pas*pas*G(:,:);

end




function [c] = rot(l1,l2,ijkl,r,ec)
% Rotate elastic constants from crystal to spacial coordinates
% Cijkl = Rip*Rjq*Rkr*Rls*Cpqrs


c = zeros(6,6);
d = zeros(6,6);

d(:,:) = ec(:,:);

% Rotate elastic constants form crystal to spacial coordinates
% Cijkl = Rip*Rjq*Rkr*Rls*Cpqrs
for m = 1:6
    i = l1(m);
    j = l2(m);
    % compute lower diagonal
    for n = 1:m
        k = l1(n);
        l = l2(n);
        x = 0;
        for lp = 1:3
            y = 0;
            for lq = 1:3
                lt = ijkl(lp,lq);
                y = y + r(j,lq)*(r(k,1)*(r(l,1)*d(lt,1) + ...
                                         r(l,2)*d(lt,6) + ...
                                         r(l,3)*d(lt,5)) ...
                               + r(k,2)*(r(l,1)*d(lt,6) + ...
                                         r(l,2)*d(lt,2) + ...
                                         r(l,3)*d(lt,4)) ...
                               + r(k,3)*(r(l,1)*d(lt,5) + ...
                                         r(l,2)*d(lt,4) + ...
                                         r(l,3)*d(lt,3)));
            end
            x = x + r(i,lp)*y;
        end
        c(m,n) = x;
        % copy to upper diagonal
        c(n,m) = x;
    end
end

end




function [cdem] = difdem2(cdem,vb,cinc,gt,dvb)
% Differential Effective Medium (DEM) based on method of
% R. McLaughlin (1977) Int. J. Engng. Sci. Vol.15, pp.237-244.
% B.E. Hornby et al. (1994) Geophysics Vol.59, pp.1570-1583.
%
% Cdem(6,6) = DEM stiffness tensor at Vb
% Cinc(6,6) = Inclusion stiffness tensor
% gt(6,6)   = Green's tensor for DEM
% dVb       = Integration step

cdem = cdem';
cinc = cinc';
gt   = gt';

iden  = zeros(6,6);
am    = zeros(6,6);
cdevm = zeros(6,6);

% Iij - cartesian idenity matrix
for i = 1:6
    iden(i,i) = 1;
end

% Transform all Voigt tensors to cartesian Kelvin matrices (or EigenTensors)
[cdemm] = Voigt2Kelvin(cdem);
[cincm] = Voigt2Kelvin(cinc);
[gm]    = Voigt2Kelvin(gt);


% I+G(Ci-Cdem)
for i = 1:6
    for j = 1:6
        am(i,j) = iden(i,j);
        for k = 1:6
            am(i,j) = am(i,j) + gm(i,k)*(cincm(k,j)-cdemm(k,j));
        end
    end
end


% Am = inv(I+G(Ci-Cdem))
am = inv(am);

% Cdevm=(Ci-Cdem)*inv(I+G(Ci-Cdem))
for i = 1:6
    for j = 1:6
        cdevm(i,j) = 0;
        for k = 1:6
            cdevm(i,j) = cdevm(i,j) + (cincm(i,k)-cdemm(i,k))*am(k,j);
        end
    end
end

% Cdemm -> Cdem
[cdev] = Kelvin2Voigt(cdevm);

% Cdem(Vb+dVb) = Cdem(Vb)+dVb*f*Cdev(i,j) - simple Euler method
f = 1/(1-vb);
for i = 1:6
    for j = 1:6
        cdem(i,j) = cdem(i,j) + dvb*f*cdev(i,j);
    end
end

cdem = cdem';

end




function [Kelvin] = Voigt2Kelvin(Voigt)
% The Voigt matrix does not preserves the norm of the tensor.
% Converts Voigt to Kelvin matrix and Kelvin to Voigt matrix
% The Kelvin elastic matrix preserves the norm of the tensor
% and its eigen-values and vectors have physical meaning.
% 
% All mathematical matrix functions can be applied to Kelvin matrix 
%
% References
% W.K. Thomson (Lord Kelvin),1856 
%   Elements of a mathematical theory of elasticity,
%   Philos. Trans. R. Soc. 146, 481-498.
% W.K. Thomson (Lord Kelvin),1856 
%   On six principal strains of an elastic solid, 
%   Phil. Trans. R. Soc., 166, 495-498.
% W.K. Thomson (Lord Kelvin),1878 
%   Mathematical theory of elasticity,
%   Encycl. Br. 7, 819-825. 
%
% INPUT
% Voigt matrix
% 
% OUTPUT
% Kelvin matrix
%
% David Mainprice (18/07/2018)
%**************************************************************************
% Matrix conversion Voigt to Kelvin = inv(Kelvin2Voigt)
Voigt2KelvinMatrix = [[  1      0      0     0      0       0];  ...
                     [   0      1      0     0      0       0];  ...
                     [   0      0      1     0      0       0];  ...
                     [   0      0      0   sqrt(2)  0       0];  ...
                     [   0      0      0     0    sqrt(2)   0];  ...
                     [   0      0      0     0      0  sqrt(2)]];
% Matrix Voigt to Kelvin method  
Kelvin = Voigt2KelvinMatrix*Voigt*Voigt2KelvinMatrix;
end



function [Voigt] = Kelvin2Voigt(Kelvin)
% The Voigt matrix does not preserves the norm of the tensor.
% Converts Voigt to Kelvin matrix and Kelvin to Voigt matrix
% The Kelvin elastic matrix preserves the norm of the tensor
% and its eigen-values and vectors have physical meaning.
% 
% All mathematical matrix functions can be applied to Kelvin matrix 
%
% References
% W.K. Thomson (Lord Kelvin),1856 
%   Elements of a mathematical theory of elasticity, 
%   Philos. Trans. R. Soc. 146, 481-498.
% W.K. Thomson (Lord Kelvin),1856 
%   On six principal strains of an elastic solid, 
%   Phil. Trans. R. Soc. 166, 495-498.
% W.K. Thomson (Lord Kelvin),1878 
%   Mathematical theory of elasticity, 
%   Encycl. Br. 7, 819-825. 
%
% INPUT
% Kelvin matrix
% 
% OUTPUT
% Voigt matrix
%
% David Mainprice (18/07/2018)
%**************************************************************************
% Matrix conversion Voigt to Kelvin = inv(Kelvin2Voigt)
Kelvin2VoigtMatrix = [[  1      0      0     0      0       0];...
                     [   0      1      0     0      0       0];...
                     [   0      0      1     0      0       0];...
                     [   0      0      0   1/sqrt(2)  0     0];...
                     [   0      0      0     0    1/sqrt(2) 0];...
                     [   0      0      0     0      0  1/sqrt(2)]];
% Matrix Kelvin to Voigt method  
Voigt = Kelvin2VoigtMatrix*Kelvin*Kelvin2VoigtMatrix;
end



function [Csat] = gasman(Cdry,Cmin,Bf,pore,Csat)
% Gassmann's (1951) poroelastic relation extended to anisotropic media
% by R.Brown and J.Korringa (1975) Geophysics vol.40, pp.608-616.
%
% Cdry(6,6) = Cij of dry porous media
% Cmin(6,6) = Cij of non-porous media
% Bf(6,6)   = Beta of pore fluid (1/Kf)
% pore      = Porosity (range 0 to 1.0)
% Csat(6,6) = Cij of saturated media

Sdry   = zeros(6,6); 
Smin   = zeros(6,6);
Ssat   = zeros(6,6); 
Sdry4  = zeros(3,3,3,3);
Smin4  = zeros(3,3,3,3);
Ssat4  = zeros(3,3,3,3); 
Sdiff4 = zeros(3,3,3,3);

Sdry(:,:) = Cdry(:,:);
Smin(:,:) = Cmin(:,:);

% invert Cij to Sij
Sdry = inv(Sdry); % GPa 
Smin = inv(Smin); % GPa

% convert compilance to 4 index
% Sij -> Sijkl
[Sdry4,~,~,~] = stiffness(Sdry4,Sdry,3);
[Smin4,~,~,~] = stiffness(Smin4,Smin,3);

% beta = 1/K for dry composite and non-porous
Bdry = 0;
Bmin = 0;

for alpha = 1:3
    for beta = 1:3
        Bdry = Bdry + Sdry4(alpha,alpha,beta,beta);
        Bmin = Bmin + Smin4(alpha,alpha,beta,beta);
    end
end

% dividing constant
s = (Bdry-Bmin)+(Bf-Bmin)*pore;

% Sijkl loop
for i = 1:3
    for j = 1:3
        % sum over i,j,alpha,alpha
        Sijdry = 0;
        Sijmin = 0;
        for alpha = 1:3
            Sijdry = Sijdry + Sdry4(i,j,alpha,alpha);
            Sijmin = Sijmin + Smin4(i,j,alpha,alpha);
        end
        for k = 1:3
            for l = 1:3
                %  sum over k,l,alpha,alpha
                Skldry = 0;
                Sklmin = 0;
                for alpha = 1:3
                    Skldry = Skldry + Sdry4(k,l,alpha,alpha);
                    Sklmin = Sklmin + Smin4(k,l,alpha,alpha);
                end
                % difference Sdryijkl - Satijkl
                Sdiff4(i,j,k,l) = ((Sijdry-Sijmin)*(Skldry-Sklmin))/s;
            end
        end
    end
end

% Ssat = Sdry - Sdiff
Ssat4(:,:,:,:) = Sdry4(:,:,:,:) - Sdiff4(:,:,:,:);


% convert compilance to 2 index
% Sijkl -> Sij (Mode = 1)
[~,Ssat,~,~] = stiffness(Ssat4,Ssat,1);

% invert Sij to Cij
Csat(:,:) = inv(Ssat);

end




function [a4,b2,mode,ierr] = stiffness(a4,b2,mode)
% Converts stiffnesses (Cij) and compliances (Sij) between notations
% a4: Matrix in 4-Subscript notation
% b2: Matrix in 2-Subscript notation
% in the four lines below the following information is given
% type of tensor, number of Subscripts: in -> out
% 
% MODE = 1  Compliance  Sijkl -> Sij
% MODE = 2  Stiffness   Cijkl -> Cij
% MODE = 3  Compliance    Sij -> Sijkl
% MODE = 4  Stiffness     Cij -> Cijkl
% IERR = 16 UNDEFINED MODE
%
% ref. K.Helbig 1994 Foundations of Anisotropy for Exploration Seismics
% Pergamon pp.123-124.

ierr = 0;

% Block off with non-existing code
if mode < 1 || mode > 4
    ierr = 16;
    return;
end

% Set output to zero
if mode <= 2
    b2 = zeros(6,6);
end

if mode >= 3
    a4 = zeros(3,3,3,3);
end

% Set subscript exchange and factors
for i = 1:3
    for j = 1:3
        if i == j
            p = i;
            f1 = 1;
        else
            p = 9-i-j;
            f1 = 2;
        end
        for k = 1:3
            for l = 1:3
                if k == l
                    q = k;
                    f2 = 1;
                else
                    q = 9-k-l;
                    f2 = 2;
                end
                
                % Conversion 4 -> 2
                if mode == 1
                    b2(p,q) = a4(i,j,k,l)*f1*f2; 
                end
                if mode == 2
                    b2(p,q) = a4(i,j,k,l); 
                end
                % Conversion 2 -> 4    
                if mode == 3
                    a4(i,j,k,l) = b2(p,q)/f1/f2; 
                end
                if mode == 4
                    a4(i,j,k,l) = b2(p,q); 
                                
                end
            end
        end
    end
end
end



function [bulk,shear,VpX,Vs1X,Vs2X,VpY,Vs1Y,Vs2Y,VpZ,Vs1Z,Vs2Z] = ...
    xyzv(c,drock,ijkl)
% Phase velocities in X,Y,Z directions
%
% X = lineation (East) = 010
% Y = normal to lineation (Origin) = 001
% Z = foliation pole (North) = 100

xi = zeros(1,3);

xi(1) = 1;
xi(2) = 0;
xi(3) = 0;

[V,~] = velo2(xi,drock,c,ijkl);

VpZ  = V(1);
Vs1Z = V(2);
Vs2Z = V(3);

xi(1) = 0;
xi(2) = 1;
xi(3) = 0;

[V,~] = velo2(xi,drock,c,ijkl);

VpX  = V(1);
Vs1X = V(2);
Vs2X = V(3);

xi(1) = 0;
xi(2) = 0;
xi(3) = 1;

[V,~] = velo2(xi,drock,c,ijkl);

VpY  = V(1);
Vs1Y = V(2);
Vs2Y = V(3);

% Directionally averaged velocities
Vp = (VpX + VpY + VpZ)/3;
Vs = (Vs1X + Vs2X + Vs1Y + Vs2Y + Vs1Z + Vs2Z)/6;

% Directionally averaged elastic constants (GPa)
c11 = drock*Vp^2;
c44 = drock*Vs^2;

% Isotropic shear and bulk moduli (GPa)
shear = c44;
bulk  = c11 - 4*c44/3;

end



function [v,eigvec] = velo2(x,rho,c,ijkl)
% Phase-velocity surfaces in an anisotropic medium
% revised April 1991
%
% x(3)        - Direction of interest
% rho         - Density [g/cm3]
% v           - Phase velocities (1,2,3 = P,S,SS) [km/s]
% c           - Elastic stiffness [GPa]
% eigvec(3,3) - Eigenvectors stored by columns

c = c';
v      = zeros(1,3);
t      = zeros(3,3);
eigvec = zeros(3,3);

% Form symmetric matrix Tik = Cijkl*Xj*Xl
for i = 1:3
    for k = 1:3
        t(i,k) = 0;
        for j = 1:3
            for l = 1:3
                m = ijkl(i,j);
                n = ijkl(k,l);
                t(i,k) = t(i,k) + c(m,n)*x(j)*x(l);
            end
        end
    end
end


% Determine the eigenvalues of symmetric Tij
[ei,eval] = jacobi(t,3,3);

for i = 1:3
    for j = 1:3
        eigvec(i,j) = ei(j,i);
    end
end

% X(I) = wave normal
% Eigenvalues  = function of wave velocity
% Eigenvectors = displacement vector
% polarization plane = contains wave normal and displacement vector
for i = 1:3
    de = eval(i);      % [GPa]
    dv = sqrt(de/rho); % [km/s]
    v(i) = dv;
end

% sort velocities into ascending order
for iop = 1:2
    nop = iop + 1;
    for inop = nop:3
        if v(inop) <= v(iop)
            continue;
        end
        value   = v(iop);
        v(iop)  = v(inop);
        v(inop) = value;
        for m = 1:3
            val = eigvec(iop,m);
            eigvec(iop,m) = eigvec(inop,m);
            eigvec(inop,m) = val;
        end
    end
end

eigvec = eigvec';

end




function [V,D] = jacobi(A,N,NP)
% Computes all eigenvalues and vectors of a real symmetric matrix A(i,j)
% A(N,N), sorted in physical array A(NP,NP)
% on output the elements of A above the diagonal are destroyed
% D(NP) eigenvalues of A in its first N elements
% V(NP,NP) columns contain normalized vectors
% NROT number of Jacobi rotations

A = A'; 
nmax = 100;
D = zeros(1,NP);
V = zeros(NP,NP);
b = zeros(1,nmax);
z = zeros(1,nmax);


for ip = 1:N
    for iq = 1:N
        V(ip,iq) = 0;
    end
    V(ip,ip) = 1;
end

for ip = 1:N
    b(ip) = A(ip,ip);
    D(ip) = b(ip);
    z(ip) = 0;
end

nrot = 0;

for i = 1:50
    sm = 0;
    for ip = 1:N-1
        for iq = ip+1:N
            sm = sm + abs(A(ip,iq));
        end
    end
    if sm == 0
        return;
    end
    if i < 4 
        tresh = 0.2*sm/N^2;
    else
        tresh = 0;
    end
    for ip = 1:N-1
        for iq = ip+1:N
            G = 100*abs(A(ip,iq));
            if i > 4 && (abs(D(ip))+G == abs(D(ip))) && (abs(D(iq))+G == abs(D(iq)))
                A(ip,iq) = 0;
            elseif abs(A(ip,iq)) > tresh
                h = D(iq) - D(ip);
                if abs(h)+G == abs(h)
                    t = A(ip,iq)/h;
                else
                    theta = 0.5*h/A(ip,iq);
                    t = 1/(abs(theta) + sqrt(1+theta^2));
                    if theta < 0
                        t = -t;
                    end
                end
                c = 1/sqrt(1+t^2);
                s = t*c;
                tau = s/(1+c);
                h = t*A(ip,iq);
                z(ip) = z(ip) - h;
                z(iq) = z(iq) + h;
                D(ip) = D(ip) - h;
                D(iq) = D(iq) + h;
                A(ip,iq) = 0;
                for j = 1:ip-1
                    G = A(j,ip);
                    h = A(j,iq);
                    A(j,ip) = G - s*(h+G*tau);
                    A(j,iq) = h + s*(G-h*tau);
                end
                for j = ip+1:iq-1
                    G = A(ip,j);
                    h = A(j,iq);
                    A(ip,j) = G - s*(h+G*tau);
                    A(j,iq) = h + s*(G-h*tau);
                end
                for j = iq+1:N
                    G = A(ip,j);
                    h = A(iq,j);
                    A(ip,j) = G - s*(h+G*tau);
                    A(iq,j) = h + s*(G-h*tau);
                end
                for j = 1:N
                    G = V(j,ip);
                    h = V(j,iq);
                    V(j,ip) = G - s*(h+G*tau);
                    V(j,iq) = h + s*(G-h*tau);
                end 
                nrot = nrot + 1;
            end
        end
    end
    
    for ip = 1:N
        b(ip) = b(ip) + z(ip);
        D(ip) = b(ip);
        z(ip) = 0;
    end
end

V = V';

disp('50 iterations should never happen');
input('press any key to continue');

end




function [bulk,shear,CVpX,CVs1X,CVs2X,CVpY,CVs1Y,CVs2Y,CVpZ,CVs1Z,CVs2Z] = ...
    xyzv3(c,drock,ijkl)
% Phase velocities in X,Y,Z directions
%
% X = lineation (East) = 010
% Y = normal to lineation (Origin) = 001
% Z = foliation pole (North) = 100

xi = zeros(1,3);

xi(1) = 1;
xi(2) = 0;
xi(3) = 0;

[V,~] = velo2(xi,drock,c,ijkl);

VpZ  = V(1);
Vs1Z = V(2);
Vs2Z = V(3);

xi(1) = 0;
xi(2) = 1;
xi(3) = 0;

[V,~] = velo2(xi,drock,c,ijkl);

VpX  = V(1);
Vs1X = V(2);
Vs2X = V(3);

xi(1) = 0;
xi(2) = 0;
xi(3) = 1;

[V,~] = velo2(xi,drock,c,ijkl);

VpY  = V(1);
Vs1Y = V(2);
Vs2Y = V(3);

% Directionally averaged velocities
Vp = (VpX + VpY + VpZ)/3;
Vs = (Vs1X + Vs2X + Vs1Y + Vs2Y + Vs1Z + Vs2Z)/6;

% Directionally averaged elastic constants (GPa)
c11 = drock*Vp^2;
c44 = drock*Vs^2;

% Isotropic shear and bulk moduli (GPa)
shear = c44;
bulk  = c11 - 4*c44/3;

% Seismic elastic moduli (GPa) in directions X,Y,Z
CVpX  = drock*VpX^2;
CVs1X = drock*Vs1X^2;
CVs2X = drock*Vs2X^2;
CVpY  = drock*VpY^2;
CVs1Y = drock*Vs1Y^2;
CVs2Y = drock*Vs2Y^2;
CVpZ  = drock*VpZ^2;
CVs1Z = drock*Vs1Z^2;
CVs2Z = drock*Vs2Z^2;

end




function [fVpX,fVpY,fVpZ] = fmax3(CVpX,CVpY,CVpZ,visc,Vpf,rhof)
% CVpX,CVpY,CVpZ - Seismic elastic moduli (GPa) in directions X,Y,Z
% visc - Viscosity of fluid (Pa.s)
% Vpf  - Vp of fluid (km/s)
% rhof - Fluid density (g/cm3)
% fVpX,fVpY,fVpZ - fmax (1/s) associated with VpX,VpY,VpZ 

% Seismic Vp -> elastic moduli of fluid in GPa
CVpf = rhof*Vpf^2;

% Seismic moduli in GPa -> Pa
Xmodd = CVpX*10^9;
Ymodd = CVpY*10^9;
Zmodd = CVpZ*10^9;
Fmodd = CVpf*10^9;

% fmax
fVpX = sqrt(Fmodd*(Fmodd+Xmodd))/visc;
fVpY = sqrt(Fmodd*(Fmodd+Ymodd))/visc;
fVpZ = sqrt(Fmodd*(Fmodd+Zmodd))/visc;

end




function [Qp,Qs1,Qs2] = qinv(rho,VpL,Vs1L,Vs2L,VpH,Vs1H,Vs2H)
% rho           - Density (g/cm3)
% VpL,Vs1L,Vs2L - Low frequency velocities [km/s]
% VpH,Vs1H,Vs2H - High frequency velocities [km/s]
% Qp,Qs1,Qs2    - 1/Qmax

% Elastic constants (GPa) associated with Vp,Vs1 and Vs2
c11l  = rho*VpL^2;
c44l1 = rho*Vs1L^2;
c44l2 = rho*Vs2L^2;
c11h  = rho*VpH^2;
c44h1 = rho*Vs1H^2;
c44h2 = rho*Vs2H^2;

% 1/Qp max
Qp = (c11h-c11l)/(2*sqrt(c11h*c11l));

% 1/Qs max
Qs1 = (c44h1-c44l1)/(2*sqrt(c44h1*c44l1));
Qs2 = (c44h2-c44l2)/(2*sqrt(c44h2*c44l2));

end




function [v,eigvec,vg1,vg2,vg3] = ...
    vpfile(c,drock,ihemi,icont,outputDir,title,textvf,ijkl)

sn1 = zeros(1,3); 
sn2 = zeros(1,3); 
sn3 = zeros(1,3); 
pm1 = zeros(1,3); 
pm2 = zeros(1,3); 
pm3 = zeros(1,3); 
c4  = zeros(3,3,3,3); 

%*****************************************************************
%     Convert from Cij to Cijkl
%*****************************************************************

[c4,~,~,~] = stiffness(c4,c,4);

%*****************************************************************
%     PLANE WAVE CALCULATION OVER WHOLE GEOGRAPHIC HEMISPHERE    *
%     USING A GRID SPACING OF 6 DEGREES                          *
%*****************************************************************
% FILE=Vfile
vfile = [outputDir,deblank(title),'_Clow',deblank(textvf),'_VpG.txt'];
%fprintf(1,' WRITING :  %s \n',vfile);
fid_30 = fopen(vfile,'w+');
fprintf(fid_30,'%s \n',[deblank(title),'_Clow',deblank(textvf),'_VpG.txt']);

for i = 1:16
    xinc = icont*(90-((i-1)*6));
    for j = 1:61
        az = 90*(j-1)/15;
        
        % Convert geographic coordinates to direction cosines
        [xi] = cart(xinc,az,ihemi);
                
        % Calculate phase velocities 
        [v,eigvec] = velo2(xi,drock,c,ijkl);
        
        % Particle motion vectors for qP,qVS1,qVS2
        for k = 1:3
            pm1(k) = eigvec(1,k);
            pm2(k) = eigvec(2,k);
            pm3(k) = eigvec(3,k);
        end
        
        % Slowness vectors SN1,2,3 in s/Km
        snm1 = 1/v(1);
        snm2 = 1/v(2);
        snm3 = 1/v(3);
        for k = 1:3
            sn1(k) = snm1*xi(k);
            sn2(k) = snm2*xi(k);
            sn3(k) = snm3*xi(k);
        end
        
        % Ray-velocity vector corresponding to the slowness vectors SN1,2,3
        [vg1] = rayvel(c4,sn1,drock);
        [vg2] = rayvel(c4,sn2,drock);
        [vg3] = rayvel(c4,sn3,drock);
        
        % Magnitude of group velocity
        vg1m = sqrt(vg1(1)^2 + vg1(2)^2 + vg1(3)^2);
        vg2m = sqrt(vg2(1)^2 + vg2(2)^2 + vg2(3)^2);
        vg3m = sqrt(vg3(1)^2 + vg3(2)^2 + vg3(3)^2);
        
        % WRITE TO FILE 15
        % Phase velocities qVp,qVS1,qVS2
        fprintf(fid_30,[repmat(' ',1,1),repmat('%12.4f',1,9),' \n'],...
                        v(1),v(2),v(3));
        
        % Particle motion vectors for qVp,qVS1,qVS2
        % XP,XS1,XS2,YP,YS1,YS2,ZP,ZS1,ZS2
        fprintf(fid_30,[repmat(' ',1,1),repmat('%12.4f',1,9),' \n'],...
                        eigvec(1,1),eigvec(2,1),eigvec(3,1),...
                        eigvec(1,2),eigvec(2,2),eigvec(3,2),...
                        eigvec(1,3),eigvec(2,3),eigvec(3,3));
        
        % Group velocity vectors for qVp,qVS1,qVS2
        fprintf(fid_30,[repmat(' ',1,1),repmat('%12.4f',1,9),' \n'],...
                        vg1(1),vg1(2),vg1(3),...
                        vg2(1),vg2(2),vg2(3),...
                        vg3(1),vg3(2),vg3(3));
                 
    end
end
%******************************************************************
%     END LOOP
%******************************************************************

% CLOSE OUTPUT FILE 30 Vp file
fclose(fid_30);

end



function [x] = cart(xinc,azm,ihemi)
% Convert from spherical to cartesian coordinates
%
% NORTH X=100  WEST Y=010 UP Z=001 DOWN -Z=00-1
% IHEMI=-1 LOWER HEMISPHERE +VE DIP
% IHEMI=+1 UPPER HEMISPHERE +VE DIP

x = zeros(1,3);

rad  = pi/180;
caz  = cos(azm*rad);
saz  = sin(azm*rad);
rinc = ihemi*xinc;
cinc = cos(rinc*rad);
sinc = sin(rinc*rad);
x(1,1) = caz*cinc;
x(1,2) = -saz*cinc;
x(1,3) = sinc;

% Normalize to direction cosines
r = sqrt(x(1,1)*x(1,1) + x(1,2)*x(1,2) + x(1,3)*x(1,3));

x(:,:) = x(:,:)/r;

end



function [vg] = rayvel(c4,sn,rho)
% Calculate the ray-velocity vector corresponding to a slowness vector
%
% c4        - Stiffness tensor
% sn        - Slowness vector
% vg        - Ray velocity vector
% f(i,k)    - Components of determinant in terms of components of 
%             the slowness vector
% df(i,j,k) - Array of derivatives of f with respect to the components of 
%             the slowness vector
% dfd(i)    - Gradient of the function f

f   = zeros(3,3);
cf  = zeros(3,3); 
df  = zeros(3,3,3); 
dfd = zeros(1,3); 
vg  = zeros(1,3); 

for i = 1:3
    for k = 1:3
        f(i,k) = 0;
        for j = 1:3
            for l = 1:3
                % EQN. 4.23
                f(i,k) = f(i,k) + c4(i,j,k,l)*sn(j)*sn(l);
            end
        end
        if i == k
            f(i,k) = f(i,k) - rho; 
        end
    end
end

% Signed cofactors of f(i,k)
cf(1,1) = f(2,2)*f(3,3) - f(2,3)^2;
cf(2,2) = f(1,1)*f(3,3) - f(1,3)^2;
cf(3,3) = f(1,1)*f(2,2) - f(1,2)^2;
cf(1,2) = f(2,3)*f(3,1) - f(2,1)*f(3,3);
cf(2,1) = cf(1,2);
cf(2,3) = f(3,1)*f(1,2) - f(3,2)*f(1,1);
cf(3,2) = cf(2,3);
cf(3,1) = f(1,2)*f(2,3) - f(1,3)*f(2,2);
cf(1,3) = cf(3,1);

% Derivatives of determinant elements 
for i = 1:3
    for j = 1:3
        for k = 1:3
            df(i,j,k) = 0;
            for l = 1:3
                df(i,j,k) = df(i,j,k) + (c4(i,j,k,l) + c4(k,j,i,l))*sn(l);
            end
        end
    end
end

% Components of gradient
for k = 1:3
    dfd(1,k) = 0;
    for i = 1:3
        for j = 1:3
            dfd(1,k) = dfd(1,k) + df(i,j,k)*cf(i,j);
        end
    end
end

% Normalize to obtain group velocity vg(i)
denom = 0;
for i = 1:3
    denom = denom + sn(i)*dfd(1,i);
end

vg(:,:) = dfd(:,:)/denom;

end







