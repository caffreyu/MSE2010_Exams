%% MSE 2010 Midterm 1

% Disclaimer: 
% This MATLAB file is purely for the review of the MSE 2010Exam 1.
% It may serve as an equation sheet, but please use it at your own risk. 

clc; clear;
close all;
format short;

%% Chapter 2. 
% Material covered: IC character; Force of attraction; Fraction of occurency. 

clc; clear

% Percent Ionic Character
% *** INPUT ***
xa = 1; % Electonegativity of the elements
xb = 1;

% *** OUTPUT ***
ic_char = (1 - exp((-0.25) * (xa - xb)^2))* 100;  % In '%'
fprintf('The ionic character is %.4f%%. \n', ic_char)

% Force of attraction
% *** INPUT ***
xa = 1;     % cation valence
xb = -1;    % anion valence
rnm = 1;    % radius in nm

r = rnm*1e-9;   % radius in m
const = 2.31e-28;   % constant for force, in N/m^2

% *** OUTPUT ***
Fa = const * abs(xa) * abs(xb) / r^2;  % Attractive force in N
fprintf('Force of attraction is %.4f N. \n', Fa)

% Fraction of occurency
% *** INPUT ***
wt_a = 78.918; % in amu
wt_b = 80.916;
wt_avg = 79.903; % Average wt in amu

syms x
eqn = x * wt_a + (1 - x) * wt_b - wt_avg;

% *** OUTPUT ***
S = double(solve(eqn == 0, x));    % Fraction of occurency for a
1 - S;   % Fraction of occurency for b
fprintf('Fraction of occurency for a is %.4f. \n', S)
fprintf('Fraction of occurency for b is %.4f. \n', 1-S)

%% Chapter 3
% Material covered: Theoretical density; Ceramic Theoretical Density

clc; clear;

% Theoretical Density
% Specify the crystal structure at the very beginning
% *** INPUT ***
BCC = true;
FCC = false;
simple = false;
r = 0.1363e-7;  % radius of the atom, cm
atomic_weight = 95.94;  % g/mol

n = getN(BCC, FCC, simple);     % atoms/unit cell
NA = 6.022e23;
l = getLength(BCC, FCC, simple, r);  % length of the unit cell, cm
Vc = l^3;   % volume of the cubic, cm^3

% *** OUTPUT ***
rho = n * atomic_weight / Vc / NA;  % g/cm^3, pay attention to the unit
fprintf('Theoretical density is %.4f g/cm^3. \n', rho)

% Get the radius of an atom based on density
% *** INPUT ***
% Specify the crystal structure at the very beginning
BCC = false;
FCC = true;
simple = false;
rho = 22.4;    % g/cm^3
atomic_weight = 192.2;  % g/mol

n = getN(BCC, FCC, simple);
syms R
l = getLength(BCC, FCC, simple, R);
eqn = n * atomic_weight / l^3 / NA;
S = double(solve(eqn == rho));
S_vector = real(S);
fprintf('The radius for the atom is %.12f cm. \n', S_vector(1))

% Ceramic theoretical density
% *** INPUT ***
% Specify the crystal structure at the very beginning
BCC = true;
FCC = false;
simple = false;
anion_radius = 0.181;   % Ionic radius, in nm
cation_radius = 0.170;  % Cation radius, in nm
anion_atomic_weight = 132.91;    % Anion atomic weight, in g/mol
cation_atomic_weight = 35.45;   % Cation atomic weight, in g/mol

n = getN(BCC, FCC, simple) / 2;
V = getCeramicV(BCC, FCC, simple, anion_radius, cation_radius);

% *** OUTPUT ***
rho_here = n * (anion_atomic_weight + cation_atomic_weight) / V / NA;
fprintf('The ceramic theoretical density is %.4f g/cm^3. \n', rho_here)


%% Chapter 4
% Material covered: Repeat unit molecular weight; Polymer molecular weight
% (number-average/weight-average); Degree of polymerization; Chain length
% and End-toend distance; Percent crystallinity; Copolymers

clc; clear;

% Repeat unit molecular weight
% *** INPUT ***
% Change the numbers below
Unit = 'C12H4N2O3';

% *** OUTPUT ***
weight = calculateWeight(Unit)
fprintf('The molecular weight of this repeat unit is %.4f g/mol. \n', weight);

% Polymer molecular weight
% *** INPUT ***
weight_range = [5000 15000; 15000 25000; 25000 35000];  % weight range in g/mol
xi = [0.50 0.40 0.10];  % for number average
wi = [0.31 0.50 0.19];  % for weight average

res = checkOne(xi, wi);
weight_avg = mean(weight_range');

% *** OUTPUT ***
number_avg_w = weight_avg * xi'
weight_avg_w = weight_avg * wi'
fprintf('The number-average molecular weight is %.4f g/mol. \n', number_avg_w);
fprintf('The weight_average molecular weight is %.4f g/mol. \n', weight_avg_w);

% Degree of polymerization
% *** INPUT ***
weight_number_avg = 1e6;    % Number-average molecular weight, g/mol
Unit_here = 'C2H4O1';

weight_here = calculateWeight(Unit_here);
DP = weight_number_avg / weight_here
fprintf('The degree of polymerization is %.1f. \n', DP)

% Chain length and End-to-end distance
% *** INPUT ***
DP = 1000;  % Degree of polymerization
N = 2;      % Number of bonds on backbone
d = 0.154;      % Length of a bond, in nm
theta = 109;    % degree

% *** OUTPUT ***
chain_len = N * DP * d * sind(theta / 2)
end_to_end_d = d * sqrt(N * DP)
fprintf('The chain length is %.4f nm. \n', chain_len)
fprintf('The end_to_end distance is %.4f nm. \n', end_to_end_d)

% Percent crystallinity
per_crys = 75;  % In %, for example, if it is 20%, just put 20
rho_sample = 1.21;    % Density of sample, in g/mol
rho_amo = 1.12;     % Density of amorphous, in g/mol

syms rho_c
eqn = rho_c * (rho_sample - rho_amo) / rho_sample / (rho_c - rho_amo) * 100;
S = double(solve(eqn == per_crys))
fprintf('The density of 100%% crystalline is %.4f g/mol. \n', S)

% Copolymers
% *** INPUT ***
Unit1 = 'C2H4';
Unit2 = 'C2H4O1';
percent_1 = 40;     % Percentage of unit1
percent_2 = 60;     % Percentage of unit2
total_weight = 32e4;    % Total weight of this copolymer, g/mol

weight_1 = calculateWeight(Unit1);
weight_2 = calculateWeight(Unit2);

% *** OUTPUT ***
m_avg = weight_1 * percent_1 / 100 + weight_2 * percent_2 / 100   % in g/mol
DP_here = total_weight / m_avg      % Degree of polymerization of this copolymer
fprintf('The average molecular weight of this copolymer is %.4f g/mol. \n', m_avg)
fprintf('The degree of polymerization of this copolymer is %.1f. \n', DP_here)

%% Chapter 5
% Material covered: Lattice sites; Vacancy concentration; Atomic % to
% weight %

clc; clear;
format short;

% Lattice sites
% *** INPUT ***
NA = 6.022e23;      % Constant, mol^-1
rho = 20.7;         % Density, g/cm^3
atomic_weight = 100.0;  % Atomic weight, g/mol

% *** OUTPUT ***
N = NA * rho / atomic_weight    % Lattice sites, atoms/m^3
fprintf('The lattice sites is %.4f per vol. \n', N)

% Vacancy concentration
% *** INPUT ***
N = 1.25e29;    % Mute this line if N is calculated from the above part
Qv = 1.0;       % Energy for vacancy formation, in eV
T = 300;        % Temp, in K
K = 8.62e-5;    % Constant, eV/K

% *** OUTPUT ***
Nv = N * exp(-Qv / K / T)       % Vacancy concentration, vacancies/m^3
fprintf('The vacancy concentration is %.4f per vol. \n', Nv)
ratio = Nv / N
fprintf('The ratio is %.20f. \n', ratio)

% Atomic % to weight %
m1 = 100;   % Mass of component 1, in g
atomic_weight_1 = 50.0;     % Atomic weight of component 1, in g/mol
m2 = 100;   % Mass of component 2, in g
atomic_weight_2 = 100.0;    % Atomic weight of component 2, in g/mol

n1 = m1 / atomic_weight_1;
n2 = m2 / atomic_weight_2;
C1 = m1 / (m1 + m2) * 100   % Weight percent of component 1
C2 = 100 - C1               % Weight percent of component 2
C11 = n1 / (n1 + n2) * 100  % Atomic percent of component 1
C22 = 100 - C11             % Atomic percent of component 2

% *** OUTPUT ***
C1_con = C11 * atomic_weight_1 / (C11 * atomic_weight_2 + C22 * atomic_weight_2) * 100     % Atomic % to weight %
C11_con = C1 * atomic_weight_2 / (C1 * atomic_weight_2 + C2 * atomic_weight_1) * 100       % Weight & to atomic %
fprintf('The composition in weight percent of component 1 is %.2f %%. \n', C1)
fprintf('The composition in atomic percent of component 1 is %.2f %%. \n', C1_con)


%% Functions 
function weight = calculateWeight(Unit)
    L = strlength(Unit);
    weight = 0;
    mass = 0;
    arr = [];
    num = 0;
    
    for i = 1 : L
        if isletter(Unit(i))
            LL = length(arr);
            if LL ~= 0
                count = 0;
                num = 0;
                for j = linspace(LL, 1, LL)
                    num = num + power(10, count) * arr(j);
                    count = count + 1;
                end
            end
            weight = weight + num * mass;
            arr = [];
            mass = molecularMass(Unit(i));
        else
            arr = [arr round(str2num(Unit(i)))];
        end
    end
    
    count = 0;
    num = 0;
    for j = linspace(LL, 1, LL)
        num = num + power(10, count) * arr(j);
        count = count + 1;
    end
    weight = weight + num * mass;
end

function mass = molecularMass(chr)
    if strcmp(chr, 'C')
        mass = 12.01;
    elseif strcmp(chr, 'H')
        mass = 1.008;
    elseif strcmp(chr, 'N')
        mass = 14.01;
    elseif strcmp(chr, 'O')
        mass = 16.00;
    else
        disp('WARNING: Character specified does not have a known molecular weight')
    end
end

% Get the length the unit cell
function l = getLength(BCC, FCC, simple, r)
    l = 1;
    if BCC
        l = 4 * r / sqrt(3);
    elseif FCC
        l = 2 * sqrt(2) * r;
    elseif simple
        l = 2 * r;
    else
        disp("WARNING: Didn't specify crystal structure")
    end
end

% Get the volume of the unit cell for ceramic 
function V = getCeramicV(BCC, FCC, simple, anion_radius, cation_radius)
    V = 1;
    anion_radius = anion_radius * 1e-7;     % radius in cm
    cation_radius = cation_radius * 1e-7;   % radius in cm
    if BCC
        V = power(2 * (anion_radius + cation_radius) / sqrt(3), 3);
    elseif FCC
        V = power(2 * (anion_radius + cation_radius) / sqrt(2), 3);
    elseif simple
        V = power((anion_radius + cation_radius), 3);
    else
        disp("WARNING: Didn't specify crystal structure")
    end
end

% Get the atoms/unit cell
function n = getN(BCC, FCC, simple)
    n = 0;
    if simple 
        n = 1;
    elseif BCC 
        n = 2;
    elseif FCC
        n = 4;
    else
        disp("WARNING: Didn't specify crystal structure")
    end        
end

function res = checkOne(xi, wi)
    if sum(xi) ~= 1
        disp('WARNING: The sum of xi is not 1')
    end
    if sum(wi) ~= 1
        disp('WARNING: The sum of wi is not 1')
    end
    res = 0;
end

