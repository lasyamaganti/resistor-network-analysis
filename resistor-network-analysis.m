clear;
close all;

resistors = read_resistances('node_resistances.txt');
disp(resistors);

fixed_voltages = read_fixed_voltages('node_voltages.txt');
disp(fixed_voltages);

% figure out number of nodes N from the largest index seen
N = max([resistors(:,1); resistors(:,2); fixed_voltages(:,1)] );

% A for the 5x5 matrix, n is 25
A = compute_A(25, resistors, fixed_voltages);

b = zeros(N,1);
if ~isempty(fixed_voltages)
    b(round(fixed_voltages(:,1))) = fixed_voltages(:,2); % put fixed voltages in b
end

[L,U,P] = lu_fact(A); % LU with partial pivoting
y = forward_sub(L, P*b); % solve Ly = P b (forward sub)
x = back_sub(U, y); % solve Ux = y (back sub)

% from node i to node j
% I{i->j} = (x(i) - x(j)) / R_ij
M = size(resistors,1);
I = zeros(M,1);
for t = 1:M
    i = round(resistors(t,1));
    j = round(resistors(t,2));
    R = resistors(t,3);
    I(t) = (x(i) - x(j)) / R;
end

% open new .txt file for results
fid = fopen('results.txt','w');

fprintf(fid, 'A (size %dx%d):\n', size(A,1), size(A,2));
for i = 1:size(A,1)
    fprintf(fid, '% .6g ', A(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, 'b (size %dx1):\n', length(b));
for i = 1:length(b)
    fprintf(fid, '% .10g\n', b(i));
end
fprintf(fid, '\n');

fprintf(fid, 'LU factorization (P*A = L*U)\n');

fprintf(fid, 'P:\n');
for i = 1:size(P,1)
    fprintf(fid, '% .6g ', P(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, 'L:\n');
for i = 1:size(L,1)
    fprintf(fid, '% .6g ', L(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, 'U:\n');
for i = 1:size(U,1)
    fprintf(fid, '% .6g ', U(i,:));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

fprintf(fid, 'Node voltages (index : value)\n');
for k = 1:length(x)
    fprintf(fid, '%3d : %.10g\n', k, x(k));
end

fprintf(fid, '\nLink currents (i  j  R[Ohm]  I[i->j][A])\n');
for t = 1:size(resistors, 1)
    fprintf(fid, '%3d %3d  %10.6g  % .10g\n', ...
        round(resistors(t,1)), round(resistors(t,2)), resistors(t,3), I(t));
end

fclose(fid);


%% Plotting Step 8
% Part i: Total power vs. uniform link resistance

k = 5; % grid size
N = k*k;
R0 = 1.0; % all links at 1 Ohm

% Build resistor list: [i j R]
links = zeros(2*k*(k-1), 3);   % edges
t = 0;
for r = 1:k
    for c = 1:k
        i = (r-1)*k + c;
        if c < k
            t = t+1; links(t,:) = [i, i+1, R0]; %horizontal edge
        end
        if r < k
            t = t+1; links(t,:) = [i, i+k, R0]; %vertical edge
        end
    end
end

% A and b with fixed voltages. V(1)=1, V(N)=0
A = zeros(N,N); 
b = zeros(N,1);
for e = 1:size(links,1)
    i = links(e,1); 
    j = links(e,2); 
    R = links(e,3); 
    g = 1/R;
    A(i,i) = A(i,i) + g;  
    A(j,j) = A(j,j) + g;
    A(i,j) = A(i,j) - g;  
    A(j,i) = A(j,i) - g;
end
A(1,:) = 0; 
A(1,1) = 1; 
b(1) = 1;
A(N,:) = 0; 
A(N,N) = 1; 
b(N) = 0;

% solve using LU and forward/back sub
[L,U,P] = lu_fact(A);
y = forward_sub(L, P*b);
x = back_sub(U, y);

% compute S = (Vi - Vj)^2
S = 0;
for e = 1:size(links,1)
    i = links(e,1); j = links(e,2);
    S = S + (x(i) - x(j))^2;
end

% Power vs. R
Rvals = linspace(0.2, 10, 200);
Pvals = S ./ Rvals;

figure; plot(Rvals, Pvals, 'LineWidth', 1);
xlabel('Uniform link resistance R (\Omega)');
ylabel('Total power dissipated P (W)');
title(sprintf('Total power vs. R for %d\\times%d grid', k, k));
grid on;

% Part ii: current into bottom-right node vs. k

Ks = [5 10 15 20 25];
I_in = zeros(size(Ks));

for idx = 1:numel(Ks)
    k = Ks(idx);
    N = k*k;
    R0 = 1.0;

    % resistor list based on grid size
    resistors = zeros(2*k*(k-1), 3);
    t = 0;
    for r = 1:k
        for c = 1:k
            i = (r-1)*k + c;
            if c < k
                t = t+1; resistors(t,:) = [i, i+1, R0];
            end
            if r < k
                t = t+1; resistors(t,:) = [i, i+k, R0];
            end
        end
    end

    % A and b with V(1)=1, V(N)=0
    A = zeros(N,N); b = zeros(N,1);
    for e = 1:size(resistors,1)
        i = resistors(e,1); 
        j = resistors(e,2); 
        R = resistors(e,3); 
        g = 1/R;
        A(i,i) = A(i,i) + g;  
        A(j,j) = A(j,j) + g;
        A(i,j) = A(i,j) - g; 
        A(j,i) = A(j,i) - g;
    end
    A(1,:) = 0; 
    A(1,1) = 1; 
    b(1) = 1;
    A(N,:) = 0; 
    A(N,N) = 1; 
    b(N) = 0;

    % solve
    [L,U,P] = lu_fact(A);
    y = forward_sub(L, P*b);
    x = back_sub(U, y);

    % incoming current at bottom-right node
    incoming = 0;
    if mod(N-1, k) ~= 0
        j = N-1; 
        incoming = incoming + (x(j) - x(N))/R0;
    end
    if N-k >= 1
        j = N-k; 
        incoming = incoming + (x(j) - x(N))/R0;
    end
    I_in(idx) = incoming;
end

figure; plot(Ks, I_in, 'LineWidth', 1, 'MarkerSize', 6);
xlabel('Grid size k');
ylabel('Incoming current at bottom-right node (A)');
title('Current vs. grid size (R=1, V_{TopLeft}=1, V_{BottomRight}=0)');
grid on;

%% Step 10

Ks = [5 10 15 20 25]; % grid sizes
Rvals = linspace(0.1, 10, 200); % uniform link resistance values
Reff_vs_R = cell(size(Ks)); % store curves for each k

for idx = 1:numel(Ks)
    k = Ks(idx);
    N = k*k;

    % k×k grid with R=1Ω on every link
    links = zeros(2*k*(k-1), 3);
    t = 0;
    for r = 1:k
        for c = 1:k
            i = (r-1)*k + c;
            if c < k
                t = t+1; links(t,:) = [i, i+1, 1.0];
            end
            if r < k
                t = t+1; links(t,:) = [i, i+k, 1.0];
            end
        end
    end

    % Laplacian for R=1Ω
    L = zeros(N,N);
    for t = 1:size(links,1)
        i = links(t,1); j = links(t,2); g = 1.0;  % R=1 => g=1
        L(i,i) = L(i,i) + g;
        L(j,j) = L(j,j) + g;
        L(i,j) = L(i,j) - g;
        L(j,i) = L(j,i) - g;
    end

    % pseudoinverse and base effective resistance between 1 and N
    Lplus = pinv(L);
    Reff1 = Lplus(1,1) + Lplus(N,N) - 2*Lplus(1,N);  % at R=1Ω

    % scaling
    Reff_vs_R{idx} = Reff1 * Rvals;
end

% plot
figure; hold on;
for idx = 1:numel(Ks)
    plot(Rvals, Reff_vs_R{idx}, 'LineWidth', 1.6, 'DisplayName', sprintf('k=%d', Ks(idx)));
end
xlabel('Link resistance R (\Omega)');
ylabel('Effective resistance R_{eff} (\Omega)');
title('Effective resistance vs. link resistance');
grid on; 
legend('Location','northwest');


%% Functions

% function to read the given txt file with resistances
function resistors = read_resistances(node_resistances)

    fid = fopen(node_resistances, 'r');
    if fid == -1
        error('Could not open file %s', node_resistances);
    end

    data = textscan(fid, '%f %f %f'); % '%f' keeps the resistances as decimals
    fclose(fid);

    resistors = [data{1}, data{2}, data{3}]; % keeps the new variable in the correct structure.
end

% function to read given fixed voltage values
function fixed = read_fixed_voltages(node_voltages)

    fid = fopen(node_voltages, 'r');
    if fid == -1
        error('Could not open file %s', node_voltages);
    end

    data = textscan(fid, '%f %f');  
    fclose(fid);

    fixed = [data{1}, data{2}];
end

% function to compute matrix using the previous read functions.
function A = compute_A(N, resistors, fixed)
% N is the number of nodes
% resistors and fixed are the data the previous two functions read
% the matrix A has dimensions of NxN and is the coefficient for Ax = b

% initialize an NxN zero matrix
A = zeros(N,N);

for t = 1:size(resistors,1) % loops through each row of the resistors array

    i = resistors(t,1); % first column of resistors array with first node number
    j = resistors(t,2); % second column of resistors array with second node number
    R = resistors(t,3); % third column of resistors array with resistor values

    g = 1.0 / R; % compute conductnce of each resistance value

    % updates using Laplacian rule. Fills A matrix with conductance
    A(i,i) = A(i,i) + g;
    A(j,j) = A(j,j) + g;
    A(i,j) = A(i,j) - g;
    A(j,i) = A(j,i) - g;
end

% checks if the fixed_voltages is nonempty and if it is, apply Dirichlet
% rows through row replacement
if ~isempty(fixed)
    fixed_nodes = fixed(:,1);  % nodes to fix
    for k = fixed_nodes.'
        A(k,:) = 0; % zero the row
        A(k,k) = 1; % put 1 on the diagonal
    end
end
end

% Compute the LU factorization
function [L,U,P] = lu_fact(A)
% L is the unit lower-triangle
% U is the upper triangle
% P is a permutation matrix that allows P*A = L*U with partial pivoting

    [m,n] = size(A);
    
    if m ~= n 
        error('A must be square.'); 
    end

    U = A;
    L = eye(n); %nxn identity matrix
    P = eye(n); %nxn identity matrix. row swaps will accumulate in the matrix

    for k = 1:n-1 % iterate through columns
        % looks everything below and including the pivot position (row k-n, column k). Find the row p with the largest magnitude in that column
        [~,p_rel] = max(abs(U(k:n,k))); 
        p = p_rel + k - 1;

        % Swap rows in U and P; and in L's left block (columns 1..k-1)
        if p ~= k % if the best pivot row is not k
            U([k p],k:n) = U([p k],k:n); % swap rows k and p of U, but only from columns k->n
            P([k p],:) = P([p k],:); % swap rows k and p of P, the permutation matrix
            if k > 1
                L([k p],1:k-1) = L([p k],1:k-1); % swap rows k and p of L to keep the corectness of the previous columns swaps
            end
        end

        % % Eliminate below the pivot
        % if U(k,k) == 0
        %     error('Zero pivot encountered at column %d.', k);
        % end

        % Gaussian elimination for column k
        for i = k+1:n
            L(i,k) = U(i,k) / U(k,k); % Compute the multiplier that zeros out anything below the pivot
            U(i,k:n) = U(i,k:n) - L(i,k) * U(k,k:n); % row operation that sets U(i,k) to 0 and updated the rest of the row
        end
    end
end

function y = forward_sub(L, b)
% Solve L*y = b for lower-triangular L.
    n = length(b);
    y = zeros(n,1);
    for i = 1:n
        y(i) = b(i) - L(i,1:i-1) * y(1:i-1); % If L is not unit-lower
    end
end

function x = back_sub(U, y)
% Solve U*x = y for upper-triangular U.
    n = length(y);
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = (y(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
    end
end

% effective resistance matrix
function Reff = effective_resistance_matrix(filename, N)

    fid = fopen(filename,'r');
    data = textscan(fid, '%f %f %f');  % i, j, R
    fclose(fid);
    links = [data{1}, data{2}, data{3}];

    L = zeros(N,N);
    for t = 1:size(links,1)
        i = round(links(t,1));
        j = round(links(t,2));
        R = links(t,3);
        g = 1.0 / R;
        L(i,i) = L(i,i) + g;
        L(j,j) = L(j,j) + g;
        L(i,j) = L(i,j) - g;
        L(j,i) = L(j,i) - g;
    end

    % pseudoinverse of Laplacian
    Lplus = pinv(L);

    % all pairs effective resistance
    d = diag(Lplus);
    Reff = d + d.' - 2*Lplus;
end

%% Understanding the results

% Substituting Ohm’s law, I = (Vi - Vj)/Rij, each node gives one linear equation. This produces a system of equations A*x = b, where:
% - A is the Laplacian (conductance) matrix of the network,
% - x is the vector of unknown node voltages,
% - b encodes the fixed (Dirichlet) boundary conditions.
%
% The resistor file defines the edges and resistance values. The fixed voltage file defines which nodes have which voltages at them.
%
% Each resistor between nodes i and j contributes to A:
% Aii += 1/Rij
% Ajj += 1/Rij
% Aij -= 1/Rij
% Aji -= 1/Rij
% If a node has a fixed voltage, its row in A is replaced with an identity row and the corresponding entry of b is set to the fixed voltage.
%
% To solve A*x = b, LU factorization with partial pivoting is used: P*A = L*U; where L is lower-triangular, U is upper-triangular, and P is a permutationmatrix. 
% The solution happens in two stages:
% L*y = P*b (forward substitution)
% U*x = y (backward substitution)
%
% Once voltages are found, currents through each resistor are computed by:
% Iij = (Vi - Vj)/Rij
%
% The results (matrix A, vector b, factors L,U,P, node voltages, and link currents) are written to a file so they can be checked and interpreted.
%
% From the voltages and currents we can also compute power and scaling:
% Total power in the network is:
% P = sum((Vi - Vj)^2 / Rij) over all edges (i,j).
% For uniform resistances, this scales like 1/R.
% Current into the bottom-right node can be studied as grid size grows (k = 5,10,15,...). Larger grids give more parallel paths, so more incoming current.
%
% Effective resistance between two nodes (top-left and bottom-right) is computed from the pseudoinverse of the Laplacian:
% Reff(a,b) = L+_aa + L+_bb - 2*L+_ab
% For uniform resistors, effective resistance scales linearly with R.
%
% Overall, this approach connects circuit laws to linear algebra: the node
% voltages come from solving A*x = b, currents follow from Ohm’s law, and
% power and effective resistance give physical insights into how the grid
% behaves.