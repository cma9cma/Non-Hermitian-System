close all
clear
clc

% parameters
rho = 1;
gamma = 0.1;
kappa = 1.0;
num_points = 100000;
t = linspace(0, 2 * pi / abs(gamma), num_points);
dt = t(2) - t(1);
g = 1 - rho * cos(gamma * t);
delta = rho * sin(gamma * t);
lda1 = zeros(1, num_points);
lda2 = zeros(1, num_points);
c1 = zeros(1, num_points);
c2 = zeros(1, num_points);

% time evolution encircling EP at g = 1, delta = 0
for i = 1:length(t)
    Ham = [-1j * g(i) - delta(i), kappa;
        kappa, 1j * g(i) + delta(i)];
    [evector, evalue] = eig(Ham);

    if i == 1

        if evalue(1, 1) < 0
            lda1(i) = evalue(1, 1);
            lda2(i) = evalue(2, 2);
            evector_base1 = evector(:, 1);
            evector_base2 = evector(:, 2);
        else
            lda1(i) = evalue(2, 2);
            lda2(i) = evalue(1, 1);
            evector_base1 = evector(:, 2);
            evector_base2 = evector(:, 1);
        end

        % choose the initial state as state 1
        current_state = evector_base1;
    else

        if abs(evalue(1, 1) - lda1(i - 1)) <= abs(evalue(2, 2) - lda1(i - 1))
            lda1(i) = evalue(1, 1);
            lda2(i) = evalue(2, 2);
        else
            lda1(i) = evalue(2, 2);
            lda2(i) = evalue(1, 1);
        end

        % the new state is governed by the schrodinger equation
        current_state = current_state + (-1j) * Ham * current_state * dt;
    end

    c1(i) = dot(evector_base1, current_state);
    c2(i) = dot(evector_base2, current_state);

end

% visualization
res = 600;
plot_weight = 1.3;
fig_eigenval = figure;
plot(real(lda1), imag(lda1), 'red', 'linewidth', plot_weight);
hold on
plot(real(lda2), imag(lda2), 'blue', 'linewidth', plot_weight);
hold off
legend('\lambda_1', '\lambda_2');
set(get(gca, 'XLabel'), 'String', 'real(\lambda)');
set(get(gca, 'YLabel'), 'String', 'imag(\lambda)');
set(get(gca, 'Title'), 'String', ['evolution of eigenvalue, ', 'gamma = ', num2str(gamma)]);
filename = strcat(['fig_eigenval_', 'gamma_', num2str(gamma)]);
filename = [strrep(filename, '.png', '') '.png'];
print(fig_eigenval, '-dpng', ['-r' num2str(res)], filename);

fig_state = figure;
plot(t, log10(abs(c1).^2), 'red', 'linewidth', plot_weight);
hold on
plot(t, log10(abs(c2).^2), 'blue', 'linewidth', plot_weight);
hold off
legend('|c_{1}|^2', '|c_{2}|^2', 'location', 'southeast');
set(get(gca, 'XLabel'), 'String', 't');
set(get(gca, 'YLabel'), 'String', 'log of amplitude');
set(get(gca, 'Title'), 'String', ['evolution of state, ', 'gamma = ', num2str(gamma)]);
filename = strcat(['fig_state_', 'gamma_', num2str(gamma)]);
filename = [strrep(filename, '.png', '') '.png'];
print(fig_state, '-dpng', ['-r' num2str(res)], filename);
