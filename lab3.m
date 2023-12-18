clear
clc;

% 1. Duomenų paruošimas
X = 0.1:1/22:1;
N = size(X, 2)

f = @(x) (1 + 0.6*sin(2*pi*x/0.7) + 0.3*sin(2*pi*x))/2;
desired_output = f(X);

figure(1)
plot(X, desired_output, 'kx')

% 2. Tinklo struktūros pasirinkimas
% 1 įėjimas; 2 neuronai paslėptajame sluoksnyje, 1 išėjimas
% naudojamos Gauso funkcijos spindulio tipo bazinėse funkcijose

% 3. Pradinių parametrų reikšmių pasirinkimas

% I sluoksnis

c1 = 0.2; r1 = 0.2;
c2 = 0.9; r2 = 0.2;




% II sluoksnis
w1 = randn(1);
w2 = randn(1);
w0 = randn(1);

n = 0.01; % mokymo žingsnis
epoch = 10000; % Mokymo epochų skaičius

% Gauso funkcija
gaussian = @(x, c, r) exp(-(x - c).^2 / (2 * r^2));

% 4. Tinklo atsako skaičiavimas
for k = 1:epoch
    for indx = 1:N
        % I sluoksnis
        y1_1 = gaussian(X(indx), c1, r1);
        y2_1 = gaussian(X(indx), c2, r2);
        
        % II sluoksnis
        y1 = w1 * y1_1 + w2 * y2_1 + w0;

        % 5. Klaidos skaičiavimas E = 1/2*e^2;
        e = desired_output(indx) - y1;

        % 6. Tinklo koeficientų atnaujinimas
        % w = w + n*delta*IN;
        delta_w1 = n * e * y1_1;
        delta_w2 = n * e * y2_1;
        delta_w0 = n * e;

        w1 = w1 + delta_w1;
        w2 = w2 + delta_w2;
        w0 = w0 + delta_w0;
    end
end

% 7. Tinklo testavimas

X_test = 0.1:1/30:1;
N_test = size(X_test, 2);
Y1 = zeros(1, N_test);

desired_output_test = f(X_test);

e = 0;

for indx = 1:N_test
    % I sluoksnis
    y1_1 = gaussian(X_test(indx), c1, r1);
    y2_1 = gaussian(X_test(indx), c2, r2);
    
    % II sluoksnis
    Y1(indx) = w1 * y1_1 + w2 * y2_1 + w0;

    % 8. Klaidos skaičiavimas E = 1/2*e^2;
    e = e + abs(desired_output_test(indx) - Y1(indx));
end

e = e / N_test;

disp("Error " + e)
w0
w1
w2

figure(2)
plot(1:N_test, Y1, 'r', 1:N_test, Y1, 'ko')
hold on 
plot(1:N_test, desired_output_test, 'b', 1:N_test, desired_output_test, 'kx')
legend('predicted', 'predicted points', 'real', 'real points')