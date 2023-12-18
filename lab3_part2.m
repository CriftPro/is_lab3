clear
clc;

% 1. Duomenų paruošimas
X = 0.1:1/22:1;
N = size(X, 2);

f = @(x) (1 + 0.6*sin(2*pi*x/0.7) + 0.3*sin(2*pi*x))/2;
desired_output = f(X);

figure(1)
plot(X, desired_output, 'kx')

% 2. Tinklo struktūros pasirinkimas
% 1 įėjimas; 2 neuronai paslėptajame sluoksnyje, 1 išėjimas
% naudojamos Gauso funkcijos spindulio tipo bazinėse funkcijose

% 3. Pradinių parametrų reikšmių pasirinkimas

% I sluoksnis
c1 = randn(1), r1 = randn(1)
c2 = randn(1), r2 = randn(1)

new_centers = [];
[idx new_centers] = kmeans(X',2);

c1 = new_centers(1);
c2 = new_centers(2);

r1 = 0.1:1

% II sluoksnis
w1 = randn(1);
w2 = randn(1);
w0 = randn(1);

n = 0.01; % mokymo žingsnis

alpha = 0.1 ;  % Parametras, nusakančias naujų centrų priėmimo greitį
center_num = 10; % parametru variantu skaicius
epoch = 1000; % Mokymo epochų skaičius

% Gauso funkcija
gaussian = @(x, c, r) exp(-(x - c).^2 / (2 * r^2));

% Parametru fixuotas pasirinkimas ( pradedam nuo 0 iki 1)

fixed = 0 % fixed 1 arba 0 kad butu randomas

if fixed == 1
    
    
    c1_array = [0];
    c2_array = [0];
    r1_array = [0];
    r2_array = [0];

    for i = 1:center_num
        var = c1_array(i) + alpha;
        c1_array = [ c1_array  var];

        var = c2_array(i) + alpha;
        c2_array = [ c2_array  var];

        var = r1_array(i) + alpha;
        r1_array = [ r1_array  var];

        var = r2_array(i) + alpha;
        r2_array = [ r2_array  var];
    end
end

if fixed == 0

    c1_array = zeros(center_num,1)';
    c2_array = zeros(center_num,1)';
    r1_array = zeros(center_num,1)';
    r2_array = zeros(center_num,1)';

end

y_temp = [];
best_error = 100;

r_epoch = 6;
center_num = 26;
alpha = 0.1;


[c1_array,c2_array,r1_array,r2_array] = generating_centers(center_num,alpha,c1,c2,r1,r2,c1_array,c2_array,r1_array,r2_array);

params = choosing_centers([c1 c2 r1  r2],best_error,r_epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0);
c1 = params(1);
c2 = params(2);
r1 = params(3);
r2 = params(4);
best_error = params(5);

params

r_epoch = 3;
center_num = 6;
apha = 0.1;
iter = 0;
previous_error = 0;

while best_error > 0.0001
    
if previous_error == best_error
    apha = alpha + 0.125;
    r_epoch = 6;
end

if iter > 12
    break;
end


[c1_array,c2_array,r1_array,r2_array] = generating_centers(center_num,alpha,c1,c2,r1,r2,c1_array,c2_array,r1_array,r2_array);

params = choosing_centers([c1 c2 r1  r2],best_error,r_epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0);
params;
c1 = params(1);
c2 = params(2);
r1 = params(3);
r2 = params(4);
best_error = params(5);

params = [];

iter = iter + 1;
previous_error = best_error;
best_error;

end


% c1_array;
% c2_array;
% r1_array;
% r2_array;

% total_error = 0;
% best_error = 100;
% 
% y_temp = [];
% params = [];
% 
% r_epoch = 3;
% 
% 
% params = choosing_centers(best_error,r_epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0);
% 
% 
% c1 = params(1);
% c2 = params(2);
% r1 = params(3);
% r2 = params(4);
% 
% alpha = 0.1;
% center_num = 10;
% r_epoch = 5;
% 
% [c1_array,c2_array,r1_array,r2_array] = generating_centers(center_num,alpha,c1,c2,r1,r2,c1_array,c2_array,r1_array,r2_array);
% 
% params = choosing_centers(best_error,r_epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0);
% 
% c1 = params(1);
% c2 = params(2);
% r1 = params(3);
% r2 = params(4);
% 
% alpha = 0.1;
% center_num = 10;
% r_epoch = 10;
% 
% [c1_array,c2_array,r1_array,r2_array] = generating_centers(center_num,alpha,c1,c2,r1,r2,c1_array,c2_array,r1_array,r2_array);
% 
% params = choosing_centers(best_error,r_epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0);
% 
% c1 = params(1);
% c2 = params(2);
% r1 = params(3);
% r2 = params(4);
% 
% alpha = 0.1;
% center_num = 10;
% r_epoch = 15;
% 
% [c1_array,c2_array,r1_array,r2_array] = generating_centers(center_num,alpha,c1,c2,r1,r2,c1_array,c2_array,r1_array,r2_array);
% 
% params = choosing_centers(best_error,r_epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0);
% 
% c1 = params(1);
% c2 = params(2);
% r1 = params(3);
% r2 = params(4);







% for i_r1 = 2:size(r1_array,2)
%     for i_r2 = 2:size(r2_array,2)
%         for i_c2 = 2:size(c2_array,2)
%             for i_c1 = 2:size(c1_array,2)
%                 [y_temp,w1,w2,w0] = rbf(y_temp,gaussian,epoch,N,X,desired_output,n,w1,w2,w0,c1_array(i_c1),c2_array(i_c2),r1_array(i_r1),r2_array(i_r2));
%         
%                 total_error = (min(desired_output - y_temp))^2;
%                 total_error;
% 
%                 if total_error < best_error
%                     best_error = total_error;
%                     i_r1;
%                     i_r2;
%                     i_c2;
%                     best_params = [ c1_array(i_c1) c2_array(i_c2) r1_array(i_r1) r2_array(i_r2)];
% 
%                     
%                 end
%             end
%         end
%     end
% end




disp("Error " + best_error)
disp("Generated coefficents : ")
disp("r1 = " + r1);
disp("r2 = " + r2);
disp("c1 = " + c1);
disp("c2 = " + c2);


[y_temp,w1,w2,w0] = rbf(y_temp,gaussian,epoch,N,X,desired_output,n,w1,w2,w0,c1,c2,r1,r2);



% w0 =
% 
%    -0.3467
% 
% 
% w1 =
% 
%     1.2767
% 
% 
% w2 =
% 
%     1.0562

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



figure(2)
plot(1:N_test, Y1, 'r', 1:N_test, Y1, 'ko')
hold on 
plot(1:N_test, desired_output_test, 'b', 1:N_test, desired_output_test, 'kx')
legend('predicted', 'predicted points', 'real', 'real points')

function [c1_array,c2_array,r1_array,r2_array] = generating_centers(center_num,alpha,c1,c2,r1,r2,c1_array,c2_array,r1_array,r2_array)

    center_pos = center_num/2;
    

    c1_array(center_pos) = [c1];
    c2_array(center_pos) = [c2];
    r1_array(center_pos) = [r1];
    r2_array(center_pos) = [r2];

    for i = center_pos:center_num-1
        
        var = c1_array(i) + alpha;
        c1_array(i+1) = var;

        var = c2_array(i) + alpha;
        c2_array(i+1) = var;

        var = r1_array(i) + alpha;
        r1_array(i+1) = var;

        var = r2_array(i) + alpha;
        r2_array(i+1) = var;
    end

    for i = center_pos:-1:2
        
        var = c1_array(i) - alpha;
        c1_array(i-1) = var;

        var = c2_array(i) - alpha;
        c2_array(i-1) = var;

        var = r1_array(i) - alpha;
        r1_array(i-1) = var;

        var = r2_array(i) - alpha;
        r2_array(i-1) = var;
    end
end



function [y_temp, w1,w2,w0,total_err] = rbf(y_temp,gaussian,epoch,N,X,desired_output,n,w1,w2,w0,c1,c2,r1,r2)

total_err = 0;

for k = 1:epoch
    for indx = 1:N
        % I sluoksnis
        y1_1 = gaussian(X(indx), c1, r1);
        y2_1 = gaussian(X(indx), c2, r2);

        % II sluoksnis
        y1 = w1 * y1_1 + w2 * y2_1 + w0;
        y_temp(indx) = y1;
        % 5. Klaidos skaičiavimas E = 1/2*e^2;
        e = desired_output(indx) - y1;
        total_err = total_err + abs(e);
        % 6. Tinklo koeficientų atnaujinimas
        % w = w + n*delta*IN;
        delta_w1 = n * e * y1_1;
        delta_w2 = n * e * y2_1;
        delta_w0 = n * e;

        w1 = w1 + delta_w1;
        w2 = w2 + delta_w2;
        w0 = w0 + delta_w0;

        %             r1_error = r1_error + e^2;
    end
end
end


function params = choosing_centers(best_params,best_error,epoch,c1_array,c2_array,r1_array,r2_array,y_temp,gaussian,n,N,X,desired_output,w1,w2,w0)
    
    total_error = 0;
    
    

    for i_r1 = 2:size(r1_array,2)
        for i_r2 = 2:size(r2_array,2)
            for i_c2 = 2:size(c2_array,2)
                for i_c1 = 2:size(c1_array,2)
                    [y_temp,w1,w2,w0,total_error] = rbf(y_temp,gaussian,epoch,N,X,desired_output,n,w1,w2,w0,c1_array(i_c1),c2_array(i_c2),r1_array(i_r1),r2_array(i_r2));

%                     total_error = (min(desired_output - y_temp))^2;
                     total_error;

                    if total_error < best_error
                        best_error = total_error;
                        i_r1;
                        i_r2;
                        i_c2;
                        best_params = [ c1_array(i_c1) c2_array(i_c2) r1_array(i_r1) r2_array(i_r2)];

                        
                    end
                end
            end
        end
    end
    
    params = [ best_params best_error];

end
