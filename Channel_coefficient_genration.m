N = 20;
lower_bound = 225;
upper_bound = 275;
c = 3 * 10^8;
fc = 1.8 * 10^9;
v = 138.33;

% Generate uniform samples for tn and theta
uniform_samples = lower_bound + (upper_bound - lower_bound) * rand(1, N);
tn = uniform_samples / c;
theta = 2 * pi * rand(1, N);

% Calculate fdn
fm = (v / c) * fc;
fdn = fm * cos(theta);

cn = ones(1, N);

t = 0:0.001:1;
r = length(t);

g = zeros(1, r); % Initialize the sum result

for j = 1:r
    for i = 1:N
        % Calculate the exponential term
        phi_n(j) = exp(-1j * 2 * pi * ((fc + fdn(i)) * tn(i) - fdn(i) * t(j)));
        % Calculate the sum
        g(j) = g(j) + cn(i) * phi_n(j);
    end
end

% Display the sum_result
disp(g);



% Calculate the modulus of g
modulus_g = abs(g);

% Plot both histograms in the same figure
figure;
subplot(1,2,1); % Subplot 1
histogram(modulus_g, 'Normalization', 'probability');
xlabel('Magnitude of g');
ylabel('Probability');
title('Distribution of Magnitude of g');

subplot(1,2,2); % Subplot 2
histogram(real(g), 'Normalization', 'probability');
xlabel('Magnitude of g');
ylabel('Probability');
title('Distribution of Real part of g');

% lags =1;
% 
% autocorr_g=[];
% for m = 1:r
%     % Calculate the autocorrelation at lag i
%     autocorr_g(m) =abs(g(m)) .* abs((g(m-lags)));
% end
% 
% % autocorrelation
% figure;
% plot(lags, autocorr_g);
% xlabel('Lag');
% ylabel('Autocorrelation');
% title('Autocorrelation of g');