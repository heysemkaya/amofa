function r = mvnplot(mu,sigma);
%MVNplot(mu,sigma) Draws the covariance of the multivariate normal distribution in 2D.
%   MU is the mean vector
%   SIGMA is a symmetric positive definite matrix with size equal to the 
%   length of MU
% use plot(r(:,1),r(:,2),'-c') with the result, where c is the colour
%
% 3 February 2003 Albert Ali Salah 

[m1 n1] = size(mu);
c = max([m1 n1]);
if m1 .* n1 ~= c
   error('Mu must be a vector.');
end

[m n] = size(sigma);
if m ~= n
   error('Sigma must be square');
end

if m ~= c
   error('The length of mu must equal the number of rows in sigma.');
end

%cholesky factorization of the covariance matrix
[T p] = chol(sigma);
if p ~= 0
  error('Sigma must be a positive definite matrix.');
end


if m1 == c
  mu = mu';
end

mu = mu(ones(201,1),:);

%unit circle:
xx = 0:pi/100:pi*2;
yx = sin(xx);
yy = sqrt(1 - yx.^2);
yy(52:150) = -yy(52:150);

%translate the unit circle to the desired location, and reshape it
r = [yx' yy'] * T + mu;

