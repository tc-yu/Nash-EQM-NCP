%TCP Nash solver for 3 person game
global m A11 A22 A33 numstrat1 numstrat2 numstrat3;

numstrat1 = 2;
numstrat2 = 3;
numstrat3 = 2;

m = numstrat1 + numstrat2 + numstrat3;
% 
%    A11 = tensor(rand(numstrat1,numstrat2,numstrat3));
%    A22 = tensor(rand(numstrat1,numstrat2,numstrat3));
%    A33 = tensor(rand(numstrat1,numstrat2,numstrat3));

A11 = zeros(2,3,2);
A22 = zeros(2,3,2);
A33 = zeros(2,3,2);

A11(:, :, 1) = [0.0605 0.5269 0.6569;
                0.3993 0.4168 0.6280];
A11(:, :, 2) = [0.2920 0.0155 0.1672;
                0.4317 0.9841 0.1062];
A22(:, :, 1) = [0.3724 0.4897 0.9516;
                0.1981 0.3395 0.9203];
A22(:, :, 2) = [0.0527 0.2691 0.5479;
                0.7379 0.4228 0.9427];
A33(:, :, 1) = [0.4177 0.3015 0.6663;
                0.9831 0.7011 0.5391];
A33(:, :, 2) = [0.6981 0.1781 0.9991;
                0.6665 0.1280 0.1711];

A11 = tensor(A11);
A22 = tensor(A22);
A33 = tensor(A33);
            
A22 = permute(A22, [2 1 3]);
A33 = permute(A33, [3 1 2]);
          
mu = 0.1;
epsilon = 0.75;
delta = 10 ^ -4;

y = (0.01) * ones(m,1);
%Fy = F(y);
s = F(y);

size(s)

z = [mu;
    y;
    s]

%phi_val = PHI(z)

h = H(z);

mag = sqrt(dot(h, h))
Beta = (mag / mu);

for a = 1:30

    if mag <= 10 ^ -6 || z(1) == 0
        break
    end
    
    H_phi = H_pi(h, Beta, z);
    smth = [(1 / Beta) * mag;
        zeros(2 * m, 1)]
    delta_z = inv(H_phi) * (smth - h)

    lamda_zk = 1
    power = 0

    while sqrt(dot(H(z + lamda_zk * delta_z), H(z + lamda_zk * delta_z))) > (1 - delta * (1 - (1 / Beta)) * lamda_zk) * mag
        power = power + 1
        lamda_zk = epsilon ^ power
        
        if power > 300
            break;
        end
        
    end

    power = power - 1

    z = z + lamda_zk * delta_z;

    
    h = H(z);

    mag = sqrt(dot(h, h))
end

z

y = transpose(z(2:m+1))
s = transpose(z(m + 2:2 * m+1))

function h = H(z)
global m;

h = [z(1);
    z(m + 2:2 * m + 1) - F(z(2:m+1));
    PHI(z) + z(1) * z(2:m+1)]
end

function y = F(y_star)

global m A11 A22 A33 numstrat1 numstrat2 numstrat3;

    strat1 = y_star(1:numstrat1);
    strat2 = y_star(numstrat1 + 1:numstrat1 + numstrat2);
    strat3 = y_star(numstrat1 + numstrat2 + 1:numstrat1 + numstrat2 + numstrat3);

    y1 = ttv(A11, strat3, 3);
    y1 = ttv(y1, strat2, 2);
    y2 = ttv(A22, strat3, 3);
    y2 = ttv(y2, strat1, 2);
    y3 = ttv(A33, strat2, 3);
    y3 = ttv(y3, strat1, 2);
   
    y = [double(y1);
        double(y2);
        double(y3);]
    y = y - 1
end

function y_phi = F_phi(y_star)

    global m A11 A22 A33 numstrat1 numstrat2 numstrat3;

    strat1 = y_star(1:numstrat1);
    strat2 = y_star(numstrat1 + 1:numstrat1 + numstrat2);
    strat3 = y_star(numstrat1 + numstrat2 + 1:numstrat1 + numstrat2 + numstrat3);
    
    y_phi = []; 

    for i = 1:numstrat1

%    y1 = ttv(A11, y_star_phi(numstrat1 + 1:numstrat1+numstrat2), 2);
%    y2 = ttv(A22, y_star_phi(1:numstrat1), 2);
        
        y_star_phi = zeros(m,1);
        y_star_phi(i) = 1;
       
        y1 = zeros(numstrat1,1);
        y2 = ttv(A22, strat3, 3);
        y2 = ttv(y2, y_star_phi(1:numstrat1), 2);
        y3 = ttv(A33, strat2, 3);
        y3 = ttv(y3, y_star_phi(1:numstrat1), 2);
        
        y = [double(y1);
            double(y2);
            double(y3)]
        y_phi = [y_phi y]
    end
    
    for i = numstrat1 + 1:numstrat1 + numstrat2 
        
        y_star_phi = zeros(m,1);
        y_star_phi(i) = 1;
       
        y1 = ttv(A11, strat3, 3);
        y1 = ttv(y1, y_star_phi(numstrat1 + 1:numstrat1 + numstrat2), 2);
        y2 = zeros(numstrat2,1);
        y3 = ttv(A33, strat1, 2);
        y3 = ttv(y3, y_star_phi(numstrat1 + 1:numstrat1 + numstrat2), 2);
        
        y = [double(y1);
            double(y2);
            double(y3)]
        y_phi = [y_phi y]
    end
    
    for i = numstrat1 + numstrat2 + 1:numstrat1 + numstrat2 + numstrat3
        
        y_star_phi = zeros(m,1);
        y_star_phi(i) = 1;
       
        y1 = ttv(A11, strat2, 2);
        y1 = ttv(y1, y_star_phi(numstrat1 + numstrat2 + 1:numstrat1 + numstrat2 + numstrat3), 2);
        y2 = ttv(A22, strat1, 2);
        y2 = ttv(y2, y_star_phi(numstrat1 + numstrat2 + 1:numstrat1 + numstrat2 + numstrat3), 2);
        y3 = zeros(numstrat3,1);
        
        y = [double(y1);
            double(y2);
            double(y3)]
        y_phi = [y_phi y]
    end
end

function PHI_val = PHI(z)

    global m;
    PHI_val = zeros(m,1)

    for i = 1:m
        PHI_val(i) = phi(z(1), z(i + 1), z(i + 1 + m))
    end
end

function phi_val = phi(mu, y_i, s_i)
    phi_val = y_i + s_i - sqrt((y_i - s_i)^2 + 4 * mu)
end

function h_pi = H_pi(H, Beta, z)

    global m;
    
    mu = z(1);
    y = z(2:m+1);
    s = z(m+2:2 * m +1);
    
    d = zeros(m, 1);
    D = zeros(m);
    E = zeros(m);
    
    for i = 1:m
        d(i) = -2 / (sqrt((y(i) - s(i)) ^ 2 + 4 * mu));
    end
    
    for i = 1:m
        D(i,i) = 1 + mu - (y(i) - s(i)) / (sqrt((y(i) - s(i)) ^ 2 + 4 * mu));
    end
    
    for i = 1:m
        E(i,i) = 1 + (y(i) - s(i)) / (sqrt((y(i) - s(i)) ^ 2 + 4 * mu));
    end
       
    h_pi = [1 zeros(1, m) zeros(1, m);
        zeros(m, 1) -F_phi(y) eye(m);
        (d + y) D E]
end