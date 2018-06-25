%TCP Nash solver for 2 person game

global m A11 A22 numstrat1 numstrat2;

numstrat1 = 2;
numstrat2 = 2;

m = numstrat1 + numstrat2;
% 
%    A11 = tensor(rand(numstrat1,numstrat2,numstrat3));
%    A22 = tensor(rand(numstrat1,numstrat2,numstrat3));
%    A33 = tensor(rand(numstrat1,numstrat2,numstrat3));

A11 = tensor([2 -1;
              -1 1]);
A22 = tensor([1 -1;
              -1 2]);
           
A22 = permute(A22, [2 1]);

          
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

for a = 1:20

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

global m A11 A22 numstrat1 numstrat2;

    strat1 = y_star(1:numstrat1);
    strat2 = y_star(numstrat1 + 1:numstrat1 + numstrat2);

    y1 = ttv(A11, strat2, 2);

    y2 = ttv(A22, strat1, 2);
   
    y = [double(y1);
        double(y2);];
    y = y - 1
end

function y_phi = F_phi(y_star)

    global m A11 A22 numstrat1 numstrat2;

    strat1 = y_star(1:numstrat1);
    strat2 = y_star(numstrat1 + 1:numstrat1 + numstrat2);
  
    y_phi = []; 

    for i = 1:m

        y_star_phi = zeros(m,1);
        y_star_phi(i) = 1;
        
       y1 = ttv(A11, y_star_phi(numstrat1 + 1:numstrat1+numstrat2), 2);
       y2 = ttv(A22, y_star_phi(1:numstrat1), 2);
        
        y = [double(y1);
            double(y2)];
        
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