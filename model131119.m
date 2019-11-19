
close all

%% Data

dt = 0.01;          % Step time
T = 3;              % Time of the experiment
N = T / dt;         % Number of steps
kv = 0.1;           
lamda = 0;
tau = 50*10^(-3);   %constante de temps
dim  = 14;        %dimension des matrices
nb_of_simu = 10;
%x = [x x' y y' Fx Fy x* x'* y* y'* Fextx Fexty] 
xi = [0 0 0 0 0 0 0 0 0.3 0 0 0];
xdev=0.1;
ydev=0.1;
xobs = [xdev ydev];


w1=10;
w2=10;        
w3=1;
w4=1;
w5=1000;
w6=1000;
 

% Initialization

A = eye(dim);
B = zeros(dim,2);

A(1,2) = dt;A(2,2) = 1-dt*kv;A(2,5) = dt;A(3,4) = dt;
A(4,4) = 1-dt*kv;A(4,6) = dt;A(5,5) = 1-dt/tau;
A(6,6) = 1-dt/tau; 

B(5,1) = dt/tau;B(5,2) = (lamda*dt)/tau;B(6,1) = (lamda*dt)/tau;B(6,2) = dt/tau;

Q = zeros(dim);
Q(1,1) = w1; Q(2,2) = w3; Q(3,3) = w2; Q(4,4) = w4;
Q(7,7) = w1; Q(8,8) = w3; Q(9,9) = w2; Q(10,10) = w4;
Q(1,7) = -w1; Q(2,8) = -w3; Q(3,9) = -w2; Q(4,10) = -w4;
Q(7,1) = -w1; Q(8,2) = -w3; Q(9,3) = -w2; Q(10,4) = -w4;
Qdev = zeros(dim);
Qdev(1,1)  = w5;Qdev(3,3)  = w6;Qdev(1,13) = -w5;Qdev(3,14) = -w6;
Qdev(13,1) = -w5;Qdev(14,3) = -w6;Qdev(13,13) = w5;Qdev(14,14) = w6;

R = eye(2)*10^(-5);    
H = eye(dim);        

S     = zeros(dim,dim,N);
L     = zeros(2,dim,N);
Sigma = zeros(dim,dim,N);
K     = zeros(dim,dim,N);                     
u     = zeros(2,N,nb_of_simu);   % Control variable matrix
x     = zeros(dim,N,nb_of_simu);
xstar   = zeros(dim,N,nb_of_simu); % x with feedback control(FB)



% Noise proportionnal to B according to the document
Omegaxsi = 0.1*(B*B');                     
Omegaw = 0.1*max(max(Omegaxsi))*eye(dim);  

% Initial condition
S(:,:,N) = Q;
s = 0;
Sigma(:,:,1) = Omegaw;



% Model computation

for i=N-1:-1:1
    L(:,:,i) = (R+B.'*S(:,:,i+1)*B)\(B.'*S(:,:,i+1)*A);
    S(:,:,i) =  A'*S(:,:,i+1)*(A-B*L(:,:,i));
    if i = 150
      S(:,:,i) = S(:,:,i) +Qdev;
      end
    end
    s = s + trace(S(:,:,i+1)+Omegaxsi);

for p = 1:nb_of_simu
    % Initialization
    xstar(:,1,p) = reshape([xi xobs],14,1,1);    
    x(:,1,p) = [xi xobs];
    
    for k=1:N-1
        % Generation of the sensory and motor noises% remettre les bruits 
        %motor_noise = mvnrnd(zeros(dim,1),Omegaxsi)';
        %sensory_noise = mvnrnd(zeros(dim,1),Omegaw)'; 
        % Kalman controller
        K(:,:,k) = (A*Sigma(:,:,k)*H')/(H*Sigma(:,:,k)*H' + Omegaw);
        Sigma(:,:,k+1) = Omegaxsi*25 + (A - K(:,:,k)*H)*Sigma(:,:,k)*A';
        u(:,k,p) = -L(:,:,k)*xstar(:,k,p);
        %Calcul de la matrice des variables d'�tats
        y = H*x(:,k,p); %+ sensory_noise;
        %Calcul de l'�tat suivant de x en consid�rant qu'il y avait un contr�le de r�troaction 
        x(:,k+1,p) = A*x(:,k,p) + B*u(:,k,p) ;%+ motor_noise;
        %Calcul de x avec la commande de r�troaction
        xstar(:,k+1,p) = A*xstar(:,k,p) + B*u(:,k,p) + K(:,:,k)*(y - H*xstar(:,k,p));
        
    end
end

plot(x(1,:,1),x(3,:,1))