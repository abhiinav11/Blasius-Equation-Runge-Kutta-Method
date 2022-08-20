L = 10.0; %length of the domain
N = 100000; %number of grid points
A = 0.0; %guess value of H at eta = 0
B = 1.0; %guess value of H at eta = 0
convergence= 0.00001; %convergence criteria
h = L/(N-1); %step size
f_array = zeros(1,N);
g_array = zeros(1,N);
h_array = zeros(1,N);
Theta = zeros(1,N);
ThetaDes = zeros(1,N);


iterations = 0;
c_old=1.0;

a=A; %These are guessed values of H provided above
b=B;

%bisection method
error = 1;
eta = linspace(0,L,N);

figure(1)
title("Plot for u/U_{\infty} vs \eta")
xlabel("\eta \rightarrow")
ylabel("u/U_{\infty} \rightarrow")
hold on

figure(2)
title("Plot for \theta vs \eta")
xlabel("\eta \rightarrow")
ylabel("\theta \rightarrow")
hold on

while(error > convergence)

    iterations = iterations + 1; %counter for number of iterations
    c_new= (a+b) / 2.0; %new guess value
    [f_a, f_array, g_array, h_array]=shootingVelocity(a,f_array,g_array,h_array,h,N); %solve using shooting velocity technique
    [f_b, f_array, g_array, h_array]=shootingVelocity(b,f_array,g_array,h_array,h,N);
    [f_c, f_array, g_array, h_array]=shootingVelocity(c_new,f_array,g_array,h_array,h,N);

    if(((f_a-1.0)*(f_c-1.0))>0) %check where root lies
        a=c_new;
    else
        b=c_new;
    end 

    error=abs((c_new-c_old)/c_old); %calculate error
    c_old=c_new; %update guess value
end 
figure(1)
plot(eta,g_array)



for Pr = [0.7 0.8 1 7]
    a=A; %Guess value of ThetaDes(0) provided
    b=B;
    c_old = 1.0;
    error = 1;
    Theta = zeros(1,N);
    ThetaDes = zeros(1,N);
    while(error>convergence)
        iterations = iterations + 1; %counter for number of iterations
        c_new=(a+b)/2.0; %new guess value
        [f_a, Theta, ThetaDes]= shootingThermal(f_array,a,Theta,ThetaDes,Pr,h,N); %%% solve using shoootingThermal technique
        [f_c, Theta, ThetaDes]= shootingThermal(f_array,c_new ,Theta,ThetaDes,Pr,h,N);
        [f_b, Theta, ThetaDes]= shootingThermal(f_array,b,Theta,ThetaDes,Pr,h,N);
        if((f_a-1.0)*(f_c-1.0)>0) %check where root lies
            a=c_new;
        else
            b=c_new;
        end
        error=abs((c_new-c_old)/c_old); %calculate error
        c_old = c_new; %update old guess
    end

    figure(2)
    plot(eta,1-Theta)   %%%%for theta = (T-T_inf)/(T_w - T_inf) = 1-(T_w - T)/(T_w - T_inf)
    hold on
end



%%%%%Function Definition

%shoots using guess h(0) value, returns g(L) value
function [g_new, f_array, g_array, h_array] = shootingVelocity(h_old,f_array,g_array,h_array,h,N)
    f_old=0; %f(0) = 0
    g_old=0; %g(0) = 0
    for i = 1:N %calculate coefficients of Runge-Kutta method
        k1=h*g_old;
        l1=h*h_old;
        m1=h*-0.5*f_old*h_old;

        k2=h*(g_old+0.5*l1);
        l2=h*(h_old+0.5*m1);
        m2=h*-0.5*(f_old+0.5*k1)*(h_old+0.5*m1);

        k3=h*(g_old+0.5*l2);
        l3=h*(h_old+0.5*m2);
        m3=h*-0.5*(f_old+0.5*k2)*(h_old+0.5*m2);

        k4=h*(g_old+l3);
        l4=h*(h_old+m3);
        m4=h*-0.5*(f_old+k3)*(h_old+m3);

        %values of f,g,h for next position
        f_new=f_old+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
        g_new=g_old+(1.0/6.0)*(l1+2.0*l2+2.0*l3+l4);
        h_new=h_old+(1.0/6.0)*(m1+2.0*m2+2.0*m3+m4);
        
        %update f,g,h
        f_old=f_new;
        g_old=g_new;
        h_old=h_new;
        
        %update f,g,h in array
        f_array(i) = f_new;
        g_array(i) = g_new;
        h_array(i) = h_new;
    end
end 


function [Y1_new, Theta, ThetaDes] = shootingThermal(f_array,Y_in,Theta,ThetaDes,Pr,h,N)
    Y_old=Y_in; %Guess value for ThetaDes(0)
    Y1_old=0.0; %Theta(0) = 0

    for i = 1:N  %calculate coefficients of Runge-Kutta method
        k1=-0.5*Pr*f_array(i)*Y_old;
        k2=-0.5*Pr*f_array(i)*(Y_old+(h/2.0)*k1);
        k3=-0.5*Pr*f_array(i)*(Y_old+(h/2.0)*k2);
        k4=-0.5*Pr*f_array(i)*(Y_old+h*k3);

        l1=Y_old;
        l2=Y_old+(h/2.0)*l1;
        l3=Y_old+(h/2.0)*l2;
        l4=Y_old+h*l3;
        
        %values of Theta, ThetaDes for next position
        Y_new=Y_old+(h/6.0)*(k1+2.0*k2+2.0*k3+k4);
        Y1_new=Y1_old+(h/6.0)*(l1+2.0*l2+2.0*l3+l4);
        
        %update Theta, ThetaDes
        Y_old=Y_new;
        Y1_old=Y1_new;
        
        %update Theta, ThetaDes in array
        Theta(i) = Y1_new;
        ThetaDes(i) = Y_new;
    end
end

