
addpath('C:\Users\wconway\Documents\MATLAB\AM217Project')
%Set Initial Parameters
Number_Time_Steps=1e4; %How many timesteps we let the code run for
Number_Hec1=10; %How many Hec1 diffuse
Prob_Fall=1e-4; %Related to k_unbind, but more annoying to do the math
Prob_Bind=0.4; %Related to k_bind, not quite clear how
Tether_Length=30; %Max distance from (0,0,0) the Hec1 can diffuse
Binding_Distance=2.5; %How close you need to be to the microtubule to try and bind)

%Preallocate 
RandSteps=rand(1,Number_Time_Steps*3*Number_Hec1); %Faster to do all the random numbres upfront
Hec1_Positions=zeros(3*Number_Hec1,Number_Time_Steps);
Hec1_Positions(1:3:Number_Hec1*3)=-18;%Location of each Hec1 molecule. All start at (0,0,0)
Hec1_Bound=zeros(1,Number_Hec1); %To be used to calculate the fraction bound. Zero is unbound, one is bound

%Pull the microtubule bending dynamics
params.t_steps=Number_Time_Steps
params.n_dimers=20;
Microtubule_Position=microtubule(params);
Microtubule_Position(:,2,:)=20*Microtubule_Position(:,2,:)/max(max(max(Microtubule_Position(:,2,:))));

% Let Hec1 Diffuse
for Time=2:Number_Time_Steps
    for Hec1=1:Number_Hec1
        %First look at unbound Hec1
        if Hec1_Bound(Hec1)==0
            %Let the free Hec1 diffuse
            if RandSteps(3*Time*Number_Hec1 - Hec1 *3 - 2)<0.5
                Hec1_Positions(3*Hec1-2,Time)=Hec1_Positions(3*Hec1-2,Time-1) + 1;
                else
                Hec1_Positions(3*Hec1-2,Time)=Hec1_Positions(3*Hec1-2,Time-1) - 1;
            end

            if RandSteps(3*Time*Number_Hec1 - Hec1 *3 - 1)<0.5
                Hec1_Positions(3*Hec1-1,Time)=Hec1_Positions(3*Hec1-1,Time-1) + 1;
                else
                Hec1_Positions(3*Hec1-1,Time)=Hec1_Positions(3*Hec1-1,Time-1) - 1;
            end

            if RandSteps(3*Time*Number_Hec1 - Hec1 *3)<0.5
                Hec1_Positions(3*Hec1,Time)=Hec1_Positions(3*Hec1,Time-1) + 1;
                else
                Hec1_Positions(3*Hec1,Time)=Hec1_Positions(3*Hec1,Time-1) - 1;
            end
            
            %If the Hec1 is outside the tether, bring it back to where it
            %was the timestep before
            if ((Hec1_Positions(3*Hec1,Time)-Hec1_Positions(3*Hec1,1))^2 + (Hec1_Positions(3*Hec1 -1,Time)-Hec1_Positions(3*Hec1-1,1))^2 + (Hec1_Positions(3*Hec1-2,Time)-Hec1_Positions(3*Hec1,1))^2) > Tether_Length^2
                Hec1_Positions(3*Hec1-2,Time)=Hec1_Positions(3*Hec1-2,Time-1);
                Hec1_Positions(3*Hec1-1,Time)=Hec1_Positions(3*Hec1-1,Time-1);
                Hec1_Positions(3*Hec1,Time)=Hec1_Positions(3*Hec1,Time-1);
            else
            end
                
            %Find the minimum distance to the microtubule
            Distance=sqrt(min((Hec1_Positions(3*Hec1-2,Time)-Microtubule_Position(:,1,Time)).^2 + ...
                (Hec1_Positions(3*Hec1-1,Time)-Microtubule_Position(:,2,Time)).^2 + ...
                (Hec1_Positions(3*Hec1,Time)).^2));
            
            %Determine if this Hec1 will bind
             if Distance < Binding_Distance
                 if rand(1) < Prob_Bind
                 Hec1_Bound(Hec1)=1;
                 else
                 end
             else
             end
             
        % Now determine if the bound Hec1 will fall off
        else 
            Hec1_Positions(3*Hec1-2,Time)=Hec1_Positions(3*Hec1-2,Time-1);
            Hec1_Positions(3*Hec1-1,Time)=Hec1_Positions(3*Hec1-1,Time-1);
            Hec1_Positions(3*Hec1,Time)=Hec1_Positions(3*Hec1,Time-1);
            if rand(1) < Prob_Fall
                Hec1_Bound(Hec1)=0;
            else
            end
        end
    end
  
     %Compute the fraction bound
     Fraction_Bound(Time)=sum(Hec1_Bound)/Number_Hec1;
end

% Plot the trajectories in the z direction (adding x and y makes the graph
% harder to read)
plot(Hec1_Positions(1:3:end,:)')
title('Hec1 Positions')
figure(2)
plot(Fraction_Bound)
title('Fraction Bound')
figure(3)
for j=1:10
plot(Microtubule_Position(:,1,1000*j),Microtubule_Position(:,2,1000*j))
hold on
end
title('Microtubule Position Every 1000 Timesteps')