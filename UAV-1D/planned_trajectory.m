function [ x_des ] = planned_trajectory(t,params)

%FIXED_SET_POINT  Outputs a constant desired state = [z_des;0] except at t = 0 where it returns [0;0]

if t == 0
  x_des = [-2;-2];
else
  x_des = [params.z_des;0];
end


end
