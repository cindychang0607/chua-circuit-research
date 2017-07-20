  function [stimpoints, jumps] = generateChuaPerturbations(istate, timestep, R, dist, stimmags,stimtime,iterations)
%stimmags = [xstim ystim zstim]
% dist: distance on attractor
% iterations: number of time distance is traveled and perturbated
% stimmags [vector] : magnitude of perturbation in [x y z] direction
% stimtime: number of steps for each perturbation
% OUTPUT
% stimpoints: vector of coordinates where purturbation occured
% jumps: 1 if jumped to other attractor, 2 if jumped to outer state, 0 is
% it stayed on the same attractor

stimpoints = zeros(iterations,3);
jumps = zeros(iterations,6);

distance = dist;

%%% initial transient
[~, nextstate] =  Chua(istate, timestep, R, distance+dist*iterations, 0, 2,[0 0 0],2000);

%%% perturbates in all 6 directions
for i = 1:iterations
    istate = nextstate;
    
    [stimpoints(i,:), jumps(i,1)] = Chuap(istate,timestep ,R ,distance, stimmags(1), stimtime, [1 0 0]);
    
    [stimpoints(i,:), jumps(i,2)] = Chuap(istate,timestep ,R ,distance,-stimmags(1), stimtime, [1 0 0]);
    
    [stimpoints(i,:), jumps(i,3)] = Chuap(istate,timestep ,R ,distance, stimmags(2), stimtime, [0 1 0]);
    
    [stimpoints(i,:), jumps(i,4)] = Chuap(istate,timestep ,R ,distance,-stimmags(2), stimtime, [0 1 0]);
    
    [stimpoints(i,:), jumps(i,5)] = Chuap(istate,timestep ,R ,distance, stimmags(3), stimtime, [0 0 1]);
    
    [stimpoints(i,:), jumps(i,6)] = Chuap(istate,timestep ,R ,distance,-stimmags(3), stimtime, [0 0 1]);
    
    %iteration_num = i; %#ok<NOPRT,NASGU>
    
    nextstate = stimpoints(i,:)';
end
end    
    
    
