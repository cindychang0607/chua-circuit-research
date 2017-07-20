function plotData(istate, timestep, R, dist, iterations,jumps, stimpoints, pltpoints)

%%%% plotting %%%%
for i = pltpoints
    figure(i)
    [attractor, stim] = Chua(istate, timestep, R, dist, 0, 2,[0 0 0], ceil(iterations/.05));
    plot3(attractor(1,:),attractor(2,:),attractor(3,:));
    hold on
    plot3(stim(1,:),stim(2,:),stim(3,:),'yx');
    for j = 1:iterations
        if jumps(j,i) == 1
            plot3(stimpoints(j,1),stimpoints(j,2),stimpoints(j,3),'go')
        elseif jumps(j,i) == 2
            plot3(stimpoints(j,1),stimpoints(j,2),stimpoints(j,3),'mo')
        else
            plot3(stimpoints(j,1),stimpoints(j,2),stimpoints(j,3),'ro')
        end
    end
end