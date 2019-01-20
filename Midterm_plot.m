%% Plot the shape of view for each of the angles
clear all;

I_0deg = [0,2;    1,2;    1,0;    2,3;    3,1;    4,2];

fig = figure;
plot(I_0deg(:,1),I_0deg(:,2));
xlabel('Length(s)');
ylabel('Accumulated Intensity');
grid on;
fig.Color = [1,1,1];


I_90deg = [0,0;
    1,3;
    2,0;
    3,1;
    3,3;
    4,4];
fig = figure;
plot(I_90deg(:,1),I_90deg(:,2));
xlabel('Length(s)');
ylabel('Accumulated Intensity');
grid on;
fig.Color = [1,1,1];

I_45deg = [0,0;    
    sqrt(2),0;
    sqrt(2),2;
    1.5*sqrt(2),2;
    2*sqrt(2),2+sqrt(2)*2;
    2*sqrt(2),sqrt(2)*2;
    2.5*sqrt(2),0;
    3*sqrt(2),0;
    3*sqrt(2),2*sqrt(2);
    4*sqrt(2),0];
fig = figure;
plot(I_45deg(:,1),I_45deg(:,2));
xlabel('Length(s)');
ylabel('Accumulated Intensity');
grid on;
fig.Color = [1,1,1];
