function main()
close all;
clear;
clc;


velocity = [-15,0]; % 20 m/s
omega = 10 * 2 * pi; % 10 Hz
radius = 0.15; % 15 cm

% PART ONE, get file (only works with two blades boomerangs)
[file,path] = uigetfile('*.jpg');
if isequal(file,0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(path,file)]);
end

% PART TWO, read the jpg and do image processing
I = imread(fullfile(path,file));
I = rgb2gray(I); % convert to grayscale
BW = imbinarize(I,0.9); % can change 0.9 to different threshold
BW = imcomplement(BW);
BW = imfill(BW,'holes');
B = bwboundaries(BW,'noholes');
boomerang = B{1}';

% PART THREE: process the matrices of coordinates and get inner and outer
% curve

% shift origin to center of mass
[r_mag, c] = find(BW == 1);
boomerang(1,:) = boomerang(1,:) - mean(r_mag);
boomerang(2,:) = boomerang(2,:) - mean(c);

NUM_POINTS= length(boomerang);

if (mod(NUM_POINTS,2) == 0)
    outer = boomerang(:,1:floor(NUM_POINTS/2));
    inner = boomerang(:, (floor(NUM_POINTS/2)+1):NUM_POINTS);
else
    % NUM_POINTS is odd so divide like so
    outer = boomerang(:,1:floor(NUM_POINTS/2) + 1);
    inner = boomerang(:, (floor(NUM_POINTS/2)+1):NUM_POINTS);
end



% sort inner and outer based on y value
inner = sortrows(inner',2)';
outer = sortrows(outer',2)';


% PART FOUR: fit curves and find blade segment

% figure 1 set up
figure(1);

% plot origin
plot(0,0,'r.','markersize',15);plot(0,0,'b.','markersize',10);

% plot inner and outer curves
subplot(1,2,1);
hold on;
axis equal;
xlim([-500 500]);
ylim([-500 500]);
plot(0,0,'r.','markersize',15);plot(0,0,'b.','markersize',10);
plot(inner(2,:), inner(1,:));
plot(outer(2,:), outer(1,:));


% fit inner and outer, and middle curves
f_inner_object=fit(inner(2,:)',inner(1,:)','fourier7');
f_inner = @(x) f_inner_object(x);
f_outer_object=fit(outer(2,:)',outer(1,:)','fourier7');
f_outer = @(x) f_outer_object(x);
middle = (inner + outer)/2;
f_middle=fit(middle(2,:)',middle(1,:)','fourier7');

plot(middle(2,:), middle(1,:));

subplot(1,2,2);
hold on;
axis equal;
xlim([-500 500]);
ylim([-500 500]);
plot(0,0,'r.','markersize',15);plot(0,0,'b.','markersize',10);
plot(f_middle);
plot(f_inner_object);
plot(f_outer_object);

%skip some outlier points at two ends
NUM_OF_POINTS = 1000; % can change this variable if wanting more or less elements
CUT_OFF_POINTS = 0;
start_point = middle(:,CUT_OFF_POINTS + 1);
end_point = middle(:,length(middle)-CUT_OFF_POINTS);

points_to_calculate = [linspace(start_point(2), end_point(2), NUM_OF_POINTS);1:NUM_OF_POINTS]; % second row is place holder
inner_points = zeros(2,NUM_OF_POINTS);
outer_points = zeros(2,NUM_OF_POINTS);

for i = 1: NUM_OF_POINTS
    points_to_calculate(2,i) = f_middle(points_to_calculate(1,i));
    slope = differentiate (f_middle, points_to_calculate(1,i));
    f_normal = @(x)  ((-1) / slope ) * (x - points_to_calculate(1,i)) + points_to_calculate(2,i);
    
    
    f_inner_middle = @(x) f_normal(x) - f_inner(x);
    xint = fzero(f_inner_middle, points_to_calculate(1,i));
    inner_points(:,i) = [xint; f_inner(xint)];
    
    
    f_outer_middle = @(x) f_normal(x) - f_outer(x);
    xint = fzero(f_outer_middle, points_to_calculate(1,i));
    outer_points(:,i) = [xint; f_outer(xint)];
end

%plot the disections
% figure(3);
% hold on;
% plot(inner(2,:), inner(1,:));
% plot(outer(2,:), outer(1,:));
% for i = 1: 10:NUM_OF_POINTS
%     plot([inner_points(1,i); outer_points(1,i)],[inner_points(2,i); outer_points(2,i)]);
% end

inner_1 = inner_points(:,1: NUM_OF_POINTS/2);
inner_2 = inner_points(:, NUM_OF_POINTS/2 + 1: NUM_OF_POINTS);

outer_1 = outer_points(:,1: NUM_OF_POINTS/2);
outer_2 = outer_points(:, NUM_OF_POINTS/2 + 1: NUM_OF_POINTS);

% see where our edges are
% scatter(outer_2(1,:), outer_2(2,:),'g');
% scatter(outer_1(1,:), outer_1(2,:),'c');
% scatter(inner_1(1,:), inner_1(2,:),'r');
% scatter(inner_2(1,:), inner_2(2,:),'b');


%PART FIVE: scale boomerang and calculate lift
% scale the boomerang to match the radius

r_mag = sqrt(points_to_calculate(1,:) .* points_to_calculate(1,:) + points_to_calculate(2,:) .* points_to_calculate(2,:));
scalingfactor = radius / max(r_mag);
points_to_calculate = points_to_calculate * scalingfactor;
inner_points = inner_points * scalingfactor;
outer_points = outer_points * scalingfactor;
boomerang = boomerang * scalingfactor;

inner_1 = inner_1 * scalingfactor;
inner_2 = inner_2 * scalingfactor;
outer_1 = outer_1 * scalingfactor;
outer_2 = outer_2 * scalingfactor;



rot90 = [0 -1; 1 0];
n1 = [points_to_calculate(1,NUM_OF_POINTS/2) - points_to_calculate(1,1),points_to_calculate(2,NUM_OF_POINTS/2) - points_to_calculate(2,1)]; %vector for first wing --- V's left wing
n2 = [points_to_calculate(1,NUM_OF_POINTS) - points_to_calculate(1,NUM_OF_POINTS/2),points_to_calculate(2,NUM_OF_POINTS) - points_to_calculate(2,NUM_OF_POINTS/2)]; %vector for second wing --- V's right wing
n1 = n1 / sqrt(n1(1)*n1(1) + n1(2)*n1(2)); % unit vector
n2 = n2 / sqrt(n2(1)*n2(1) + n2(2)*n2(2)); % unit vector
n1_norm = (rot90 * n1')';
n2_norm = (rot90 * n2')';

r = zeros(NUM_OF_POINTS,1,2);
for i=1: NUM_OF_POINTS
    rw(i,:) = (rot90 * [points_to_calculate(1,i),points_to_calculate(2,i)]' * omega)';
    v_tot(i,:) = velocity + rw(i,:);
    if ( i < NUM_OF_POINTS/ 2 ) %left wing
        V_N(i,:) = dot(n1_norm, v_tot(i,:)) * n1_norm;
    else
        V_N(i,:) = dot(n2_norm, v_tot(i,:)) * n2_norm;
    end
    
end




% caculate lift based on middle points using average data for V and w
% r_mag = sqrt(points_to_calculate(1,:) .* points_to_calculate(1,:) + points_to_calculate(2,:) .* points_to_calculate(2,:));
% cos_phi = points_to_calculate(2,:) ./ r_mag;

lift = V_N(:,1) .* V_N(:,1) + V_N(:,2) .* V_N(:,2);


%PART SIX: display animation figure
figure(2);
figure('position', [70, 70, 900, 900])
%plot the origin and graph setup
hold on;
axis equal;
colorbar;
xlim([-500 500]);
ylim([-500 500]);
plot(0,0,'r.','markersize',15);
plot(0,0,'b.','markersize',10);
% switch x and y axis
boomerang([1 2],:)=boomerang([2 1],:);

boomerangHandle = fill(boomerang(1,:),boomerang(2,:),[0,0,0]+0.5);

%plot animation
%plot on correct side
toPlotCoord = zeros(2, NUM_OF_POINTS);
%             for i=length(inner_1)+1: length(lift)
%                 if (lift(i) < 0)
%                     toPlotCoord(:,i) = outer_2(:,i- length(inner_1));
%                 else
%                     toPlotCoord(:,i) = inner_2(:,i- length(inner_1));
%                 end
%             end
%
%             for i= 1: length(inner_1)
%                 if (lift(i) < 0)
%                     toPlotCoord(:,i) = inner_1(:,i );
%                 else
%                     toPlotCoord(:,i) = outer_1(:,i );
%                 end
%             end
for i = 1:NUM_OF_POINTS
    if (V_N(i,1) < 0)
        if inner_points(1,i) < outer_points(1,i)
            toPlotCoord(:,i) =  inner_points(:,i );
        else
            toPlotCoord(:,i) =  outer_points(:,i );
        end
    else
        if inner_points(1,i) > outer_points(1,i)
            toPlotCoord(:,i) =  inner_points(:,i );
        else
            toPlotCoord(:,i) =  outer_points(:,i );
        end
    end
end


correctSideHandle = scatter(toPlotCoord(1,:), toPlotCoord(2,:),30, lift,'filled');

xlim manual;
ylim manual;
xlim([-0.2 0.2]);
ylim([-0.2 0.2]);
colorbar;
caxis([-300 500]);
button = uicontrol('Style','togglebutton','String', 'stop/start');
button.Callback = @buttonPushed;
stop = true;

video_frame = 1;


    function buttonPushed(src,event)
        
        stop = ~stop;
        while(~stop)
            
            
            delete(boomerangHandle);
            delete(correctSideHandle);
            
            theta = pi /60;
            rotationFactor = [cos(theta), -sin(theta);
                sin(theta), cos(theta)];
            
            
            boomerang = rotationFactor * boomerang;
            middle = rotationFactor * middle;
            points_to_calculate = rotationFactor * points_to_calculate;
            inner_1 = rotationFactor * inner_1;
            inner_2 = rotationFactor * inner_2;
            outer_1 = rotationFactor * outer_1;
            outer_2 = rotationFactor * outer_2;
            inner_points = rotationFactor * inner_points;
            outer_points = rotationFactor * outer_points;
            n1_norm = (rotationFactor * n1_norm')';
            n2_norm = (rotationFactor * n2_norm')';
            
            
            r = zeros(NUM_OF_POINTS,1,2);
            for i=1: NUM_OF_POINTS
                rw(i,:) = (rot90 * [points_to_calculate(1,i),points_to_calculate(2,i)]' * omega)';
                v_tot(i,:) = velocity + rw(i,:);
                if ( i < NUM_OF_POINTS/ 2 ) %left wing
                    V_N(i,:) = dot(n1_norm, v_tot(i,:)) * n1_norm;
                else
                    V_N(i,:) = dot(n2_norm, v_tot(i,:)) * n2_norm;
                end
                
            end
            
            
            
            
            % caculate lift based on middle points using average data for V and w
            % r_mag = sqrt(points_to_calculate(1,:) .* points_to_calculate(1,:) + points_to_calculate(2,:) .* points_to_calculate(2,:));
            % cos_phi = points_to_calculate(2,:) ./ r_mag;
            
            lift = V_N(:,1) .* V_N(:,1) + V_N(:,2) .* V_N(:,2);
            r = zeros(NUM_OF_POINTS,1,2);
            for i=1: NUM_OF_POINTS
                rw(i,:) = (rot90 * [points_to_calculate(1,i),points_to_calculate(2,i)]' * omega)';
                v_tot(i,:) = velocity + rw(i,:);
                if ( i < NUM_OF_POINTS/ 2 ) %left wing
                    V_N(i,:) = dot(n1_norm, v_tot(i,:)) * n1_norm;
                else
                    V_N(i,:) = dot(n2_norm, v_tot(i,:)) * n2_norm;
                end
                
            end
            
            
            
            % caculate lift based on middle points using average data for V and w
            % r_mag = sqrt(points_to_calculate(1,:) .* points_to_calculate(1,:) + points_to_calculate(2,:) .* points_to_calculate(2,:));
            % cos_phi = points_to_calculate(2,:) ./ r_mag;
            
            lift = V_N(:,1) .* V_N(:,1) + V_N(:,2) .* V_N(:,2);
            colorbar;
            
            boomerangHandle = fill(boomerang(1,:),boomerang(2,:),[0,0,0]+0.5);
            
            %plot animation
            %plot on correct side
            toPlotCoord = zeros(2, NUM_OF_POINTS);
            %             for i=length(inner_1)+1: length(lift)
            %                 if (lift(i) < 0)
            %                     toPlotCoord(:,i) = outer_2(:,i- length(inner_1));
            %                 else
            %                     toPlotCoord(:,i) = inner_2(:,i- length(inner_1));
            %                 end
            %             end
            %
            %             for i= 1: length(inner_1)
            %                 if (lift(i) < 0)
            %                     toPlotCoord(:,i) = inner_1(:,i );
            %                 else
            %                     toPlotCoord(:,i) = outer_1(:,i );
            %                 end
            %             end
            for i = 1:NUM_OF_POINTS
                if (V_N(i,1) < 0)
                    if inner_points(1,i) < outer_points(1,i)
                        toPlotCoord(:,i) =  inner_points(:,i );
                    else
                        toPlotCoord(:,i) =  outer_points(:,i );
                    end
                else
                    if inner_points(1,i) > outer_points(1,i)
                        toPlotCoord(:,i) =  inner_points(:,i );
                    else
                        toPlotCoord(:,i) =  outer_points(:,i );
                    end
                end
            end
            
            
            
            correctSideHandle = scatter(toPlotCoord(1,:), toPlotCoord(2,:),30, lift,'filled');
            
            F(video_frame) = getframe(gcf) ;
            video_frame = video_frame + 1;
            pause(0.2);
            
        end
%         % create the video writer with 1 fps
%         writerObj = VideoWriter('myVideo');
%         writerObj.FrameRate = 10;
%         % set the seconds per image
%         % open the video writer
%         open(writerObj);
%         % write the frames to the video
%         for i=1:length(F)
%             % convert the image to a frame
%             frame = F(i) ;
%             writeVideo(writerObj, frame);
%         end
%         % close the writer object
%         close(writerObj);
        
    end




end
