function drawBallnBeam(X_history, Theta_history)
    %Creates an animation for the Ball n Beam Control Project.
    % Inputs:
    %       X_history       : Time history of the ball position on the beam
    %       Theta_history   : Time history of the beam angle
    % Example Usage:
    %       >> t = linspace(0,100,500);
    %       >> drawBallnBeam(sin(t)/4,cos(t)/5);

    % Author: berkguler20@ku.edu.tr
    % Koc Uni - Mech 304 Term Project
    % Lecturer: Prof. Cagatay Basdogan
    % May 2022; Last revision: 23-May-2022


    assert(length(X_history) == length(Theta_history), "Data Lengths are not the same !!!");
    % === Configurations ===
    beam_length = 0.6; %60cm
    beam_width = 0.011; % 1.1 cm
    ball_radius = 0.05; %5cm
    % === Figure Properties ===;
    set(gca,'fontname','times');
    axis square;
    xlim([-1 1]);
    xlabel("X [m]");
    ylim([-1 1]);
    ylabel("Y [m]");
    title("MECH 304 | The Ball & Beam Project")
    hold on;
    grid on;
    
    beam_x = ([-beam_length/2, beam_length/2, beam_length/2, -beam_length/2]);
    beam_y = ([-beam_width/2, -beam_width/2, beam_width/2, beam_width/2]);
    theta = 0;
    X = 0;
    HTM = zeros(2,4);
    beam_rotation = ([cos(theta) -sin(theta); sin(theta) cos(theta)]);
    for k = 1:4
        HTM(:,k) = beam_rotation*[beam_x(k); beam_y(k)];
    end
    beam_plt = fill(HTM(1,:), HTM(2,:), "k");
    ball_plt = rectangle('Position', [-ball_radius - X*cos(theta), beam_width/2-X*sin(theta), ball_radius*2, ball_radius*2],"Curvature",[1 1],'EdgeColor',"r",'FaceColor', 'r');
    
    for i = 1:length(X_history)
        delete(beam_plt);
        delete(ball_plt);
        theta = Theta_history(i);
        X = X_history(i);
        beam_rotation = ([cos(theta) -sin(theta); sin(theta) cos(theta)]);
        for k = 1:4
            HTM(:,k) = beam_rotation*[beam_x(k); beam_y(k)];
        end

    
        beam_plt = fill(HTM(1,:), HTM(2,:), "k");
        ball_plt = rectangle('Position', [-ball_radius - X*cos(theta), beam_width/2-X*sin(theta), ball_radius*2, ball_radius*2],"Curvature",[1 1],'EdgeColor',"r",'FaceColor', 'r');
        pause(0.02);
    end
end