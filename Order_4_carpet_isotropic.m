% Program made by Gabriel María Ybarra Marcaida 
% under Luis Vega González's supervision.

% Memory clear (just in case)
clear
clearvars;
clearvars -global;

% ------ DEFINITION OF THE VARIABLES ------

% PHYSICAL PARAMETERS

A = 1;                  % Wave's amplitude
c = 1;                  % Speed of light
lambda = 0.1;           % Wavelength
omega = 2*pi*c / lambda;% Angular frequency
w = 0.5*lambda;         % Grating width.
d = 1;                  % Grating separation.
zt = 2*d^2/lambda;      % Talbot distance.


% CANVAS PARAMETERS

X_max = d/2;
Z_max = zt;

extra_steps = 1.52;     % The total number of steps is zt*extra_steps/h_z

% DISCRETIZATION PARAMETERS

N = 30;                             % Space discretization, with h = λ/N.
M = sqrt(2)*N_z;                    % Time discretization, with τ = f/M.
r = ( 1/M / (1/N + 1/N_z) )^2;      % Must be smaller than 1/2 if we want a stable algorithm.
a_discrete = ceil( w/2/lambda*N );  % Grating's discrete width.

tau = lambda/c/M;       % Temporal step, τ = f/M.
h = lambda/N;           % Spatial X step (every matrix element represents an area of hxh). hX = λ/N.

n_x = ceil( X_max/h );  % Number of points in the x direction.
n_z = ceil( Z_max/h );  % Number of points in the z direction.



% OTHER VARIABLES

t = 0;                                      % Current time
x = double((-3*X_max:double(h):3*X_max));   % Matrix with the x-axis positions
z = double((0:double(h):Z_max));            % Matrix with the z-axis positions


tstep = 1;                                  % tstep counter
modulo = 3;                                 % Number of tsteps between photos






% -------- BEGINING OF THE MAIN PROGRAM -----------


% We check that r is smaller than 1/2
if r > 1/2
    disp("r = " + num2str(r) + " which is less than 1/2!");
    disp("The algorithm isn't stable!");
    return;
end



% Definition of the inicial matrices
u_ahora = zeros(n_x + 5,floor(n_z*extra_steps) + 4);
u_despues = u_ahora;
u2 = zeros(n_x+1,n_z+1);

% CALCULATION OF U AT N=1
t = tau;

% WE SET THE SOURCE
u_despues(3:a_discrete+3,1) = A*sin(omega*(t-tau));
u_despues(3:a_discrete+3,2) = A*sin(omega*t);

% WE IMPOSE PERIODICITY
u_despues(1,:) = u_despues(5,:);
u_despues(2,:) = u_despues(4,:);

u_despues(end-4,:) = u_despues(end,:);
u_despues(end-3,:) = u_despues(end-1,:);

% We update the matrices
u_antes = u_ahora;
u_ahora = u_despues;

for j=1:tstep
    u2(:,j) = u_ahora(3:end-2,j).^2;
end


% WE GENERATE THE IMAGE

fig = figure;
fig.WindowState = 'maximized';

set(gca,'FontSize',20)


colormap('turbo');      % set colormap
colorbar;               % set colorbar

foldername = ['images_d=',num2str(d),'_w=',num2str(w),'_l=',num2str(lambda),'_N=',num2str(N_z)];
destdirectory = fullfile(pwd,foldername);
mkdir(destdirectory);

generate_image(x,z,tstep,u2,A,fig,destdirectory);




% Calculamos el estado del sistema en el tiempo t = k*l*tau
    
    for tstep = 2:n_z
        disp("Step " + num2str(tstep) + " of " + num2str(n_z*extra_steps) + ".");
        t = t + tau;
        
        
        % WE SET THE SOURCE
        u_despues(3:a_discrete+3,1) = A*sin(omega*(t-tau));
        u_despues(3:a_discrete+3,2) = A*sin(omega*t);
        
        
        % We update every point
        for m = 3:tstep
            for l = 3:n_x+3
                u_despues(l,m) = 2*u_ahora(l,m) - u_antes(l,m) +...
                    (c*tau)^4 /12 * laplacian_squared(u_ahora,l,m,h) +...
                    (c*tau)^2 * laplacian(u_ahora,l,m,h);
            end
        end
        
        % WE IMPOSE PERIODICITY
        u_despues(1,1:tstep) = u_despues(5,1:tstep);
        u_despues(2,1:tstep) = u_despues(4,1:tstep);

        u_despues(end,1:tstep) = u_despues(end-4,1:tstep);
        u_despues(end-1,1:tstep) = u_despues(end-3,1:tstep);
        
        
        % We update the matrices
        u_antes = u_ahora;
        u_ahora = u_despues;
        

        
        % We generate the image
        if mod(tstep,modulo) == 1 || modulo==1

            % We compute the intensity
            for j=1:n_z+1
                u2(:,j) = u_ahora(3:end-2,j+1).^2;
            end

            generate_image(x,z,tstep,u2,A,fig,destdirectory);
        end
        
    end

    for tstep = n_z+1:floor(n_z*extra_steps)
        disp("Step " + num2str(tstep) + " of " + num2str(n_z*extra_steps) + ".");
        t = t + tau;
        
        
        % WE SET THE SOURCE
        u_despues(3:a_discrete+3,1) = A*sin(omega*(t-tau));
        u_despues(3:a_discrete+3,2) = A*sin(omega*t);
        
        
        % We update every point
        for m = 3:tstep
            for l = 3:n_x+3
                u_despues(l,m) = 2*u_ahora(l,m) - u_antes(l,m) +...
                    (c*tau)^4 /12 * laplacian_squared(u_ahora,l,m,h) +...
                    (c*tau)^2 * laplacian(u_ahora,l,m,h);
            end
        end
        
        % WE IMPOSE PERIODICITY
        u_despues(1,1:tstep) = u_despues(5,1:tstep);
        u_despues(2,1:tstep) = u_despues(4,1:tstep);

        u_despues(end,1:tstep) = u_despues(end-4,1:tstep);
        u_despues(end-1,1:tstep) = u_despues(end-3,1:tstep);
        
        
        % We update the matrices
        u_antes = u_ahora;
        u_ahora = u_despues;
        

        % We generate the image
        if mod(tstep,modulo) == 1 || modulo==1

            % We compute the intensity
            for j=1:n_z+1
                u2(:,j) = u_ahora(3:end-2,j+1).^2;
            end

            generate_image(x,z,tstep,u2,A,fig,destdirectory);
        end

    end

    

% We generate the video
generate_video(n_z*extra_steps,modulo,destdirectory);



fprintf("\n\nProgram made by Gabriel María Ybarra Marcaida under the supervision of Luis Vega González.\n");
return;




% --------- FUNCTIONS ---------


% Computes the laplacian of u at (l,m)
% A fourth order scheme is used
function lapl = laplacian(u,l,m,h)

    lapl = -5 * u(l,m) + ...
           4*( u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1))/3  + ...
           - ( u(l+2,m) + u(l-2,m) + u(l,m+2) + u(l,m-2) )/12;
    lapl = lapl / h^2;
end

% Computes the squared laplacian of u at (l,m)
% A second order scheme is used
function lapl2 = laplacian_squared(u,l,m,h)
    
    lapl2 = ( u(l-2,m) - 4*u(l-1,m) + 6*u(l,m) - 4*u(l+1,m) + u(l+2,m) )  +...
            ( 4*u(l,m) + ( u(l+1,m+1) + u(l+1,m-1) + u(l-1,m+1) + u(l-1,m-1) ) - 2*( u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1) ) );
    lapl2 = 2*lapl2 / h^4; 
end


% Generates the image

function generate_image(x,z,tstep,u2,A,fig,destdirectory)

    imagen = [flip(u2(3:end-3,:)); u2(3:end-2,:); flip(u2(4:end-3,:)); u2(3:end-2,:); flip(u2(4:end-3,:)); u2(3:end-2,:)];

    imagesc(x,z,imagen, [0 A^2]);        % draw image and scale colormap to values range
    colorbar; % set colorbar

    set(gca,'FontSize',20)

    Z = length(z)-1;
    yticks([0 z(floor(Z/6)) z(floor(Z/3)) z(floor(Z/2)) z(floor(2*Z/3)) z(floor(5*Z/6)) z(Z)])
    yticklabels({'-3/2 d','- d','-1/2 d','0','1/2 d','d','3/2 d'})
    
    
    X = length(x)-1;
    xticks([x(1) x(floor(X/4)) x(floor(X/2)) x(floor(3*X/4)) x(X)])
    xticklabels({'0','1/4 Z_T','1/2 Z_T','3/4 Z_T','Z_T'})

    drawnow
    frame = getframe(fig);
    im = frame2im(frame);
    thisimage = ['img_',num2str(tstep),'.png'];
    fulldestination = fullfile(destdirectory,thisimage);
    imwrite(im,fulldestination);
end

function generate_video(tstep,modulo,destdirectory)
    % Create a video writer object
    writerObj = VideoWriter(fullfile(destdirectory,'aa'));

    % Set frame rate
    writerObj.FrameRate = 30;

    % Open video writer object and write frames sequentially
    open(writerObj)

    for i = 1:modulo:tstep                   % Some number of frames
         % Read frame
         frame = sprintf('img_%d.png', i);
         frame = fullfile(destdirectory, frame);
         input = imread(frame);

         % Write frame now
         writeVideo(writerObj, input);
    end

    % Close the video writer object
    close(writerObj);

    % 'Video.avi' will be created in the folder that contains the code.
end

