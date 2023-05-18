% Program made by Gabriel María Ybarra Marcaida 
% under Luis Vega González's supervision.

%Constants
lambda = [0.01];    % Wavelength. Several values are accepted, and the carpet will be plotted for each case.
w = 5*lambda;       % Grating width
d = 1;              % Distance in between gratings
zt = lambda./(1 - sqrt(1 - (lambda./d).^2)); % Talbot distance


pixel_x = 1080*2;   % Number of pixels in the X direction
pixel_z = 1920*2;   % Number of pixels in the X direction

destdirectory = fullfile(pwd,'images');    % Destination directory
mkdir(destdirectory);

fig = figure();

for k = 1:length(lambda)

    X_max = 1.5*d;  % Sets the X range of the image so there are three slits in display
    Z_max = zt(k);  % Sets the Z range of the image so there is only one Talbot plane

    escala_x = X_max/(pixel_x*2);
    escala_z = Z_max/(pixel_z);
    
    X = -X_max:escala_x:X_max;
    Z = 0:escala_z:Z_max;
    intensidades = zeros(length(X),length(Z));
    

    % We compute the intensity at each pizel

    for i = 1:length(X)
        for j = 1:length(Z)
            intensidades(i,j) = abs(T(X(i), Z(j), w(k), d, lambda(k)))^2;
        end
    end
    

    % We generate the image

    fig.Position = [0 0 pixel_z+1 2*pixel_x+1];    
    s = surf(Z,X,intensidades);
    fig.WindowState = 'maximized';

    
    set(gca,'FontSize',20)
    xticks([0 Z_max/4 Z_max/2 3*Z_max/4 Z_max])
    xticklabels({'0','1/4 Z_T','1/2 Z_T','3/4 Z_T','Z_T'})
    
    yticks([-X_max -2*X_max/3 -X_max/3 0 X_max/3 2*X_max/3 X_max])
    yticklabels({'3/2 d','d','1/2 d','0','-1/2 d','-d','-3/2 d'})
    set(gca,'XLim',[0 Z_max],'YLim',[-X_max X_max])
    
    colormap("turbo");
    colorbar;
    view(2);
    s.EdgeColor = 'none';
    s.LineStyle = 'none';
    s.FaceColor = 'interp';
    drawnow;

    frame = getframe(fig);
    im = frame2im(frame);
    thisimage = ['d_λ=',num2str(d/lambda(k)),'_w_λ=',num2str(w(k)/lambda(k)),'_carpet.png'];
    fulldestination = fullfile(destdirectory,thisimage);
    imwrite(im,map,fulldestination);

end

fprintf("\n\nProgram made by Gabriel María Ybarra Marcaida under the supervision of Luis Vega González.\n");
return;




function amplitud = T(x,z,w,d,lambda)
    
    amplitud = exp(-2*pi*1i/lambda*z);

    for n = 1:floor(d/lambda)
        amplitud = amplitud + 2 * sin(n*pi*w/d)/(n*pi*w/d) * exp(-2*pi*1i*sqrt(1/lambda^2 - (n/d)^2)*z) * cos(2*pi*n*x/d);
    end

    for n = ceil(d/lambda):11*ceil(d/lambda)
        amplitud = amplitud + 2 * sin(n*pi*w/d)/(n*pi*w/d) * exp(-2*pi*sqrt(-1/lambda^2 + (n/d)^2)*z) * cos(2*pi*n*x/d);
    end


    amplitud = amplitud*w/d; 
end
