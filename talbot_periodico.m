% ------ DEFINICIÓN DE VARIBLES ------

n = int64(1799); % Puntos verticales
m = int64(399); % Puntos horizontales
t = 0; % Tiempo
A = 1; % Amplitud de la onda

N = 20; % "Discretizacion" espacial. Ver h.
M = sqrt(2)*N; % "Discretizacion" temporal. Debe ser mayor que N para que sea estable.
l = 10; % Anchura de las rendijas
w = pi*12*10^14; % Frecuencia angular de la onda para lambda=500 nm
tau = 5/3*10^(-15)/M; % Tiempo discreto, f/M
lambda = 5*10^(-7); % Longitud de onda
h = lambda/N; % Longitud del espacio discreto (cada "elemento" de la matriz representa un area de hxh)
r = (N/M)^2; % Debe ser menor que 1 para asegurar la estabilidad del algoritmo

R = 3; % Número de rendijas

x = double((0:m))*double(h); % Posiciones en X
y = double((0:n))*double(h); % Posiciones en Y
contador = 1; % Contador de pasos
modulo = 2; % Cada cuantos pasos se toma captura





% -------- EMPIEZA EL PROGRAMA -----------


% Comprobamos que r es menor que 1
disp(r);

% Mostramos la distancia de Talbot
disp(['Zt = ',num2str(5*10^(-7)/(1-sqrt(1-(N/(double(m+1)/R))^2)))] );
disp(double(m+1)/R - double(2*l));
disp(['Distancia entre rendijas = ',num2str(h*(double(m+1)/R - double(2*l)))]);
disp(['Distancia horizontal = ',num2str(h*double(n))]);

% Defino las matrices con los datos iniciales
u_ahora = zeros(n+1,m+1);
u_despues = zeros(n+1,m+1);

% Calculamos u en t=1
t = tau;

% Ponemos la fuente
u_despues(1,1:m+1) = A*sin(w*t);

% Ponemos la pared y las rendijas
u_despues(1,1:1/(2*R)*(m+1)-l) = 0;
for c=1:R-1
    u_despues(1,(2*c-1)/(2*R)*(m+1)+l:(2*c+1)/(2*R)*(m+1)-l) = 0;
end
u_despues(1,(2*R-1)/(2*R)*(m+1)+l:m+1) = 0;

for i = 2:n
    % Actualizamos cada punto
    for j = 2:m
        u_despues(i,j) = 0.5*r*(u_ahora(i+1,j) + u_ahora(i,j+1) + u_ahora(i-1,j) + u_ahora(i,j-1)) + (1-2*r)*u_ahora(i,j);
    end
end
% Derivada nula en los extremos
u_despues(1:n+1,1) = u_despues(1:n+1,2);
u_despues(1:n+1,m+1) = u_despues(1:n+1,m);
    
% Consideramos que apenas hay cambio al fondo
u_despues(n+1,1:m+1) = u_despues(n,1:m+1);

% Actualizamos las matrices
u_antes = u_ahora;
u_ahora = u_despues;

% Generamos la imagen
%contour(x,y,u_ahora,10);
fig = figure;
fig.WindowState = 'maximized';
colormap('gray');   % set colormap
imagesc(y,x,u_ahora.'.^2, [0 A^2]);        % draw image and scale colormap to values range
colorbar; % set colorbar
drawnow
frame = getframe(fig);
im = frame2im(frame);
[imind,map] = rgb2ind(im,256);
imwrite(imind,map,['prueba',num2str(contador),'.jpg'])

% Calculamos el estado del sistema en el tiempo t = k*l*tau
for repeticiones = 1:100
    disp(repeticiones)
    for k = 1:M
        contador = contador + 1;
        t = t + tau;
        
        % Ponemos la fuente en el origen
        u_despues(1,1:m+1) = A*sin(w*t);
        
        % Ponemos la pared y las rendijas
        u_despues(1,1:1/(2*R)*(m+1)-l) = 0;
        for c=1:R-1
            u_despues(1,(2*c-1)/(2*R)*(m+1)+l:(2*c+1)/(2*R)*(m+1)-l) = 0;
        end
        u_despues(1,(2*R-1)/(2*R)*(m+1)+l:m+1) = 0;
        
        % Actualizamos cada punto
        for i = 2:n
            for j = 2:m
                u_despues(i,j) = r*(u_ahora(i+1,j) + u_ahora(i,j+1) + u_ahora(i-1,j) + u_ahora(i,j-1)) + (2-4*r)*u_ahora(i,j) - u_antes(i,j);
            end
        end
        % Derivada nula en los extremos
        u_despues(1:n+1,1) = u_despues(1:n+1,2);
        u_despues(1:n+1,m+1) = u_despues(1:n+1,m);

        % Consideramos que en el fondo vale 0
        u_despues(n+1,1:m+1) = 0;
        
        % Actualizamos las matrices
        u_antes = u_ahora;
        u_ahora = u_despues;
        
        
        % Generamos la imagen
        
        imagesc(y,x,u_ahora.'.^2, [0 A^2]);        % draw image and scale colormap to values range
        colorbar; % set colorbar
        drawnow
        frame = getframe(fig);
        im = frame2im(frame);
        [imind,map] = rgb2ind(im,256);
        imwrite(imind,map,['prueba',num2str(contador),'.jpg'])
    end
end



% Create a video writer object
writerObj = VideoWriter('prueba.avi');

% Set frame rate
writerObj.FrameRate = 30;

% Open video writer object and write frames sequentially
open(writerObj)

for i = 1:contador                   % Some number of frames
     % Read frame
     frame = sprintf('prueba%d.jpg', i);
     input = imread(frame);

     % Write frame now
     writeVideo(writerObj, input);
end

% Close the video writer object
close(writerObj);

% 'Video.avi' will be created in the folder that contains the code.

disp("¡Terminado!");

