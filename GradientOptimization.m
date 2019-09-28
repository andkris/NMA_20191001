

function GradientOptimization
    % clear all; 
    close all; clc; close all;

    % preparing data
    n = 10;  optNodes = 6;
    positions = randi([-5 5], 2*n,1);
    searchedPositions = randi([-5 5], 2*optNodes,1);


    % preparing figure for visualization
    figure(1); hold on; grid on; axis([-5 5 -5 5]);
    plot(positions(1:2:end),positions(2:2:end),'ro',...
     'MarkerFaceColor',[1 0 0],'MarkerSize',8);
    plot(searchedPositions(1:2:end), searchedPositions(2:2:end),'go',...
    'MarkerSize',6);
    h = [];

    % Gradient params
    step = 0.1;   eps = 1e-9;   maxIt = 1000;
    x0 = searchedPositions;
    f0 =  objectiveF(x0, positions, n, optNodes);

    prec = 1;
    for i = 1:maxIt
      grad = quasiGradient(x0, positions, n,optNodes, step/1000);
      if norm(grad) < eps % if we do not find 
        fprintf(1,'gradientas %6.5f\n',norm(grad));
        break;
      end
      grad = grad/norm(grad);
      x1 = x0-step*grad;
      f1 =  objectiveF(x1, positions, n,optNodes);
      if(f1>f0)
          step = step / 10;
      else
          prec = norm(x1-x0)/norm(x1+x0);
          f0 = f1; 
          x0 = x1;
      end
      fprintf(1,'iteration %d, objective function %6.5f, step %6.8f, prec %6.8f\n', i,f0, step, prec);

      % visualization
      if ~isempty(h), 
          delete(h);
      end
      h = visualization(x0,optNodes,positions);
      pause(0.00001);
      if prec < eps
          break;
      end
    end
end

% Objective function
function f = objectiveF(x0, positions, n, m)

    x = positions(1:2:end);  y = positions(2:2:end);
    
    f = 0; 
    D = [];
    
    for i = 1:m 
        xi = x0(2*i-1); 
        yi = x0(2*i); 
        D = [D; calculateDistances([xi, yi], [x y] )]; 
    end

    x0other = x0(1:2:end); 
    y0other = x0(2:2:end);
    for i = 1:m-1
        xi = x0(2*i-1); 
        yi = x0(2*i); 
        D = [D; calculateDistances([xi, yi], [x0other(i+1:end) y0other(i+1:end)] )];
    end
    f = sum((D-mean(D)).^2);
end


% distance between points
function D = calculateDistances(point, array)
    D = sqrt((point(1)-array(:,1)).^2+(point(2)-array(:,2)).^2);
end

% approximated derivatives
function df = quasiGradient(positions, givenPositions, n,optNodes, dx)
    df = zeros(size(positions)); 
    f0 = objectiveF(positions, givenPositions, n,optNodes); 
    for i = 1:numel(positions)
        XNew = positions; 
        XNew(i) = XNew(i)+dx;
        f1 = objectiveF(XNew, givenPositions, n,optNodes);
        df(i) = (f1-f0)/dx;
    end

end



function h = visualization(positions,n,givenPoints)
    x = positions(1:2:end);     y = positions(2:2:end);
    h = [];
    for i = 1:n-1
        for j = i+1:n
            h = [h; plot([x(i) x(j)],[y(i) y(j)],'Color',[0.5 0.5 0.5])];
        end
    end
    xGiven = givenPoints(1:2:end); 
    yGiven = givenPoints(2:2:end); 
    for i = 1:n
        for j = 1:numel(xGiven)
            h = [h; plot([x(i) xGiven(j)],[y(i) yGiven(j)],'Color',[0.5 0.5 0.5])];
        end
    end
    h = [h;  plot(x(1:end),y(1:end),'go','MarkerFaceColor',[0 0.5 1],'MarkerSize',8);];
end