function [] = partII()

    % generate the data

    rng(1); 
    r = sqrt(rand(100,1)); 
    t = 2*pi*rand(100,1);  
    data1 = [r.*cos(t), r.*sin(t)]; 

    r2 = sqrt(3*rand(100,1)+1); 
    t2 = 2*pi*rand(100,1);      
    data2 = [r2.*cos(t2), r2.*sin(t2)]; 

    % plot the data

    figure;
    plot(data1(:,1),data1(:,2),'r.','MarkerSize',15)
    hold on
    plot(data2(:,1),data2(:,2),'b.','MarkerSize',15)
    axis equal
    hold on

    % work on class 1
    [a1, R1] = calcRandCentre(data1);

    % work on class 2
    [a2, R2] = calcRandCentre(data2);

    % plot centre and radius for class 1
    plot(a1(1), a1(2), 'rx', 'MarkerSize', 15);
    viscircles(a1', R1, 'Color', 'r', 'LineWidth', 1);
    hold on

    % plot centre and radius for class 2
    plot(a2(1), a2(2), 'bx', 'MarkerSize', 15);
    viscircles(a2', R2, 'Color', 'b', 'LineWidth', 1);

end

function [a, R] = calcRandCentre(data)

    C = 0.45;

    H = data*data';
    
    f = -diag(H);
    A = ones(1,length(data));
    c=ones(1,1);
    g_l = zeros(1,length(data));
    g_u = C*ones(1,length(data));

    g = quadprog (H, f, A,c,A,c,g_l, g_u);
    a = (g'*data)';
    
    R_squared = 0;
    counter = 0;
    epsilon = 0.0000000001;
    for i = 1:length(data)
        if g(i) > epsilon & g(i) < C - epsilon
            R_squared = R_squared + (data(i)-a)'*(data(i)-a);
            counter = counter+1;
        end
    end
    R = sqrt(R_squared/counter);
end