function [R, error] = multicluster_linear_regression(A, n_trials, if_plot)
% ------------------------ no PCA------------------------------------------
% A = A';
% N = round(size(A,2)/n_trials);
% R = zeros(N,N);
% error = ones(N,N);
% for n = 1:N
%     for m = 1:N
%         if n==m
%             continue
%         else
%             downsample_idx = 1:10:size(A,1);
%             X = [A(downsample_idx,(n-1)*n_trials+1:n*n_trials), A(downsample_idx,(m-1)*n_trials+1:m*n_trials)];
%             if isempty(find(X~=0))
%                 continue
%             end
%             X = X';
%             Y = [ones(n_trials,1);zeros(n_trials,1)];
%             w_opt = pinv(X'*X)*X'*Y;
%             y_est = X*w_opt;
%             error(n,m) = norm(y_est-Y)^2;
%             a = corrcoef(y_est,Y);
%             R(n,m) = a(2,1)^2;
%         end
%     end
% end
% ------------------------ no PCA------------------------------------------
if nargin < 3, if_plot = 0; end

% n_pc = 6; % number of principle component
A = A';
% A = A-mean(A);
N = round(size(A,2)/n_trials);
R = zeros(N,N);
error = ones(N,N);
for n = 1:N
    for m = 1:N
        if n==m
            continue
        else
            X = [A(:,(n-1)*n_trials+1:n*n_trials), A(:,(m-1)*n_trials+1:m*n_trials)];
            if isempty(find(X~=0))
                continue
            end
            X = X';
            [U,S,V] = svd(X);
            lambda = diag(S(1:20,1:20));
            [~,n_pc] = min(abs(lambda-lambda(1)/2)); % number of principle component
            X = X*V(:,[1:n_pc]);
            X = X-mean(X);
            Y = [ones(n_trials,1);-1*ones(n_trials,1)];
            Y = Y./norm(Y);
            w_opt = pinv(X'*X)*X'*Y;
            y_est = X*w_opt;
            error(n,m) = norm(y_est-Y)^2;
%             a = corrcoef(y_est,Y);
%             R(n,m) = a(2,1)^2;
            R(n,m) = (1-norm(y_est-Y)^2)/(norm(Y-mean(Y))^2);
        end
    end
end
for n = 1:N
    for m = 1:N
        if isnan(R(n,m))&&(~isnan(R(m,n)))
            R(n,m) = R(m,n);
        end
        if R(n,m)>R(m,n)
            R(m,n) = R(n,m);
            error(m,n) = error(n,m);
        end
    end
end
if if_plot
    figure
    colors = [[0, 145, 69];[247,147,30]]/255;
    scatter(X(1:n_trials,1),X(1:n_trials,2), 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'none')
    hold on
    scatter(X(n_trials+1:n_trials*2,1),X(n_trials+1:n_trials*2,2), 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'none')
    hold on
    d = 0.01;
    [x1Grid,x2Grid] = meshgrid(min(X(:,1))-0.02:d:max(X(:,1))+0.02,...
    min(X(:,2))-0.02:d:max(X(:,2))+0.02);
    xGrid = [x1Grid(:),x2Grid(:)];
    Mdl1 = fitcsvm(X(:,1:2),Y,'KernelFunction','linear','Standardize',true);
    hold on
    [~,scores1] = predict(Mdl1,xGrid);
    contour(x1Grid,x2Grid,reshape(scores1(:,2),size(x1Grid)),[0 0],'k');
    box on
end
end