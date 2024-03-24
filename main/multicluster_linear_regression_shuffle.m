function [R, error] = multicluster_linear_regression_shuffle(A, n_trials)
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
            X = X(:, randperm(size(X, 2)));
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
end