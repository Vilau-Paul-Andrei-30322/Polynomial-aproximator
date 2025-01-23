%% Project Part 1 : Fitting a unknown function in a polynomial approximation
    
% This script approximates an unknown function using a polynomial regression model
% It calculates the mean squared error (MSE) for various polynomial degrees
% and determines the optimal polynomial degree for fitting

    clc; clear; close all;
    load ('proj_fit_12.mat');
    
    m = input('Enter the degree of your polynomial: ');
    
    % Declaring every variable needed for better performance
    % e.g. preallocating errors to avoid resizing
    X_id = id.X;
    Y_id = id.Y;
    X_val = val.X;
    Y_val = val.Y;
    n = length(X_id{1});
    N = length(X_val{1});
    Y_flat = Y_id(:);
    Y_val_flat = Y_val(:);
    errors = zeros(1,m);
    errors_val = zeros(1,m);
    x1 = X_val{1};
    x2 = X_val{2};
    
    % Creating Regressor Matrix
    R = Polynomial(X_id{1}, X_id{2}, m, n);
    R_val = Polynomial(X_val{1}, X_val{2}, m, N);
    
    % Computing the approximation model, coefficients of the
    % polynomial and mean squared errors (MSE) for each degree
    for degree = 1:m
        % these R's are created to ensure the traversing of 
        % the matrix on each row, also the formula (d+1)*(d+2)/2 
        % calculates the number of terms in a polynomial
        R_degree = R(:, 1:(degree+1)*(degree+2)/2);
        Rval_degree = R_val(:, 1:(degree+1)*(degree+2)/2);
    
        theta = R_degree\Y_flat;
        Y_hat = R_degree*theta;
        Y_hat_val = Rval_degree*theta;
    
        errors(degree) = mean((Y_flat - Y_hat).^2);
        errors_val(degree) = mean((Y_val_flat - Y_hat_val).^2);
    end

    % Taking the optimal error and degree for which the function fits
    [min_error, optimum_degree] = min(errors_val);
    disp(min_error);
    disp(optimum_degree);
    % Plotting the evolution of the error on 
    figure
    plot(1:m,errors_val, 'g'); hold on;
    plot(optimum_degree, min_error, 'mh', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
    xlabel('Polynomial Degree');
    ylabel('Mean Squared Error');
    
    %%
    % Using only the optimal values
    R_optimal = R(:, 1:(optimum_degree+1)*(optimum_degree+2)/2);
    Rval_optimal = R_val(:, 1:(optimum_degree+1)*(optimum_degree+2)/2);
    theta_optimal = R_optimal \ Y_flat;
    Y_hat_optimal = R_optimal * theta_optimal;
    Y_hat_val_optimal = Rval_optimal * theta_optimal;
    Y_hat_val_optimal_reshaped = reshape(Y_hat_val_optimal, size(Y_val));
    
    % Comparing the true validation data and the predicted values using the optimal polynomial
    figure
    mesh(X_val{1}, X_val{2}, Y_val, 'FaceColor', 'c'); hold on;
    colormap(cool)
    mesh(X_val{1}, X_val{2}, Y_hat_val_optimal_reshaped,'FaceColor', 'interp');
    colormap(cool)
    title('Validation Data: True vs Predicted');
    legend('True', 'Predicted');


%%
% This function creates the matrix for polynomial regression
% Each row corresponds to a point (x1, x2) in the grid, and the columns
% are the monomials of the polynomial up to the chosen degree
function matrix = Polynomial(x1, x2, m, n)
    matrix = [];
    for i = 1:n
        for j = 1:n
            polynomial_terms = 1;
            for degree = 1:m
                for k = 0:degree
                    l = degree - k;
                    monomial = x2(i).^k * x1(j).^l;
                    polynomial_terms = [polynomial_terms, monomial];
                end
            end
            matrix = [matrix; polynomial_terms];
        end
    end
end
