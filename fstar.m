function [fstar] = fstar(f, f_eq, t)
%% Function [fstar] = fstar(f, f_eq, t) computes the post-collision PDF based on the given PDF, equalibrium PDF, and relaxation time
%% f is the PDF and must be a row vector
%% f_eq is the equalibrium PDF and must be a row vector
%% t is the time relaxation constant and must be a scalar greater than 0.5
%% fstar is the post-collision PDF and will be a row vector

% Calculate fstar
fstar = f - (1/t) * (f - f_eq);