clear; close all; clc;

[~,~,a]=laplacian([10 10],{'NN','NN'});
a = full(a);

a_cpp = load('Test.txt');

isequal(a,a_cpp)
