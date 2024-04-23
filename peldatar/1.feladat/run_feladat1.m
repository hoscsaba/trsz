clear all, clc, clf

%% konyvtarak beallitasa
addpath('/Users/hoscsaba/program/trsz')
addpath('/Users/hoscsaba/program/trsz/lib')
addpath('/Users/hoscsaba/program/trsz/elemek')

%% Szamitas futtatasa, tmax=30s
trsz('feladat1.tpr',30)

%% Eredmeny rajzolasa
trsz_rajz('feladat1.tpr','cso','ev','p')