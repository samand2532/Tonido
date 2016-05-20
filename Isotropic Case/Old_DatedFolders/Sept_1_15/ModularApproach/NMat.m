clc; clear all;

N=6;
AllBras = csvread('AllBrasLLL.csv');
AllKets = AllBras';

NMatFinal = (eye(length(AllBras)).*N)