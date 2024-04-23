function trag1 = tranziens_agelem_1csp(varargin)

trag1.nev = varargin{1};
trag1.csp = varargin{2};
trag1.Q   = varargin{3};
trag1.Qr  = varargin{3};
trag1.p   = 1e5;

trag1.folynev = 'viz';
trag1.ro      = 1e3;
trag1.nu      = 1e-6;
trag1.mu      = 1e-3;
trag1.B       = 1e9;
trag1.resfile = 0;

%% 2009. 08. 24. Hos Csaba
% A fignum-ot nem szamoljuk, hanem veletlengeneratorral gyartjuk
trag1.fignum = round(rand(1)*1000);

%% Fignum kiszamitasa
% chartable=char(32:127);
% charname=char(varargin{1});
% fignumero=1;
% for i=1:length(chartable)
%     for j=1:length(charname)
%         if strcmp(chartable(i),charname(j))
%             fignumero=fignumero+10*abs(round(10*sin(i)))+abs(round(10*cos(j)));
%         end
%     end
% end
% trag2.fignum=fignumero;

%% Letrehozas
trag1 = class(trag1,'tranziens_agelem_1csp');
