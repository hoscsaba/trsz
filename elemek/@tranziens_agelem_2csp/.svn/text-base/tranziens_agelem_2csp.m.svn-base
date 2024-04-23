function trag2 = tranziens_agelem_2csp(varargin)

trag2.nev = varargin{1};
trag2.csp = [varargin{2} varargin{3}];
trag2.Q   = varargin{4};
trag2.Qr  = varargin{4};
trag2.p   = [1e5 1e5];

trag2.folynev = 'viz';
trag2.ro      = 1e3;
trag2.nu      = 1e-6;
trag2.mu      = 1e-3;
trag2.B       = 1e9;
trag2.resfile = 0;

%% 2009. 08. 24. Hos Csaba
% A fignum-ot nem szamoljuk, hanem veletlengeneratorral gyartjuk
trag2.fignum = round(rand(1)*1000);

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
trag2 = class(trag2,'tranziens_agelem_2csp');
