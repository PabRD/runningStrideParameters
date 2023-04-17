%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                      NAME: Traitement Accelerometre                     %
%                      AUTHOR: PabDawan                                   %
%                      DATE: Février/Mai 2022                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Description: Ouverture des fichiers bruts et extraction des cadences,
%temps de contact, temps de vol, raideur verticale, raideur de jambe
tic
clc                                                                         %Clear Command Window
clear                                                                       %Clear Workspace
close all                                                                   %Close figure

load rawdata

%% Structuration des données

acc_1=rawdata(1).physilogData(1).data(:,1:3); %Slave
acc_2=rawdata(2).physilogData(1).data(:,1:3); %Master

gyro_1=rawdata(1).physilogData(2).data(:,1:3);
gyro_2=rawdata(2).physilogData(2).data(:,1:3);

%TEMPS
t_G=rawdata(1).physilogData(1).timestamps;   
t_D=rawdata(2).physilogData(1).timestamps;   

% plot(t_G,acc_1)
% hold on
% plot(t_D,acc_2)
% legend(['x'; 'y'; 'z'])


acc_result_G=sqrt(acc_1(:,1).^2+acc_1(:,2).^2+acc_1(:,3).^2);     %Gauche Slave                       % Je calcule l'accélération résultante
acc_result_D=sqrt(acc_2(:,1).^2+acc_2(:,2).^2+acc_2(:,3).^2);     %Droite Master                       % Je calcule l'accélération résultante
gyro_result_G=sqrt(gyro_1(:,1).^2+gyro_1(:,2).^2+gyro_1(:,3).^2);
gyro_result_D=sqrt(gyro_2(:,1).^2+gyro_2(:,2).^2+gyro_2(:,3).^2);

% plot(t_G,acc_result_G*100)
% hold on
% plot(t_D,acc_result_D*100)
% plot(t_D,gyro_result_D)
% plot(t_G,gyro_result_G)

%Préparation filtre
Sig{1}=acc_result_G;
Sig{2}=acc_result_D;
Sig{3}=gyro_result_G;
Sig{4}=gyro_result_D;
%plot(t,[gxbrut,gybrut,gzbrut])
%plot(t,gxbrut)
%hold on
%plot(t,aybrut)
%legend('gyrox','accy','peak start')
%% Découper le signal pour ne traiter que 20 secondes de course à pied/ extraction de sous matrice

mark1=75;                                                                   % marqueur auquel je veux commencer mon traitement
mark2=105;                                                                  % marqueur auquel je veux finir mon traitement
index1=find(abs(t_G-mark1)==min(abs(t_G-mark1)));                           % index1 sert a stocker la position de la valeur la plus proche de 80sec
index2=find(abs(t_G-mark2)==min(abs(t_G-mark2)));                           % pareil: je ne suis pas sur de trouver la valeur exacte 80 ou 100 donc je cherche la plus proche
t_G=t_G(index1:index2,:);                                                   % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
Signal_coupe{1}=Sig{1}((index1:index2),1);                                  % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
Signal_coupe{3}=Sig{3}((index1:index2),1);                                  % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet

index1=find(abs(t_D-mark1)==min(abs(t_D-mark1)));                           % index1 sert a stocker la position de la valeur la plus proche de 80sec
index2=find(abs(t_D-mark2)==min(abs(t_D-mark2)));                           % pareil: je ne suis pas sur de trouver la valeur exacte 80 ou 100 donc je cherche la plus proche
t_D=t_D(index1:index2,:); 
Signal_coupe{2}=Sig{2}((index1:index2),1);                                  % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet
Signal_coupe{4}=Sig{4}((index1:index2),1);                                  % apres avoir pris les indexes, on crée une sous matrice avec les 20 secondes d'interet



%% Lowpass sur GYRO car ce que je veux observer
sample_rate=rawdata(1).physilogData(1).Fs;                                  % Fréquence d'acquisition  
S = Signal_coupe;                                                           % Signal à traiter

fc=50;                                                                      % fréquence de coupure
fs=sample_rate;
[b_l,a_l] = butter(3,fc/(fs/2),'low');
S_filt_ex{1}=filtfilt(b_l,a_l,S{3});                                        % Gauche
S_filt_ex{2}=filtfilt(b_l,a_l,S{4});                                        % Droite

%% Application du Butterworth Lowpass 50 Hz a phase nulle sur l'accelero

fc=50;                                                                      %fréquence de coupure
fs=sample_rate;
[b_l,a_l] = butter(3,fc/(fs/2),'low');
S_filt_lp = cellfun(@(x) filtfilt(b_l,a_l,x),S,'uni',0);

%% Détéction de pics qui va me servir à trouver les pics de fin de contact sur les vitesses angulaires
% Sur cette partie on peut jouer sur 'Prominence', 'Width', 'Height' si le
% traitement ne fonctionne pas afin de detecter les bons pics
% col=cbrewer('div','BrBG',11); %palette de couleur graphe (attention toolbox cbrewer nécéssaire)
col = [    0.3294    0.1882    0.0196
    0.5490    0.3176    0.0392
    0.7490    0.5059    0.1765
    0.8745    0.7608    0.4902
    0.9647    0.9098    0.7647
    0.9608    0.9608    0.9608
    0.7804    0.9176    0.8980
    0.5020    0.8039    0.7569
    0.2078    0.5922    0.5608
    0.0039    0.4000    0.3686
         0    0.2353    0.1882];

[y_fin_G,locs_fin_G]=findpeaks(-S_filt_ex{1,1},'MinPeakProminence',5,'MinPeakWidth',5,'MinPeakHeight',-600,'MinPeakDistance',10);
[y_fin_D,locs_fin_D]=findpeaks(-S_filt_ex{1,2},'MinPeakProminence',1,'MinPeakWidth',1,'MinPeakHeight',-600,'MinPeakDistance',10);


%% Trouver les pics de fin de contact parmi tous les pics détéctés
% Pour cela; les pics de fins de contacts arrivent après une periode de
% très faible vitesse angulaire: j'isole ce pic puis je séléctionne le pic d'après

%% Pied Gauche
longPics=length(locs_fin_G);                                                % Nombre de pics détéctés
picsClassmt=sort(-S_filt_ex{1}(locs_fin_G));                                % Classer la hauteur des pics dans l'ordre croissant
a=findchangepts(picsClassmt,"MaxNumChanges",2);                             %trouver le/les points de changement radical
%nombre de pics interressants - marge de sécu
nb_pics=longPics-max(a);                                                    % Le deuxième gros changement correspond au passage de pics intermédiaires aux pics à vitesse très basse
low_bnd=mean(picsClassmt(max(a):end));                                      % Fenetre basse fixée au point le plus petit des ces pics à vitesse faible
high_bnd=max(picsClassmt(max(a):end));                                      % Fenetre haute fixée au point le plus grand des ces pics à vitesse faible

% Detection des pics à plus faible vitesse angulaire (phase ou le pied est posé au sol: la phase qui suit (propulsion) sera le pic qui nous interresse

%pied gauche
locsPrePics=find(locs_fin_G(y_fin_G>=low_bnd));                             % Parmi tous les pics détéctés, trouve moi ceux qui sont supérieur à ma fenetre basse

picVzero=strfind([y_fin_G>=low_bnd]',[1 0]);                                %passage d'un pic a vitesse très faible à un pic qui n'est pas dans la fenetre (passage de 1 à 0 en booléen)
locs_locs=picVzero+1;                                                       %Séléction de tous les pics suivant

%locs_locs=find(y_fin_G>=low_bnd)+1;
locs_locs(end-5:end)=[];                                                    % j'enleve 5 pics de fin pour pas me retoruver avec une valeur qui excede la valeur possible car index+1
locs_fin_G=locs_fin_G(locs_locs);

%% Pied Droit
longPics=length(locs_fin_D);
picsClassmt=sort(y_fin_D);
a=findchangepts(picsClassmt,"MaxNumChanges",2);
%nombre de pics interressants - marge de sécu
nb_pics=longPics-max(a);
low_bnd=mean(picsClassmt(longPics-nb_pics:end));
high_bnd=max(picsClassmt(longPics-nb_pics:end));

%Pied droit
locsPrePics=find(locs_fin_D(y_fin_D>=low_bnd));

picVzero=strfind([y_fin_D>=low_bnd]',[1 0]); %reflexion solution car moins on filtre, plus on voit deux gros pics a basse velocité
locs_locs=picVzero+1;

%locs_locs=find(y_fin_D>=low_bnd)+1;
locs_locs(end-5:end)=[]; 
locs_fin_D=locs_fin_D(locs_locs);

%% Lorsque je peux pas detecter je reprends la main:
% plot(S_filt_ex{1})
% [loc_manu,~]=ginput;
% loc_manu=floor(loc_manu);
% loc_manu_classe=sort(loc_manu);

%% detection pic d'impact
S_filt=S_filt_lp;
maxPeakHeight=max(S_filt{1});
[~,pics_acc]=cellfun(@(x) findpeaks(x,'MinPeakHeight',0.4*maxPeakHeight,'MinPeakProminence',5,'MinPeakDistance',fs*200/1000*2),S_filt,'UniformOutput',false);

%[~,pics_acc{1}]=findpeaks(S_filt{2},'MinPeakHeight',0.4*maxPeakHeight,'MinPeakProminence',5,'MinPeakDistance',fs*200/1000*2)

FPA_mean(1)=mean(S_filt{1}(pics_acc{1}));                                   %mean Foot Peak Acceleration
FPA_mean(2)=mean(S_filt{2}(pics_acc{2}));                                   %mean Foot Peak Acceleration

FPA_sd(1)=std(S_filt{1}(pics_acc{1}));                                      %standrad deviation
FPA_sd(2)=std(S_filt{2}(pics_acc{2}));                                      %standrad deviation

%% supprimer les points en trop entre deux pics de decceleration

NewLocsFinG=nan(1,length(locs_fin_G));                                      % Pré allocation
NewLocsFinG(1)=locs_fin_G(1);                                               % Pré allocation

for kk=1:length(pics_acc{1})-1
    xdeb=pics_acc{1}(kk);                                                   %Gauche
    xfin=pics_acc{1}(kk+1);
    xInRange = (locs_fin_G>=xdeb) & (locs_fin_G<=xfin);                     % ya t'il plus que un 1 entre deux pics d'accel?
    RangeNumbPic=sum(xInRange==1);                                          %combien de locs fin entre mes pics d'accel?
%     disp(RangeNumbPic)
    premsUN=find(xInRange,1,"first");                                       % Selectionner uniquement le premier pic
    NewLocsFinG(premsUN)=locs_fin_G(premsUN);
end
NewLocsFinG(kk+2:end)=[];
NewLocsFinG=rmmissing(NewLocsFinG);

NewLocsFinD=nan(1,length(locs_fin_D));
NewLocsFinD(1)=locs_fin_D(1);

for kk=1:length(pics_acc{2})-1
    xdeb=pics_acc{2}(kk);                                                   %Droite
    xfin=pics_acc{2}(kk+1);

xInRange = (locs_fin_D>=xdeb) & (locs_fin_D<=xfin);                         % ya t'il plus que un 1 entre deux pics d'accel?
RangeNumbPic=sum(xInRange==1);                                              %combien de locs fin entre mes pics d'accel?
% disp(RangeNumbPic)
premsUN=find(xInRange,1,"first");
NewLocsFinD(premsUN)=locs_fin_D(premsUN);
end
NewLocsFinD(kk+2:end)=[];
NewLocsFinD=rmmissing(NewLocsFinD);
%% Visualisation: cette partie permet de controler visuellement si il y a un problème dans le traitement pour ensuite jouer sur certains facteurs ci dessus
% Visualisation de detection des pics de fin de contact pour vérification
plot(t_G,S_filt_ex{1},'Color',col(10,:))
hold on
scatter(t_G(locs_fin_G),S_filt_ex{1}(locs_fin_G),'filled','CData',[1 0 0],'HandleVisibility','off')
scatter(t_G(NewLocsFinG),S_filt_ex{1}(NewLocsFinG),'filled','CData',col(9,:))
plot(t_G,S_filt_lp{1}*100,'color',col(8,:))

% Visualisation (suite)
plot(t_D,S_filt_ex{2},'Color',col(4,:))
scatter(t_D(locs_fin_D),S_filt_ex{2}(locs_fin_D),'filled','CData',[1 0 0],'HandleVisibility','off')
scatter(t_D(NewLocsFinD),S_filt_ex{2}(NewLocsFinD),'filled','CData',col(4,:))
plot(t_D,S_filt_lp{2}*100,'color',col(5,:))

legend(["gyro G","finContact G","acc G","gyro D","finContact D","acc D"])
ax=gca;
ax.XLim=[75 85];
%Tracer les pics d'accélérations
scatter(t_G(pics_acc{1}),S_filt{1}(pics_acc{1})*100,'filled','CData',col(9,:),'HandleVisibility','off')
scatter(t_D(pics_acc{2}),S_filt{2}(pics_acc{2})*100,'filled','CData',col(4,:),'HandleVisibility','off')
hold off

%% calcul temps de contact et temps de vol sur 15 pas
newpics_acc=cell(1,2);

newpics_acc{1}=[pics_acc{1};nan(50-length(pics_acc{1}),1)];                 % je remplis tout avec des NaN pour avoir des vecteurs de même taille
newpics_acc{2}=[pics_acc{2};nan(50-length(pics_acc{2}),1)];                 % je remplis tout avec des NaN
newpics_acc{3}=[pics_acc{3};nan(50-length(pics_acc{3}),1)];                 % je remplis tout avec des NaN
newpics_acc{4}=[pics_acc{4};nan(50-length(pics_acc{4}),1)];                 % je remplis tout avec des NaN

locs_fin_D=NewLocsFinD';
locs_fin_G=NewLocsFinG';

%% Le calcul des temps de contacts se fait via les deux signaux: l'ordre dans lequel les pics détéctés arrivent dans la fenetre est donc important
% Ici je dois trouver les cas ou le premier pic est un pic d'impact ou un pic de fin de contact
locs_fin_D=[locs_fin_D;nan(50-length(locs_fin_D),1)];
locs_fin_G=[locs_fin_G;nan(50-length(locs_fin_G),1)];

% PIED DROIT
if locs_fin_D(1)<pics_acc{2}(1)                                             % Si le pic de fin arrive avant le pic de début
    newlocs_fin_D=locs_fin_D(2:16);                                         %Alors je prends à partir du deuxième pic de fin, et je prends 15 pas
    newpics_acc{2}=pics_acc{2}(1:16-1);                                     %Alors je prends les 15 premiers pic de début
else                                                                        %Si le premier pic de fin arrive bien après le premier pic de début
    newlocs_fin_D=locs_fin_D(1:15);                                         %Alors je prends bien les 15 premiers pics de fin
    newpics_acc{2}=pics_acc{2}(1:15);                                       %Alors je prends bien les 15 premiers pics de début aussi
end
%PIED GAUCHE
if locs_fin_G(1)<pics_acc{1}(1) %G
    newlocs_fin_G=locs_fin_G(2:16);
    newpics_acc{1}=pics_acc{1}(1:16-1);
else
    newlocs_fin_G=locs_fin_G(1:15);
    newpics_acc{1}=pics_acc{1}(1:15); %G
end



%% Calcul des paramètres biomecaniques initiaux

% Cadence (ppm)
t_pics=[t_G(pics_acc{1}(1:15)) t_D(pics_acc{2}(1:15))];
cad=60./abs(diff(t_pics,1,2));                                              %diff de colonnes absolue
cad_moyenne=mean(cad);                                                     	% La cadence c'est 60/(tContact+tVol) donc bien le temps d'un cycle complet , c'est à dire d'un contact gauche à un contact droite
cad_SD=std(cad);

%temps de contact a droit et gauche (s)
tc_G=diff([t_G(newpics_acc{1}) t_G(newlocs_fin_G)],1,2);                    %G
tc_D=diff([t_D(newpics_acc{2}) t_D(newlocs_fin_D)],1,2);                    %D

t_contact_moyen_G=mean(tc_G);
t_contact_moyen_D=mean(tc_D);

%% Temps de vol (s)
% Attention: cette fois le temps de vol est calculé à l'aide des deux signaux gauche et droite
% Je commence par le premier pic d'acceleration qui arrive dans ma fenetre
% puis en fonction, je prends le premier point de fin contact avec le premier point de debut de contact de l'autre pied, etc
% Pour l'autre signal (celui qui arrive en deuxième), je prends le premier point de fin de contact mais le deuxième point de debut de contact de l'autre pied 

%booléens
bool=newpics_acc{1}(1)<newpics_acc{2}(1); %est ce que le pic gauche de départ arrive avant le pic droit de départ

switch bool
    case 1 %Si premier pic de la fenetre est un pic du pied gauche
        t_vol_G=abs(t_G(newlocs_fin_G)-t_G(pics_acc{2}(1:length(newlocs_fin_G))));    %temps de vol apres impact gauche
        t_vol_D=abs(t_D(newlocs_fin_D)-t_D(pics_acc{1}(2:length(newlocs_fin_D)+1)));     %temps de vol apres impact droit

    case 0 %Si premier pic de la fenetre est un pic du pied droit
        t_vol_G=abs(t_G(newlocs_fin_G)-t_G(pics_acc{2}(2:length(newlocs_fin_G)+1)));    %temps de vol apres impact gauche
        t_vol_D=abs(t_D(newlocs_fin_D)-t_D(pics_acc{1}(1:length(newlocs_fin_D))));      %temps de vol apres impact droit
end

t_vol_moyen_G=mean(t_vol_G);
t_vol_moyen_D=mean(t_vol_D);

t_vol_SD_G=std(t_vol_G);
t_vol_SD_D=std(t_vol_D);

%% Rigidité de la jambe (kN/m) & Rigidité verticale (kN/m) sur la base de Morin et al. (2005)
t_contact_G=tc_G;
t_contact_D=tc_D;

%Parametres d'entrée
g=9.81;                         %acceleration terrestre (m/s²)
m=50;                           %masse (kg)
h=1.71;                         %hauteur du sujet (m)
tf=[t_vol_G t_vol_D];           %temps de vol (s)
tc=[t_contact_G t_contact_D];   %temps de contact (s)
L=0.53*h;                       %Longueur de jambe (m)                      %according to the anthropometric equations of Winter (1979)
v=10;                           %vitesse de course (kmh)


%Force de réaction au sol maximale (modélisée)
%1 Newton= 1 kg m s−2
Fmax_G=(m*g*(pi/2)*((tf(:,1)./tc(:,1))+1))/1000; % kN
Fmax_D=(m*g*(pi/2)*((tf(:,2)./tc(:,2))+1))/1000; % kN

%Deplacement vertical du CM (modélisé)
delta_yc_G=(((Fmax_G.*tc(:,1).^2)/(m*pi^2))+g.*((tc(:,1).^2)/8)); %m
delta_yc_D=(((Fmax_D.*tc(:,2).^2)/(m*pi^2))+g.*((tc(:,2).^2)/8)); %m

%Rigidité verticale (modélisée)
k_vert_G=Fmax_G./delta_yc_G; %kN.m-1
k_vert_moyen_G=mean(k_vert_G);
k_vert_SD_G=std(k_vert_G);

k_vert_D=Fmax_D./delta_yc_D; %kN.m-1
k_vert_moyen_D=mean(k_vert_D);
k_vert_SD_D=std(k_vert_D);


%rigidité de la jambe (modélisée)
v_ms=v/3.6; %m/s

delta_L_G=L-sqrt((L^2)-((v_ms.*tc(:,1))/2).^2)+delta_yc_G;
delta_L_D=L-sqrt((L^2)-((v_ms.*tc(:,2))/2).^2)+delta_yc_D;

k_leg_G=Fmax_G./delta_L_G;
k_leg_D=Fmax_D./delta_L_D;

k_leg_moyen_G=mean(k_leg_G);
k_leg_SD_G=std(k_leg_G);

k_leg_moyen_D=mean(k_leg_D);
k_leg_SD_D=std(k_leg_D);

%% extraction fichier traité
% On extrait un tableau avec l'ensemble des paramètres, on retire les
% outliers (causé par un fail de la part du sujet, ou un fail de
% traitement)
rapport_final=[cad,t_contact_G,t_contact_D,t_vol_G,t_vol_D,k_vert_G,k_vert_D,k_leg_G,k_leg_D,S_filt{1}(newpics_acc{1}(1:length(t_contact_G))),S_filt{2}(newpics_acc{2}(1:length(t_contact_D)))];
% outliers=isoutlier(rapport_final,"quartiles");                            %1,5 fois le range interquartile est considéré comme outlier
outliers=isoutlier(rapport_final);                                          %MAD = median(∣Ai−median(A)∣) median absolute derivation
rapport_final(outliers)=nan;
tab=array2table(rapport_final);                                             %Transformation d'une matrice à un tableau
% tab=table(cad,t_contact_G,t_contact_D,t_vol_G,t_vol_D,k_vert_G,k_vert_D,k_leg_G,k_leg_D,S_filt{1}(newpics_acc{1}(1:length(t_contact_G))),S_filt{2}(newpics_acc{2}(1:length(t_contact_D))));
tab.Properties.VariableNames={'Cadence (ppm)','temps contact gauche (s)','temps contact droite (s)','temps vol gauche (s)','temps vol droit (s)','raideur verticale gauche (kN/m)','raideur verticale droite (kN/m)','raideur de jambe gauche (kN/m)','raideur de jambe droite (kN/m)','FPA gauche (g)','FPA droite (g)'}



%% A faire:
% Calcul de l'angle d'attaque pour determiner si ataque talon, midfoot , forefoot

toc
