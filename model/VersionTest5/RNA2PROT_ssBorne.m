% Modèle simple interaction RNA-Protéines
%
% Lecture des données: nom, (DPA_i,r_i,p_i) i=1..n
% Normalisation des données par rapport à la moyenne
% P_i=p_i/p0 et R_i=r_i/r0 où p0=mean(p_i) et r0=mean(r_i)
% lissage des données RNA par régression polynomiale ordre 3
% R(t)=polynome de degré 3
% Equation différentielle pour modéliser P
% dP/dt=Ks*R(t)-kd*P
% d'où P(t)=poly deg 3 + a*exp(-kd*t)
% calibrage des paramètres Ks et kd: moindres carrés
% S(Ks,kd)=sum(P(DPA_i)-P_i)^2
% sensibilité de P par rapport à (Ks,kd)

%provisoire
global ptrans3 ptrans6 St6 mt6 St3 mt3 plogtrans3 Slogt3 mlogt3 degre h fitR
global xi yi
global DPA_R DPA_P R
global fitPoids parv parg parc parl pare pard Wp3 Sp3 mp3

% repertoire courant;
chemin=pwd();
% ajout au path du répertoire stixbox
addpath([chemin,'\stixbox']);

% Lecture des données dans fichier(s)
% RNA table de structures et PROT table de structures contenant
% NOM cell
% DPA matrice colonne
% data matrice colonne

% Options à modifier à la main 
% calcul des dommaines de confiance oui= 1 non= 0
% calcul de la matrice de sensibilité/parametres
stat=1;
% sauvegarde des résultats et graphiques oui= 1 non= 0
sauve=1;


close all;
repertoire_a_modeliser ='donnees_a_modeliser/';
%repertoire_a_modeliser ='donnees_a_modeliser/';
% nom des fichiers par défaut
%nomfich1 = '50PROT';
%nomfich2 = '50RNA';
%nomfich3= 'POIDS';
%nomfich1 = '538P1_model';
%nomfich2 = '538RNA';
%nomfich3= 'POIDS';
nomfich1 = '37P';
nomfich2 = '37R';
nomfich3= 'POIDS';
nomfich3='PoidsCerise';
% Saisie nom répertoire et fichier des données (0=existence-mode lecture)
[File1] = LectureNomFichier(repertoire_a_modeliser,[nomfich1 '.mat'],0);
[File2] = LectureNomFichier(repertoire_a_modeliser,[nomfich2 '.mat'],0);
[File3] = LectureNomFichier(repertoire_a_modeliser,[nomfich3 '.mat'],0);

% lecture des fichiers: PROT,RNA et POIDS
load(File1);load(File2);load(File3);

% ajout de données manuelles: (0,0)
%[PROT,RNA,poids]=ajoutdata(PROT,RNA,poids);
xi=poids(:,1);yi=poids(:,2);

% choix de la courbe de croissance (fitPoids variable globale)
[errPoids,messagePoids]=fit_testpoids(xi,yi);
cpt_fen=4;%compteur fenetre

stop=0
if stop==1
% sauvegarde figures poids et fermeture fenetres
if sauve==1
    for i=1:cpt_fen
        nomfich=['graphes/Poids'  num2str(i),'.pdf']
        saveas(figure(i),nomfich)
        % courbe pour publi
        saveas(figure(5000),'graphes/PoidsCourbe.svg')
    end
    close all %fermeture fenetre
    cpt_fen=0;
end

% choix fit de R
[ fitR,messageR] = choix_fitRNA;

% nombre de couples (PROT/RNA) des fichiers
nb_indiceP=numel(PROT);
% choix des couples à modéliser
[ indiceP ] = choix_data(nb_indiceP)

% ks fixé: Ks=6*0.486*3*3600*24/Lp (jour-1)
% valeur minimum de ks
ksmin=3*4*3*3.6*24;%ksmax=6*10*3*3.6*24=5*ksmin

% pour chaque couple (RNA,proteine) d'indice contenu dans indiceP
cpt_boucle_fin=numel(indiceP);
for cpt_boucle=1:cpt_boucle_fin
    
    i=indiceP(cpt_boucle)
    cpt_fen=cpt_fen+1;
    % compteur fenetre couple (rna, prot) en cours
    cpt_fen_debut=cpt_fen;
    % temps (DPA) des mesures en gFW
    DPA_P=PROT(i).DPA;RNAg=RNA(i).data;
    DPA_R=RNA(i).DPA;PROTg=PROT(i).data;
    % tri par ordre croissant des DPA
    [DPA_P,indiceDPA]=sort(DPA_P);PROTg=PROTg(indiceDPA);
    [DPA_R,indiceDPA]=sort(DPA_R);RNAg=RNAg(indiceDPA);
    % nom solyc prot/rna
    NOM=RNA(i).NOM;
    % score optimisation  initialisation
    score=0;
    
    % CALCUL SUR DONNEES/gFW
    % normalisation données brutes et du ks correspondant
    % titre=[NOM '(/gFW)'];
    titre=[NOM ];
    [Rg,Pg,ks_nmin,R0g,P0g]=normalisation(cpt_fen,titre,DPA_R,RNAg,DPA_P,PROTg,ksmin);
    kdmin=ks_nmin
    cpt_fen=cpt_fen+1;
    % choix de la fonction t->RNA(t)(données brutes en gFW)
    R=Rg;
    errR=fit_testRNA(cpt_fen,titre,DPA_R,R);
    % score du fit RNA
    %score=ScoreErreur( errR,score );
    
    % détermination du couple(kd,ks) à partir des données brutes normalisées(R,P)
    % en utilisant la fonction t->mu(t) taux de croissance
    % dp/dt=ks*r-(kd+mu)*p (données en gFW)
    [kparg,fprotg,errg,scoreOF,messageOF]=solg_ssBorne(NOM,DPA_P,Pg,ks_nmin,score);
    %kparg
    %errg
    % affichage des résultats /gFW
    cpt_fen=cpt_fen+1;
    fit_Prot(cpt_fen,titre,DPA_P,Pg,kparg,DPA_P,fprotg,errg,R0g,P0g,ks_nmin,scoreOF);
    
    if stat==1
        cpt_fen=cpt_fen+1;figure(cpt_fen);
        [X,Y,Z]=matrice_sensibilite(DPA_P,kparg(:,1));
        rho=matrice_correlation(X,numel(Pg),(errg(1)*norm(Pg,2))^2)
        hx=plot(unique(DPA_P),[X(:,1),X(:,2),-X(:,3)]);hold on;
        title(['sensibilité ' NOM],'FontWeight','bold')
        xlabel('DPA ')
        ylabel('sensibilité/paramètres (kmin)')
        legend('dp/dp0','dp/dks','dp/dkd');
        limaxis=axis;
        xPos=limaxis(1);
        yl = limaxis(4);
        yPos = yl*.7;
        texte=text(xPos,yPos,sprintf('correlation=%8.4f ',rho(2,3)),'FontSize',8);
        % si la correlation est proche de 1
        % tentative de reparametrisation
        options = optimset('Algorithm','sqp','TolX',1e-5,'TolFun', 1e-5,'Tolcon',1e-5,'MaxFunEvals',1000);
        [P0,Pmin,Pmax]=initiale(DPA_P,Pg);lb=[Pmin;0;0];
        ub=[Pmax;inf;inf];
        % fminunc : min sans contrainte
        para=fminunc(@(par)fonction_ecart(par,DPA_P,Pg),kparg(:,1))
        if find(para<0)% si l'un des paramètres est négatif
            % optimisation avec contraintes de positivité
            [para,fval2,exitflag2,output2] = fmincon(@(par)fonction_ecart(par,DPA_P,Pg),kparg(:,1),[],[],[],[],lb,ub);
        end
        
        % courbe de niveau de la fonction ecartMC(ks,kd) au voisinage du min
        cpt_fen=cpt_fen+1;
        [c1,seuil2,Aire,scoreS,messageS]=courbe_niveau_ellipse(titre,X,kparg,errg,Pg,cpt_fen);
        
    end
    
    %%%%%%%%%%%%%%%% SAUVEGARDE FICHIERS
    if sauve==1
        % sauve(filename,i,NOM,ks_n,kd(i-1),par,DPA,P,sol);
        nom_Fichier_mat=['resultats/' NOM '.mat'];
        %fid = fopen(mon_Fichier, 'a');%ouverture: écriture à la fin du fichier
        % fprintf(fid, '%5.2f %10.4f\n', [x ; exp(x)]);
        %save(nom_Fichier_mat,'kparg','P0g','R0g','errg','fitR','errR','fitPoids','errPoids','score') ;
        
        % définition de la structure pour sauvegarde
        if stat~=1
            % résultat sans stat et sans sensibilité
            RESULT=struct('NOM',[],'Poids',[],'Rna',[],'Prot',[]);
            
        else
            % résultat avec stat et sensibilité
            RESULT=struct('NOM',[],'Poids',[],'Rna',[],'Prot',[],'Stat',[],'Sens',[]);
            RESULT.Stat{1,1}={'domaine confiance' c1 'seuil' seuil2 'aire' Aire 'score' scoreS messageS};
            RESULT.Sens{1,1}={'sensibilité' X rho};
        end
        RESULT.NOM=NOM;
        
        RESULT.Poids{1,1}={fitPoids;messagePoids;errPoids;'erreur relative fit'};
        RESULT.Poids{1,2}=[{'verhulst' parv 'gompertz' parg 'contois' parc 'logistic' parl 'empirique' pare 'double sigmoide' pard 'poly' Wp3 Sp3 mp3 'noyau' h}];
        RESULT.Poids{1,3}=[xi,yi];
        
        RESULT.Rna{1,1}={fitR;messageR;errR;'erreur relative fit'};
        RESULT.Rna{1,2}=[{'poly deg3' ptrans3 St3 mt3 'poly deg 6' ptrans6 St6 mt6  'log poly deg3' plogtrans3 Slogt3 mlogt3 degre h}];
        RESULT.Rna{1,3}=[DPA_R,Rg];
        RESULT.Rna{1,4}=[{'moyenne',R0g}];
        
        RESULT.Prot{1,1}=[{'P0' kparg(1,1) 'ks' kparg(2,1)*P0g/R0g 'Ks' kparg(2,1) 'kd' kparg(3,1)}]
        RESULT.Prot{1,2}=[{{'erreur relative fit',errg(1)}; {messageOF(1) scoreOF(1);messageOF(2) scoreOF(2)}}];
        RESULT.Prot{1,3}=[DPA_P,Pg];
        RESULT.Prot{1,4}=[{'moyenne',P0g}];
        
        save(nom_Fichier_mat,'RESULT') ;
        % sauvegarde figures
        for j=cpt_fen_debut:cpt_fen
            nomfich=['graphes/' NOM num2str(j),'.pdf'];
            saveas(figure(j),nomfich);
        end
        % courbe pour publi
              nomfich=['graphes/' NOM 'RNA.svg'];
             saveas(figure(5001),nomfich)
             nomfich=['graphes/' NOM 'PROT.svg'];
             saveas(figure(5002),nomfich)
             %%%%%%%%%%%%%%%%
        
         prompt = 'fermeture O/N ?';
    dlg_title = 'fenetre';
    num_lines = 1;
    defaultans = {'O'};
    choix = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    
    if choix{1}=='O'
         close all %fermeture fenetre
        cpt_fen=0;            
    end
        
    end
    
end
end

