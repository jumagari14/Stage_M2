function [ fitR,message ] = choix_fitRNA
% choix de l'utilisateur du fit RNA
% par d�faut choix=5 fit RNA par exp(poly)

fitR=-1;
while fitR==-1
    prompt = 'Fit 1-poly degr�3 2-poly degr� 6 3-Noyau h=8 4-Noyau h=6 5-exp(poly) ';
    dlg_title = 'Choix courbe RNA';
    num_lines = 1;
    defaultans = {'5'};
    choix = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    switch choix{1}
        case '1'
            fitR=3;message='poly degr� 3';%poly degr�3
        case '2'
            fitR=6;message='poly degr� 6';%poly degr� 6
        case '3'
            fitR=0;message='estimateur � noyau h=8';%estimateur � noyau h=8
        case '4'
            fitR=1;message='estimateur noyau h=6';%estimateur � noyau h=6
        otherwise % exp(poly)
            fitR=2;message='exp(poly deg 3)';
    end;
end;



