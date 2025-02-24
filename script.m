% Questo script creato da lcava000
clear all; clc;

% Richiesta dell'ordine del sistema
n = input('Inserisci l''ordine del sistema: ');

% Richiesta della matrice A
A = zeros(n, n);
for i = 1:n
    for j = 1:n
        A(i, j) = input(sprintf('Elemento A(%d,%d): ', i, j));
    end
end

% Richiesta della matrice B
B = zeros(n, 1);
for i = 1:n
    B(i, 1) = input(sprintf('Elemento B(%d,1): ', i));
end

% Richiesta della matrice C
C = zeros(1, n);
for j = 1:n
    C(1, j) = input(sprintf('Elemento C(1,%d): ', j));
end

%---------------------------

clc;
disp('Matrice A:');
disp(A);
disp('Matrice B:');
disp(B);
disp('Matrice C:');
disp(C);

correct = input('Le matrici sono corrette? (1 = Sì, 0 = No): ');
if correct == 0
    disp('Riavviare lo script per correggere le matrici.');
    return;
end

%---------------------------
while true
    clc;
    disp('MENU');
    disp('1. Verifica la stabilità del sistema');
    disp('2. Verifica la controllabilità del sistema');
    disp('3. Verifica l''osservabilità del sistema');
    disp('4. Progettazione del controllore usando Kalman');
    disp('5. Progettazione dell osservatore usando Kalman');
    disp('6. Progettazione del controllore usando Ackerman');
    disp('7. Progettazione dell osservatore usando Ackerman');
    disp('8. Calcola la funzione di Trasferimento G(S)');
    disp('9. Diagramma di Bode sulla G(s) calcolata');
    disp('10. Calcola la risposta al gradino');
    disp('11. Esci');
    scelta = input('Seleziona un''opzione: ');
    
    switch scelta
        case 1
            clc;
            eigenvalues = eig(A);
            disp('Autovalori di A:');
            disp(eigenvalues);
            if all(real(eigenvalues) < 0)
                disp('Il sistema è stabile.');
            else
                disp('Il sistema NON è stabile.');
            end
            input('Premi Invio per tornare al menu principale...');
        
        case 2
            clc;
            % Matrice di controllabilità
            Mc = ctrb(A, B);
            rankMc = rank(Mc);
            disp('Matrice di controllabilità Mc:');
            disp(Mc);
            if rankMc == n
                disp('Il sistema è controllabile.');
            else
                disp('Il sistema NON è completamente controllabile.');
            end
            input('Premi Invio per tornare al menu principale...');

        case 3
            clc;
            % Matrice di osservabilità
            Obsv = obsv(A, C);
            rankObsv = rank(Obsv);
            disp('Matrice di osservabilità:');
            disp(Obsv);
            if rankObsv == n
                disp('Il sistema è osservabile.');
            else
                disp('Il sistema NON è completamente osservabile.');
            end
            input('Premi Invio per tornare al menu principale...');
        
        
        case 4
            clc;
            % Verifica rango della matrice di controllabilità
            Mc = ctrb(A, B);
            rankMc = rank(Mc);
            n = size(A,1);
            disp('Matrice di controllabilità Mc:');
            disp(Mc);
            
            if rankMc == n
                disp('Il sistema è completamente controllabile, procediamo.');
                Ak2 = A;
                Bk2 = B;
                Ck2 = C;
            else
                disp('Il sistema NON è completamente controllabile. Separiamo la parte controllabile e non controllabile.');
                
                % Selezioniamo le colonne linearmente indipendenti di Mc
                [R, pivotCols] = rref(Mc);
                Xc = Mc(:, pivotCols);
                base = null(Xc', 'r');
                Xnc = base(:, 1);
                
                % Costruzione della matrice T
                T = [Xc Xnc];
                T_inv = inv(T);
                AK = T_inv * A * T;
                BK = T_inv * B;
                CK = C * T;
                
                % Estrazione delle sottomatrici controllabili
                Ak2 = AK(1:rankMc, 1:rankMc);
                Bk2 = BK(1:rankMc, :);
                Ck2 = CK(:, 1:rankMc);
            end
            
            disp('Sottomatrici controllabili:');
            disp('Ak2 ='); disp(Ak2);
            disp('Bk2 ='); disp(Bk2);
            disp('Ck2 ='); disp(Ck2);
            
            % Inserimento del polinomio desiderato
            poly_desired = input('Inserisci il polinomio caratteristico desiderato come vettore di coefficienti es: [1 10 25]: ');
            
            if length(poly_desired) - 1 ~= size(Ak2, 1)
                disp('Errore: Il grado del polinomio desiderato deve essere uguale alla dimensione di Ak2.');
            else

                % Creazione dinamica delle variabili simboliche in base alla dimensione di Ak2
                num_vars = size(Ak2, 1); % Numero di variabili = dimensione di Ak2
                syms_vars = sym('k', [1 num_vars]); % Crea k1, k2, ..., kn dinamicamente
                K = syms_vars; % Assegna le variabili simboliche a K
                
                
                % Costruzione corretta di Kc
                Kc = Ak2 - Bk2 * K;
                
                disp('Matrice Kc = Ak2 - Bk2 * K:');
                disp(Kc);
                
                % Polinomio caratteristico di Kc
                char_poly_Kc = charpoly(Kc);
                disp('Polinomio caratteristico di Kc:');
                disp(char_poly_Kc);
                
                % Risoluzione delle equazioni per k1 e k2
                eqs = char_poly_Kc == poly_desired;
                sol = solve(eqs, syms_vars); % Ora risolve automaticamente per tutti i k1, k2, ..., kn
                
                % Costruzione della matrice Kc con ordine n
                Kc = [double(sol.k1), double(sol.k2), zeros(1, n - 2)];
                
                % Calcolo di K finale
                if rankMc < n  % Se il sistema NON è completamente controllabile
                    disp('Mostro Kc:');
                    disp(Kc);
                    K = Kc * T_inv;
                    disp('Mostro K finale:');
                    disp(K);
                else  % Se il sistema è completamente controllabile
                    disp('Il sistema è completamente controllabile. Mostro direttamente Kc:');
                    disp(Kc);
                end
  
            end
            input('Premi Invio per tornare al menu principale...');

        case 5
            clc;
            Oc = obsv(A, C);
            rankOc = rank(Oc);
            n = size(A,1);
            disp('Matrice di osservabilità Oc:');
            disp(Oc);
            
            if rankOc == n
                disp('Il sistema è completamente osservabile.');
                Ao2 = A;
                Bo2 = B;
                Co2 = C;
                T = eye(n);
            else
                disp('Separazione della parte osservabile.');
                
                % Trova le righe linearmente indipendenti di Oc
                [~, pivotRows] = rref(Oc);  
                XoT = Oc(pivotRows, :); % Prende le righe linearmente indipendenti
                
                % Trova il vettore mancante per rendere T invertibile
                Xno = null(XoT, 'r'); % Trova una base per il complemento ortogonale
                
                if isempty(Xno)
                    error('Errore: impossibile trovare una base per X_NO.');
                end
                
                % Se Xno ha più colonne, ne seleziono una valida
                Xno = Xno(:, 1); 
                
                % Costruzione della matrice T
                T = [XoT; Xno']; 
                
                % Verifica che T sia invertibile
                if det(T) == 0
                    error('Errore: la matrice T è singolare, il sistema potrebbe essere mal condizionato.');
                end
                
                T_inv = inv(T);
                
                disp('Matrice T corretta:');
                disp(T);
                disp('Matrice T inversa:');
                disp(T_inv);
                
                % Trasformazione di coordinate
                Ao = T_inv * A * T;
                Bo = T_inv * B;
                Co = C * T;
                
                % Estrazione delle sottomatrici osservabili
                Ao2 = Ao(1:rankOc, 1:rankOc);
                Bo2 = Bo(1:rankOc, :);
                Co2 = Co(:, 1:rankOc);
            end

            disp("Ao1");
            disp(Ao2);
            disp("Co2");
            disp(Co2);
       
          poly_desired = input('Inserisci il polinomio caratteristico desiderato per l\osservatore es: [1 10 25]: ');

        if length(poly_desired) - 1 ~= size(Ao2, 1)
            disp('Errore: Il grado del polinomio desiderato deve essere uguale alla dimensione di Ao2.');
        else
            syms l [size(Ao2,1) 1]; % Vettore colonna simbolico
            Lc = Ao2 - l * Co2; % Mantiene la moltiplicazione corretta
            
            char_poly_Lc = charpoly(Lc);
            eqs = char_poly_Lc == poly_desired;
            sol = solve(eqs, l); % Risolve automaticamente per l1, l2, ..., ln
            
            % Costruzione di Lc come vettore colonna
            Lc = zeros(size(Ao2,1), 1); 
            for i = 1:size(Ao2,1)
                fieldname = sprintf('l%d', i); % Nome del campo simbolico (l1, l2, ...)
                Lc(i) = double(sol.(fieldname)); % Assegna il valore numerico
            end
            
            % AGGIUNGI ZERI IN MODO DINAMICO SE SERVE
            diff_size = n - size(Lc,1); % Differenza tra n e righe attuali di Lc
            if diff_size > 0
                Lc = [Lc; zeros(diff_size, 1)]; % Aggiunge il numero corretto di zeri
            end
            
            % Calcolo di L finale
            if rankOc < n
                L = T * Lc; % Moltiplicazione corretta
                disp('L finale:');
                disp(L);
            else
                disp('Lc:');
                disp(Lc);
            end
        end
        input('Premi Invio per tornare al menu principale...');
        
case 6
            clc;
            disp('Progettazione del regolatore con la formula di Ackermann.');

            % Calcola la matrice di controllabilità
            Mc = ctrb(A, B);
            rankMc = rank(Mc);
            disp('Matrice di controllabilità Mc:');
            disp(Mc);
            
            % Controllabilità del sistema
            if rankMc ~= size(A,1)
                disp('Errore: Il sistema NON è completamente controllabile. Ackermann non è applicabile.');
                input('Premi Invio per tornare al menu principale...');
                return;
            end
            
            % Chiede il polinomio desiderato all'utente
            poly_desired = input('Inserisci il polinomio caratteristico desiderato es: [1 10 75 125]: ');

            if length(poly_desired) - 1 ~= size(A,1)
                disp('Errore: Il grado del polinomio desiderato deve essere uguale alla dimensione di A.');
                input('Premi Invio per tornare al menu principale...');
                return;
            end
            
            % Calcola Pd(A) = A^n + c1*A^(n-1) + ... + cn*I
            n = size(A,1);
            Pd_A = zeros(n);
            for i = 1:length(poly_desired)
                Pd_A = Pd_A + poly_desired(i) * (A^(length(poly_desired)-i));
            end
            
            % Calcola la matrice di Ackermann
            try
                Mc_inv = inv(Mc); % Inversa di Mc
            catch
                disp('Errore: Impossibile invertire Mc, il sistema potrebbe essere mal condizionato.');
                input('Premi Invio per tornare al menu principale...');
                return;
            end
            
            % Costruzione della matrice [0 0 1] per l'ultima riga
            vec = zeros(1, n);
            vec(end) = 1; % Ultima colonna = 1

            disp("Pd(A)");
            disp(Pd_A);
            disp("Mc Inversa");
            disp(Mc_inv);
            
            % Applicazione della formula di Ackermann
            K = vec * Mc_inv * Pd_A;

            
            % Mostra il risultato
            disp('Matrice di guadagno K:');
            disp(K);

            input('Premi Invio per tornare al menu principale...');


        case 7
            clc;
            disp('Progettazione dell\osservatore con la formula di Ackermann.');

            % Calcola la matrice di osservabilità
            Mo = obsv(A, C);
            rankMo = rank(Mo);
            disp('Matrice di osservabilità Mo:');
            disp(Mo);
            
            % Osservabilità del sistema
            if rankMo ~= size(A,1)
                disp('Errore: Il sistema NON è completamente osservabile. Ackermann non è applicabile.');
                input('Premi Invio per tornare al menu principale...');
                return;
            end
            
            % Chiede il polinomio desiderato all'utente
            poly_desired = input('Inserisci il polinomio caratteristico desiderato es: [1 10 75 125]: ');

            if length(poly_desired) - 1 ~= size(A,1)
                disp('Errore: Il grado del polinomio desiderato deve essere uguale alla dimensione di A.');
                input('Premi Invio per tornare al menu principale...');
                return;
            end
            
            % Calcola Pd(A) = A^n + c1*A^(n-1) + ... + cn*I
            n = size(A,1);
            Pd_A = zeros(n);
            for i = 1:length(poly_desired)
                Pd_A = Pd_A + poly_desired(i) * (A^(length(poly_desired)-i));
            end
            
            % Calcola la matrice di Ackermann per l'osservatore
            try
                Mo_inv = inv(Mo); % Inversa di Mo
            catch
                disp('Errore: Impossibile invertire Mo, il sistema potrebbe essere mal condizionato.');
                input('Premi Invio per tornare al menu principale...');
                return;
            end
            
            % Costruzione della matrice [0; 0; 1] per l'ultima riga
            vec = zeros(n, 1);
            vec(end) = 1; % Ultima riga = 1

            disp("Pd(A)");
            disp(Pd_A);
            disp("Mo Inversa");
            disp(Mo_inv);
            
            % Applicazione della formula di Ackermann per l'osservatore
            L = Pd_A * Mo_inv * vec;

            % Mostra il risultato
            disp('Matrice di guadagno L:');
            disp(L);

            input('Premi Invio per tornare al menu principale...');

        case 8
            % Calcolo della funzione di trasferimento passo per passo
            clc; 
            syms s;
            
            % Calcolo della funzione di trasferimento simbolicamente
            I = eye(size(A));
            SI_A = s*I - A;
            disp('Matrice (sI - A):');
            disp(SI_A);
            
            inv_SI_A = simplify(inv(SI_A));
            disp('Matrice inversa di (sI - A):');
            disp(inv_SI_A);
            
            TF_matrix = simplify(C * inv_SI_A * B);
            disp('Matrice della funzione di trasferimento:');
            disp(vpa(TF_matrix, 4));
            
            % Creazione del sistema dinamico e calcolo della funzione di trasferimento
            sys = ss(A, B, C, 0);
            FDT = tf(sys);
            
            % Estrazione dei coefficienti del numeratore e denominatore
            [num, den] = tfdata(FDT, 'v'); 
            
            % Stampa della funzione di trasferimento in formato leggibile
            disp('Funzione di trasferimento finale:');
            fprintf('G(s) = (%s) / (%s)\n', poly2str(num, 's'), poly2str(den, 's'));
            
            input('Premi Invio per tornare al menu principale...');

       case 9
            clc; 
            if exist('FDT', 'var')
                disp('Funzione di trasferimento già calcolata. Mostrando il diagramma di Bode...');
                figure;
                bode(FDT);
                grid on;
                title('Diagramma di Bode della Funzione di Trasferimento');
                input('Premi Invio per tornare al menu principale...');
            else
                disp('Errore: La funzione di trasferimento G(s) non è stata ancora calcolata.');
                input('Premi Invio per tornare al menu principale...');
            end

        case 10
            clc; % Pulisce la console
            
            % Verifica se la variabile FDT esiste
            if exist('FDT', 'var')
                syms t s; % Definisce le variabili simboliche t e s
                
                disp('Calcolo della risposta al gradino...');
                
                % Visualizza la risposta al gradino utilizzando la funzione step
                step(FDT);
                
                % Costruisce la funzione di trasferimento G(s) utilizzando i coefficienti num e den
                Gs = poly2sym(num, s) / poly2sym(den, s);
                
                % Definisce l'ingresso al gradino nel dominio delle frequenze (U(s) = 1/s)
                U_s = 1/s;
                
                % Calcola la uscita nel dominio delle frequenze (Y(s) = G(s) * U(s))
                Y_s = simplify(Gs * U_s);
                
                % Decomponi Y(s) in frazioni parziali per facilitare la trasformata inversa
                Y_s = simplify(partfrac(Y_s));
                
                % Calcola la risposta temporale y(t) tramite la trasformata inversa di Laplace
                Y_t = ilaplace(Y_s, s, t);
                
                % Mostra la risposta al gradino in funzione del tempo
                disp('Risposta al gradino in funzione del tempo:');
                fprintf('y(t) = %s\n', char(vpa(Y_t, 4)));
                
                % Attende un input dell'utente prima di tornare al menu principale
                input('Premi Invio per tornare al menu principale...');
            else
                % Gestisce il caso in cui la funzione di trasferimento G(s) non sia stata calcolata
                disp('Errore: La funzione di trasferimento G(s) non è stata ancora calcolata.');
                input('Premi Invio per tornare al menu principale...');
            end


        case 11
            disp('Uscita dal programma.');
            break;
        
        otherwise
            disp('Opzione non valida, riprova.');
            pause;
    end
end