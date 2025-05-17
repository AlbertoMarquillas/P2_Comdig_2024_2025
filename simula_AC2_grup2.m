% simula_AC2_grup2.m
% Simulació BER per codi lineal i codi Hamming en canal BSC

% CONSTANTS
% Definim el nombre de bits que simularem per cada valor de la probabilitat d’error del canal p. 
% Es divideix el total (4e5) en paquets de 100.000 bits per fer la simulació més eficient i modular.

nbitsPaquet = 100000;         % Nombre de bits per paquet
nbitsTotal = 4e5;             % Mínim de bits a simular per valor de p
nPaquets = ceil(nbitsTotal / nbitsPaquet); % Nombre de paquets a simular

% PROBABILITATS DE TRANSICIÓ (de 0.1 a 0.01, 9 valors)
% Calculem els 9 valors de la probabilitat d’error del canal BSC. Són valors logarítmicament 
% distribuïts entre 0.01 i 0.1. Aquesta p s’utilitzarà per simular els errors de transmissió.

k_vec = 0:8;
p = 10.^((k_vec - 16)/8);

% PARÀMETRES CODIS
% Es defineixen els paràmetres de codificació (k, n) per a cadascun dels dos codis:
% Codi lineal de (5,2)
% Codi Hamming de (7,4)

k_lin = 2;
n_lin = 5;
k_ham = 4;
n_ham = 7;

% MATRÍUS CODIS
% G_lineal és la matriu generadora del codi lineal general donada per l’enunciat.
% G_hamming s’obté amb la funció hammgen(3), que genera la matriu del codi Hamming amb paràmetre m=3 (i per tant (7,4)).

G_lineal = [1 0 1 1 0; 1 1 1 0 1]; % matriu G proporcionada
[~, G_hamming] = hammgen(3);      % codi Hamming (7,4)

% TAULA DE SÍNDROMES (només per codi lineal)
% Es calcula la taula de síndromes per al codi lineal a partir de la seva matriu generadora.
% Aquesta taula optimitza la descodificació fent-la més ràpida (lookup directe).

synd_lineal = syndtable(gen2par(G_lineal));

% RESULTATS BER
% Inicialitzem vectors per guardar la BER (bit error rate) de cada codi per a cada valor de p.
BER_hamming = zeros(1, length(p));
BER_lineal = zeros(1, length(p));

% INICI SIMULACIÓ
% Comencem el bucle principal de simulació. Iterem per cada valor de la probabilitat d’error p(i).
% Inicialitzem els comptadors d’errors a zero per aquest valor de p.
for i = 1:length(p)
    errH = 0;
    errL = 0;
    % Inicialitzem els comptadors d’errors a zero per aquest valor de p.
    % Es creen tantes paraules com calgui per arribar a nbitsPaquet en cada paquet.
    % Cada fila de mH i mL és una paraula de missatge de k bits.
    for j = 1:nPaquets
        % GENERACIÓ DE MISSATGES
        mH = randi([0 1], ceil(nbitsPaquet / k_ham), k_ham);
        mL = randi([0 1], ceil(nbitsPaquet / k_lin), k_lin);

        % CODIFICACIÓ
        % Codifiquem cada paraula utilitzant el codi corresponent.
        % encode retorna una matriu on cada fila és una paraula codificada de n bits.
        % Pel codi lineal s’especifica la matriu G manualment.
        cH = encode(mH, n_ham, k_ham, 'hamming/binary');
        cL = encode(mL, n_lin, k_lin, 'linear/binary', G_lineal);

        % TRANSMISSIÓ PEL CANAL BSC
        % Simulem el canal BSC: cada bit té una probabilitat p(i) de ser invertit.
        % Obtenim les paraules rebudes (amb errors) rH i rL.
        rH = bsc(cH, p(i));
        rL = bsc(cL, p(i));

        % DESCODIFICACIÓ
        % Descodifiquem els codis rebuts.
        % Per Hamming, MATLAB ja coneix la taula de síndromes.
        % Per al codi lineal passem explícitament la taula de síndromes (synd_lineal) per millorar la velocitat.
        dH = decode(rH, n_ham, k_ham, 'hamming/binary');  % sense taula
        dL = decode(rL, n_lin, k_lin, 'linear/binary', G_lineal, synd_lineal);

        % COMPTATGE D'ERRORS
        % Comptem quants bits descodificats no coincideixen amb els bits originals.
        % biterr compara les matrius i retorna el total d’errors.
        errH = errH + biterr(mH, dH);
        errL = errL + biterr(mL, dL);

        % Final del bucle que repeteix la simulació nPaquets vegades per acumular errors suficients per a una bona estimació estadística.
    end

    % BER FINAL PER AQUEST VALOR DE p
    % Es calcula la BER dividint els errors acumulats pel nombre total de bits transmesos.
    % BER_hamming(i) i BER_lineal(i) guarden el resultat per a aquest p(i).
    bits_H = ceil(nbitsPaquet / k_ham) * k_ham * nPaquets;
    bits_L = ceil(nbitsPaquet / k_lin) * k_lin * nPaquets;

    BER_hamming(i) = errH / bits_H;
    BER_lineal(i) = errL / bits_L;

    % Fi del bucle principal: s’ha simulat i calculat la BER per a tots els valors de p.
end

% PLOT FINAL
% Es representa la BER en un gràfic log-log.
% Es visualitzen les dues corbes amb colors diferents.
% Permet comparar visualment el rendiment dels dos codis a mesura que augmenta l’error del canal.
figure;
loglog(p, BER_hamming, '-ob', p, BER_lineal, '-or', 'LineWidth', 2);
xlabel('Probabilitat d''error del canal (p)');
ylabel('Probabilitat d''error després de decodificar');
legend('Hamming (7,4)', 'Lineal (5,2)', 'Location', 'southwest');
title('Simulació BER per codis de bloc en canal BSC');
grid on;

disp('BER Hamming:');
disp(BER_hamming);
disp('BER Lineal:');
disp(BER_lineal);

saveas(gcf, 'grafica_BER_AC2.png');