clear; clc; close all;
fprintf('=== Tree Crown Segmentation (Erikson FAST) — START ===\n');
tTotal = tic;

%% 1. CARICAMENTO E PRE-PROCESSING
%  ================================================================================================================================================

% Path dei dati: Sostituisci con i tuoi percorsi locali
% Carica l'immagine aerea multispettrale originale
imgPath = '/Users/mariocassano/Desktop/POLIBA/Image Processing/Tree-Crown-Segmentation/data/esempio3/2019_WREF_3_579000_5082000_image.tif';

chmPath = '/Users/mariocassano/Desktop/POLIBA/Image Processing/Tree-Crown-Segmentation/data/esempio3/NEON_D16_WREF_DP3_579000_5082000_CHM.tif';

shpPath = '/Users/mariocassano/Desktop/POLIBA/Image Processing/Tree-Crown-Segmentation/data/esempio3/2019_WREF_3_579000_5082000_image.shp';

img = imread(imgPath);
[origH, origW, bands] = size(img); % Estrae le dimensioni spaziali e il numero di bande

imgDouble = double(img); % Converte in double per prevenire overflow durante le operazioni matematiche

% Carica il Canopy Height Model (CHM) che contiene i dati di altezza della vegetazione
chmFull = imread(chmPath);
chmFull(chmFull < 0) = 0; % Elimina i valori anomali (es. nodata negativi o artefatti del lidar)
chmDouble = double(chmFull);

% ridimensiona il CHM per farlo combaciare pixel-per-pixel con l'immagine ottica.
% L'uso di 'nearest' (nearest-neighbor) evita di interpolare e creare falsi valori di altezza.
chmDouble_Resized = imresize(chmDouble, [origH, origW], 'nearest');

% Crea una maschera logica per ignorare i pixel di background (completamente neri in tutte le bande)
% Li ignori creando una maschera che vale 1 se il pixel è > 0 e vale 0 se è tutto nero
validMask = sum(imgDouble, 3) > 0; 
if sum(validMask(:)) == 0
    error('L''immagine è completamente nera!'); % Controllo di sicurezza
end

% Estrae solo i pixel validi per calcolare statistiche robuste
validPixels = imgDouble(repmat(validMask, [1 1 bands]));

% Calcola il 1° e il 99° percentile per eseguire un contrast stretching che ignori gli outlier
low_in = prctile(validPixels, 1);
high_in = prctile(validPixels, 99);

% Loop su ogni banda per normalizzare i valori nel range [0, 1]
for b = 1:bands
    band = imgDouble(:,:,b);
    band(band < low_in) = low_in;   % Clipping dei valori inferiori al 1° percentile
    band(band > high_in) = high_in; % Clipping dei valori superiori al 99° percentile
    imgDouble(:,:,b) = (band - low_in) / (high_in - low_in); % Normalizzazione Min-Max
end

% --- CROP INTERATTIVO ---
useSubset = true; 
if useSubset
    fprintf('Seleziona un''area di foresta densa disegnando un rettangolo sull''immagine...\n');
    
    % Mostra l'immagine normalizzata in una finestra massimizzata per facilitare la selezione
    figCrop = figure('Name', 'Seleziona l''area di test', 'WindowState', 'maximized');
    imshow(imgDouble);
    title('DISEGNA UN RETTANGOLO su un''area di foresta (poi fai doppio click per confermare)');
    
    rect = getrect(figCrop); % Acquisisce le coordinate del rettangolo disegnato dall'utente
    close(figCrop);
    
    % Ritaglia contemporaneamente l'immagine originale, quella normalizzata e il CHM
    % garantendo il perfetto allineamento spaziale dei dati multispettrali e altimetrici
    img = imcrop(img, rect);       
    imgDouble = imcrop(imgDouble, rect); 
    chmCrop = imcrop(chmDouble_Resized, rect); % Ora usiamo il CHM già ingrandito!
    
    fprintf('Crop applicato: %d x %d pixel.\n', size(img,1), size(img,2));
else
    % Fallback nel caso in cui si decida di elaborare l'intera immagine
    fprintf('Elaborazione dell''INTERA IMMAGINE in corso...\n');
    chmCrop = chmDouble_Resized;
    cMin = 1; 
    rMin = 1; 
    cMax = origW;
    rMax = origH;
end

[rows, cols, ~] = size(imgDouble); % Aggiorna le dimensioni dopo il potenziale ritaglio

% Applica un filtro Gaussiano con deviazione standard 0.5 per ridurre il rumore ad alta frequenza
imgDouble = imgaussfilt(imgDouble, 0.5);

% Estrazione della prima banda (Near-Infrared), fondamentale per i calcoli successivi
% nir = imgDouble(:,:,1);

%% --- INTEGRAZIONE ICELab  ---
% Converto l'immagine normalizzata in Lab
imgLab = rgb2lab(imgDouble); 

% Estraggo la componente 'a' (canale 2)
% in MATLAB rgb2lab restituisce 'a' nel range circa [-100, 100]
comp_a = imgLab(:,:,2); 

comp_a_filtered = imgaussfilt(comp_a, 1.5); % Leviga il rumore spettrale

% Normalizzazione invertita
% Utilizzo 'comp_a'
comp_a_filtered_norm = (max(comp_a(:)) - comp_a) / (max(comp_a(:)) - min(comp_a(:)));

% Creazione banda ibrida per il growing
% Uso la banda 2 (Green) di imgDouble e la componente 'a'
imgForGrowing = (0.7 * imgDouble(:,:,2)) + (0.3 * comp_a_filtered_norm);

%% 2. DEFINIZIONE DEI PARAMETRI DI SEGMENTAZIONE
%  ========================================================================
ALPHA = 0.10;       
SIGMA1 = 0.16;
MIN_AREA = 30;      
MAX_BOX_RADIUS = 40; 

MIN_SEED_DIST = 20; 
ALTEZZA_MINIMA = 2.0; 
sigma_seeds = 2; 
sigma2_init = 4; 
ratioTreeLimit = 0.85; 
%% 3. INDIVIDUAZIONE SEED POINTS (PICCHI LUMINOSI)
%  ========================================================================
fprintf('--- Finding Seed Points ---\n');

% A Maschera vegetazione (Otsu)
% Calcola la soglia ottimale in modo automatico invece di usare un valore fisso
% è la linea di confine tra gli oggetti chiari, gli alberi, e gli oggetti scuri come lo sfondo, ecc
level = graythresh(imgForGrowing); 
nirThresh = level * 0.9; % Rilassa leggermente la soglia per non perdere i bordi delle chiome, che sono più scuri
binaryTree = imgForGrowing > nirThresh; % Binarizzazione di ogni pixel: 1 per la vegetazione potenziale, 0 per il resto
binaryTree = imopen(binaryTree, strel('disk', 1)); % Operazione morfologica per rimuovere piccoli rumori isolati che possono ingannare l'algoritmo

% B) MASCHERA ALTIMETRICA DAL CHM
maskAltezza = chmCrop >= ALTEZZA_MINIMA; % Mantiene solo i pixel che superano la soglia di altezza (es. 2m)

% FUSIONE DELLE MASCHERE: Deve essere verde (Otsu) E alto almeno 2m (CHM)
binaryTree = binaryTree & maskAltezza; % Intersezione logica (AND) per una robustezza estrema
% ----------------------------------------------

% Creiamo una versione molto sfocata dell'infrarosso (NIR)
% Questo serve a fondere le piccole variazioni all'interno di una chioma in un unico grande "blob"
% trasforma i tanti picchi luminosi sull'albero mescolandoli e
% raggruppandoli in pochi blob sfocati (molto luminosi al centro della
% matrice, sempre più sfocati verso l'esterno)
nir_smoothed = imgaussfilt(imgForGrowing, sigma_seeds);

% Trova i picchi (Cime degli alberi)
% imextendedmax trova i massimi regionali, filtrando i picchi spuri inferiori alla soglia 0.06
seedsBinary = imextendedmax(nir_smoothed, 0.06); 
seedsBinary = seedsBinary & binaryTree; % Mantiene solo i picchi che cadono dentro la maschera (Verdi + Alti)

% Estrae le coordinate (righe, colonne) dei picchi trovati
[seedRows, seedCols] = find(seedsBinary);

if ~isempty(seedRows)
    seedVals = sub2ind(size(imgForGrowing), seedRows, seedCols);
    seedVals = imgForGrowing(seedVals); 
    
    % Ordina i semi dal più luminoso (più probabile sia una cima) al meno luminoso
    [~, sortIdx] = sort(seedVals, 'descend');
    seedList = [seedRows(sortIdx), seedCols(sortIdx)];
else
    seedList = [];
    warning('Nessun seed trovato. Prova ad abbassare il valore in imextendedmax.');
end
fprintf('  Seeds trovati: %d\n', length(seedList));

% --- FILTRO SPAZIALE: DISTANZA MINIMA TRA I SEMI ---
filteredSeedList = [];
fprintf('  Applicazione filtro spaziale (Distanza minima: %d px)...\n', MIN_SEED_DIST);

% Itera sulla lista ordinata per applicare una "Non-Maximum Suppression" spaziale
for i = 1:size(seedList, 1)
    pt = seedList(i, :);
    if isempty(filteredSeedList)
        filteredSeedList = [filteredSeedList; pt]; %#ok<AGROW,AGROW> % Accetta il primo seme (il più luminoso in assoluto)
    else
        % Calcola la distanza euclidea tra il seme corrente e tutti i semi già accettati
        dists = sqrt((filteredSeedList(:,1) - pt(1)).^2 + (filteredSeedList(:,2) - pt(2)).^2);
        
        % Se il seme corrente è sufficientemente lontano da TUTTI quelli già accettati, lo mantiene
        if min(dists) >= MIN_SEED_DIST
            filteredSeedList = [filteredSeedList; pt]; %#ok<AGROW,AGROW>
        end
    end
end
seedList = filteredSeedList;
nSeeds = size(seedList, 1);
fprintf('  Seeds rimasti dopo il filtro spaziale: %d\n', nSeeds);

% 4. LOOP PRINCIPALE
%  ========================================================================
% Mappa finale per le etichette dei singoli alberi (Connected Components)
finalLabelMap = zeros(rows, cols);
currentLabel = 0;

% Maschera logica per tracciare i pixel già assegnati ed evitare sovrapposizioni
processedMask = false(rows, cols); 
hWait = waitbar(0, 'Segmentazione in corso...');
nSeeds = size(seedList, 1);

for k = 1:nSeeds
    if mod(k, 20) == 0
        waitbar(k/nSeeds, hWait, sprintf('Albero %d / %d', k, nSeeds));
    end
    
    startPt = seedList(k, :);
    
    % Se il punto di partenza cade in un albero già processato, saltalo
    if processedMask(startPt(1), startPt(2))
        continue;
    end
    
    % --- L'ALBERO APPROSSIMATO ---
    regionInit = growRegionFast(imgForGrowing, startPt, SIGMA1, sigma2_init, ALPHA, processedMask, MAX_BOX_RADIUS);
    
    if sum(regionInit(:)) < 10
        continue; 
    end
    
    % --- STEP 3: STIMA DEI PARAMETRI ---
    % Semplificazione vettorializzata: stima del raggio rho usando il diametro equivalente

    % calcola diametro e centroide delle macchie iniziali (regionInit)
    props = regionprops(regionInit, 'EquivDiameter', 'Centroid');
    if isempty(props)
        continue; 
    end

    % stima del raggio rho
    rho = props(1).EquivDiameter / 2;
    
    % Calcolo di sigma_2 inversa derivato dall'Eq. 6 del paper
    denom = sqrt(-2 * log(ALPHA));
    sigma2_est = (2 * rho) / denom;
    sigma2_est = min(sigma2_est*1.3, 40); % Tappo di sicurezza per non far esplodere la regione
    
    % --- RICERCA NUOVI STARTING POINTS ---
    % Invece di calcolare la normale al contorno, si crea una griglia di ricerca attorno al centroide
    cent = props(1).Centroid; 
    searchRad = max(2, round(rho * 0.5));
    
    [cGrid, rGrid] = meshgrid( -searchRad:2:searchRad, -searchRad:2:searchRad );
    candC = round(cent(1)) + cGrid(:);
    candR = round(cent(2)) + rGrid(:); 
    
    % Pulisce i punti fuori dall'immagine
    valid = candR>0 & candR<=rows & candC>0 & candC<=cols;
    candR = candR(valid);
    candC = candC(valid);
    
    % Ottimizzazione delle performance: testa al massimo 15 candidati casuali
    MAX_ATTEMPTS = 15;
    if length(candR) > MAX_ATTEMPTS
        idx = randperm(length(candR), MAX_ATTEMPTS);
        candR = candR(idx);
        candC = candC(idx);
    end
    
    bestMetric = Inf;
    bestRegion =[];
    foundCandidate = false;
    
    % --- STEP 5: CREAZIONE REGIONI CANDIDATE ---
    for i = 1:length(candR)
        candPt = [candR(i), candC(i)];
        
        if processedMask(candPt(1), candPt(2))
            continue;
        end
        
        candRegion = growRegionFast(imgForGrowing, candPt, SIGMA1, sigma2_est, ALPHA, processedMask, MAX_BOX_RADIUS);
        
        % --- STEP 6: VERIFICA DELLE "TREE CONDITIONS" ---
        area = sum(candRegion(:));
        if area < MIN_AREA
            continue; % Condizione 1 adattata ai propri dati
        end
        
        intersectTree = candRegion & binaryTree;
        ratioTree = sum(intersectTree(:)) / area;
        if ratioTree < ratioTreeLimit
            continue;
        end
        
        % --- STEP 6: SELEZIONE DELLA REGIONE MIGLIORE ---
        % Calcola la media dei valori Infrarosso (NIR) sul perimetro della regione
        perim = bwperim(candRegion);
        meanContourVal = mean(imgForGrowing(perim));
        
        % Vince la regione con il contorno più scuro
        if meanContourVal < bestMetric
            bestMetric = meanContourVal;
            bestRegion = candRegion;
            foundCandidate = true;
        end
    end
    
    % RIMOZIONE DEI PUNTI E AGGIORNAMENTO MASCHERA ---
    if foundCandidate
        currentLabel = currentLabel + 1;
        finalLabelMap(bestRegion) = currentLabel;
        processedMask(bestRegion) = true; % Rimuove l'albero trovato dal pool
    else
        processedMask(startPt(1), startPt(2)) = true; % Invalida solo il seme fallato
    end
end
close(hWait);
fprintf('  Segmentazione completata. Alberi trovati: %d\n', currentLabel);

%% 5. VISUALIZZAZIONE
%  ========================================================================
% Crea una finestra sufficientemente larga per ospitare comodamente i 3 subplot
figure('Name', 'Erikson Fast Results', 'Position',[100 100 1200 500]);

% --- PANNELLO 1: IMMAGINE ORIGINALE ---
subplot(1,3,1);
imshow(img); title('Originale (Crop)');
imwrite(img, 'img_big.png'); % Salva il dato grezzo per report o paper

% --- PANNELLO 2: MAPPA DELLE ETICHETTE (SEGMENTI) ---
subplot(1,3,2);
% Converte la mappa dei label (dove ogni albero ha un ID intero) in un'immagine RGB.
% 'jet' è la colormap, 'k' assegna il nero allo sfondo (label 0),
% 'shuffle' mescola i colori per massimizzare il contrasto tra alberi adiacenti.
rgb = label2rgb(finalLabelMap, 'jet', 'k', 'shuffle');
imshow(rgb); title(['Segmentazione (N=', num2str(currentLabel), ')']);
imwrite(finalLabelMap, 'segmentation_img_big.png'); 

% --- PANNELLO 3: CONTORNI SOVRAPPOSTI ---
subplot(1,3,3);
% Estrae i bordi esterni della maschera binaria (alberi > 0)
% e li sovrappone all'immagine originale usando il colore rosso RGB [1 0 0]
B = imoverlay_custom(img, bwperim(finalLabelMap > 0), [1 0 0]);
imshow(B); title('Contorni');
imwrite(B, 'contorni_img_big.png');

%% --- PREPARAZIONE DATI PER LA VALIDAZIONE ---
fprintf('Preparazione dei dati per la validazione...\n');

% Se abbiamo fatto il crop, ricalcoliamo le coordinate assolute (bounding box) 
% rispetto all'immagine originale per mantenere il riferimento spaziale.
if useSubset
    cMin = max(1, round(rect(1)));
    rMin = max(1, round(rect(2)));
    cMax = cMin + size(finalLabelMap, 2) - 1;
    rMax = rMin + size(finalLabelMap, 1) - 1;
else
    cMin = 1;
    rMin = 1;
    cMax = size(finalLabelMap, 2);
    rMax = size(finalLabelMap, 1);
end

% Percorsi dei file: immagine raster originale e Ground Truth (GT) vettoriale

% Estrae l'oggetto di referenziazione spaziale (R) necessario per allineare pixel e coordinate geografiche
[~, R] = readgeoraster(imgPath);
truthShapes = shaperead(shpPath); % Carica i poligoni tracciati a mano (Ground Truth)

% Salva l'ambiente di validazione per analisi offline o debug successivi
save('dati_validazione.mat', 'finalLabelMap', 'cMin', 'rMin', 'cMax', 'rMax', 'R', 'truthShapes', 'chmCrop');
fprintf('Dati salvati con successo in dati_validazione.mat!\n');

%% 6. VALIDAZIONE AUTOMATICA SULL'AREA APPENA RITAGLIATA
fprintf('\n=== AVVIO VALIDAZIONE ===\n');

% Estrae i centroidi delle chiome stimate dal tuo algoritmo
props = regionprops(finalLabelMap, 'Centroid', 'Area');
alberiValidi = [props.Area] >= 30; % Filtro di sicurezza aggiuntivo (coerente con MIN_AREA)
predCentroids = cat(1, props(alberiValidi).Centroid);
predX = predCentroids(:,1); 
predY = predCentroids(:,2);

% COSTRUZIONE GROUND TRUTH LOCALE (Sporco vs Pulito)
dirtyGTPolygons = []; cleanGTPolygons = [];
for i = 1:length(truthShapes)
    polyX_W = truthShapes(i).X; polyY_W = truthShapes(i).Y;
    
    % Converte le coordinate geografiche (es. UTM) in coordinate intrinseche (pixel)
    [xPix_full, yPix_full] = worldToIntrinsic(R, polyX_W, polyY_W);
    xPix_full(isnan(xPix_full)) = []; yPix_full(isnan(yPix_full)) = [];
    
    % Controlla se il poligono del GT cade all'interno del crop che hai selezionato
    if any(xPix_full >= cMin & xPix_full <= cMax) && any(yPix_full >= rMin & yPix_full <= rMax)
        
        % Traslazione delle coordinate: adatta il poligono alla mini-immagine croppata
        xPix_crop = xPix_full - cMin + 1;
        yPix_crop = yPix_full - rMin + 1;
        
        poligono = polyshape(xPix_crop, yPix_crop);
        dirtyGTPolygons = [dirtyGTPolygons; poligono]; % Salva nel GT "Sporco" (tutti i poligoni)
        
        % Estrae il bounding box del poligono per analizzarne l'altezza sul CHM
        boxMinX = max(1, floor(min(xPix_crop))); boxMaxX = min(size(chmCrop, 2), ceil(max(xPix_crop)));
        boxMinY = max(1, floor(min(yPix_crop))); boxMaxY = min(size(chmCrop, 1), ceil(max(yPix_crop)));
        chmPatch = chmCrop(boxMinY:boxMaxY, boxMinX:boxMaxX);
        
        % IL FILTRO "PULITO": Se dentro il poligono manuale c'è almeno un pixel alto >= 2 metri, 
        % allora è un VERO albero e viene salvato nel GT "Pulito". Altrimenti era un errore di annotazione.
        if max(chmPatch(:)) >= 2.0
            cleanGTPolygons = [cleanGTPolygons; poligono];
        end
    end
end

% Stampa a schermo il conteggio
fprintf('Poligoni "Sporchi" (Totali): %d\n', length(dirtyGTPolygons));
fprintf('Poligoni "Puliti" (Filtrati): %d\n', length(cleanGTPolygons));

% Verifica automatica
if length(dirtyGTPolygons) == length(cleanGTPolygons)
    disp('RISULTATO: Il filtro non ha eliminato nulla.');
else
    scartati = length(dirtyGTPolygons) - length(cleanGTPolygons);
    disp(['RISULTATO: Il filtro ha scartato ', num2str(scartati), ' poligoni.']);
end

% CALCOLO METRICHE
fprintf('\n--- RISULTATI AREA ANALIZZATA ---\n');
% Calcola metriche usando il GT non filtrato
[TP_d, FP_d, FN_d, tpPts_d, fpPts_d, matchedPolys_d, fnPolys_d, Met_d] = calcolaMetriche(predX, predY, dirtyGTPolygons);
fprintf('SPORCA -> Precision: %.1f%% | Recall: %.1f%% | F1: %.1f%%\n', Met_d(1)*100, Met_d(2)*100, Met_d(3)*100);

% Calcola metriche usando il GT filtrato dal CHM
[TP_c, FP_c, FN_c, tpPts_c, fpPts_c, matchedPolys_c, fnPolys_c, Met_c] = calcolaMetriche(predX, predY, cleanGTPolygons);
fprintf('PULITA -> Precision: %.1f%% | Recall: %.1f%% | F1: %.1f%%\n', Met_c(1)*100, Met_c(2)*100, Met_c(3)*100);



% GRAFICI FINALI
figRes = figure('Name', 'Confronto Validazione', 'WindowState', 'maximized');
ax1 = subplot(1, 2, 1);
disegnaMappa(ax1, img, 1, 1, matchedPolys_d, fnPolys_d, tpPts_d, fpPts_d, ...
    sprintf('SPORCA (F1: %.1f%%)', Met_d(3)*100));
ax2 = subplot(1, 2, 2);
disegnaMappa(ax2, img, 1, 1, matchedPolys_c, fnPolys_c, tpPts_c, fpPts_c, ...
    sprintf('PULITA (F1: %.1f%%)', Met_c(3)*100));
linkaxes([ax1, ax2], 'xy'); % Sincronizza lo zoom/pan sui due subplot

%% ========================================================================
%%  FUNZIONI HELPER
%% ========================================================================
function mask = growRegionFast(img, seed, s1, s2, alpha, globalMask, maxRad)
    [H, W, channels] = size(img); % <-- Controlla i canali
    r0 = seed(1); 
    c0 = seed(2);
    
    rMin = max(1, r0 - maxRad); 
    rMax = min(H, r0 + maxRad);
    cMin = max(1, c0 - maxRad); 
    cMax = min(W, c0 + maxRad);
    
    imgCrop = img(rMin:rMax, cMin:cMax, :);
    
    lr0 = r0 - rMin + 1;
    lc0 = c0 - cMin + 1;
    
    [colsGrid, rowsGrid] = meshgrid(1:size(imgCrop,2), 1:size(imgCrop,1));
    distSq = (rowsGrid - lr0).^2 + (colsGrid - lc0).^2;
    mu2 = exp(-0.5 * distSq / (s2^2));
    
    % --- LOGICA ADATTIVA PER I CANALI ---
    if channels == 3
        seedColor = reshape(img(r0, c0, :), [1, 1, 3]);
        diffCol = imgCrop - seedColor;
        colDistSq = sum(diffCol.^2, 3);
    else
        seedColor = img(r0, c0);
        diffCol = imgCrop - seedColor;
        colDistSq = diffCol.^2;
    end
    % ------------------------------------
    
    mu1 = exp(-0.5 * colDistSq / (s1^2));
    
    muTotal = mu1 .* mu2;
    localMask = muTotal > alpha;
    
    localMask = bwselect(localMask, lc0, lr0, 8);
    localMask = imfill(localMask, 'holes');
    
    globalCrop = globalMask(rMin:rMax, cMin:cMax);
    localMask = localMask & ~globalCrop;
    
    mask = false(H, W);
    mask(rMin:rMax, cMin:cMax) = localMask;
end

function out = imoverlay_custom(in, mask, color)
    in = im2double(in);
    mask = logical(mask);
    out = in;
    for k = 1:3
        channel = out(:,:,k);
        channel(mask) = color(k);
        out(:,:,k) = channel;
    end
end

function [TP, FP, FN, tpPts, fpPts, matchedPolys, fnPolys, Metrics] = calcolaMetriche(predX, predY, gtPolygons)
    numPred = length(predX); numGT = length(gtPolygons);
    matchedGT = false(numGT, 1); isTP = false(numPred, 1);
    
    for g = 1:numGT
        vX = gtPolygons(g).Vertices(:,1); vY = gtPolygons(g).Vertices(:,2);
        minVx = min(vX); maxVx = max(vX); minVy = min(vY); maxVy = max(vY);
        
        cands = find(predX >= minVx & predX <= maxVx & predY >= minVy & predY <= maxVy);
        if ~isempty(cands)
            in = inpolygon(predX(cands), predY(cands), vX, vY);
            matchedIdx = cands(in);
            if ~isempty(matchedIdx)
                for m = 1:length(matchedIdx)
                    pidx = matchedIdx(m);
                    if ~isTP(pidx)
                        isTP(pidx) = true; matchedGT(g) = true; 
                        break;
                    end
                end
            end
        end
    end
    
    TP = sum(isTP); FP = numPred - TP; FN = sum(~matchedGT);
    tpPts = [predX(isTP), predY(isTP)];
    fpPts = [predX(~isTP), predY(~isTP)];
    matchedPolys = gtPolygons(matchedGT);
    fnPolys = gtPolygons(~matchedGT);
    
    P = TP / max(1, (TP + FP)); R = TP / max(1, (TP + FN)); F1 = 2 * (P * R) / max(eps, (P + R));
    Metrics = [P, R, F1];
end

function disegnaMappa(ax, imgSfondo, offsetX, offsetY, matchedPolys, fnPolys, tpPts, fpPts, titolo)
    axes(ax); imshow(imgSfondo); hold on;
    title(titolo, 'FontSize', 14);
    
    % Box Gialli (Veri Positivi)
    for i = 1:length(matchedPolys)
        vX = matchedPolys(i).Vertices(:,1) - offsetX + 1;
        vY = matchedPolys(i).Vertices(:,2) - offsetY + 1;
        plot(polyshape(vX, vY), 'FaceColor', 'none', 'EdgeColor', 'y', 'LineWidth', 1.5);
    end
    
    % Box Rossi (Falsi Negativi)
    for i = 1:length(fnPolys)
        vX = fnPolys(i).Vertices(:,1) - offsetX + 1;
        vY = fnPolys(i).Vertices(:,2) - offsetY + 1;
        plot(polyshape(vX, vY), 'FaceColor', 'none', 'EdgeColor', 'r', 'LineWidth', 1.5);
    end
    
    % Pallini Verdi (TP) e Magenta (FP)
    if ~isempty(tpPts)
        plot(tpPts(:,1) - offsetX + 1, tpPts(:,2) - offsetY + 1, '.g', 'MarkerSize', 15);
    end
    if ~isempty(fpPts)
        plot(fpPts(:,1) - offsetX + 1, fpPts(:,2) - offsetY + 1, '.m', 'MarkerSize', 15);
    end
end

