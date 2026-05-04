% -- ALGORITMO DI SEGMENTAZIONE DELLE CHIOME (basato sul paper di Erikson 2003) --
% Crop Interattivo + Semi sui Picchi + Integrazione CHM

% pulisco l'ambiente di lavoro
clear; clc; close all;

fprintf('=== Tree Crown Segmentation ===\n');
tTotal = tic;

%% 1. CARICAMENTO E PRE-PROCESSING (CON CROP INTERATTIVO)
% leggo l'immagine rgb
img = imread('ABBY_rgb/2019_ABBY_3_555000_5067000_image.tif');
[origH, origW, bands] = size(img);
% converto in double
imgDouble = double(img);

% Leggo il file chm (Canopy Height Model)
chmFull = imread('ABBY_chm/NEON_D16_ABBY_DP3_555000_5067000_CHM.tif');
% Rimuovo i valori 'NoData' (spesso -9999 nei dati NEON) e le anomalie:
% porto tutti i valori negativi a 0 per farli corrispondere al livello del suolo.
chmFull(chmFull < 0) = 0; 
% converto in double
chmDouble = double(chmFull);

% --- ALLINEAMENTO ---
% Ridimensiono l'INTERO CHM per farlo combaciare con i pixel dell'immagine
% ottica, prima del crop
chmDouble_Resized = imresize(chmDouble, [origH, origW], 'nearest');

% --- NORMALIZZAZIONE GLOBALE ---
% sommo i valori dei colori dell'immagine.
% se la somma è < 0, allora non è una sezione completamente nera (come i
% bordi). Creo una mappa in bianco e nero dove il bianco è la foresta e il
% nero è il bordo vuoto.
validMask = sum(imgDouble, 3) > 0; 
if sum(validMask(:)) == 0
    error('L''immagine è completamente nera!');
end
% estraggo solo i pixel validi, ignorando i bordi neri. Cerco il valore min
% e max di luminosità nella scena utilizzando i percentili e non i valori
% assoluti per via del rumore del sensore: taglio via l'1% dei pixel
% estremi e mi concentro sul 98% dei dati reali della foresta.
validPixels = imgDouble(repmat(validMask, [1 1 bands]));
% trova la soglia sotto cui si trova solo l'1% dei pixel più scuri
low_in = prctile(validPixels, 1);   
% trova la soglia sopra cui si trova solo l'1% dei pixel più luminosi
high_in = prctile(validPixels, 99); 
% rendo tutti i valori dell'immagine compresi tra 0.0 e 1.0
for b = 1:bands
    band = imgDouble(:,:,b);
    % imposta tutti i pixel minori della soglia minima al valore minimo
    band(band < low_in) = low_in;
    % imposta tutti i pixel maggiori della soglia massima al valore massimo
    band(band > high_in) = high_in;
    % normalizza tutto tra 0.0 e 1.0
    imgDouble(:,:,b) = (band - low_in) / (high_in - low_in);
end

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

% Visualizzazione diagnostica ICELab
figure('Name', 'Analisi Componente ICELab', 'Position', [100 100 1200 400]);

subplot(1,3,1);
imshow(imgDouble); 
title('Originale (RGB)');

subplot(1,3,2);
imshow(comp_a_filtered_norm, []); 
colormap(gca, 'parula'); % Usa una mappa di calore per vedere meglio i gradienti
title('Componente a (Normalizzata)');
colorbar;

subplot(1,3,3);
imshow(imgForGrowing, []); 
title('Banda Ibrida (Growing Input)');

% --- CROP INTERATTIVO ---
% variabile booleana per decidere se utilizzare un subset dell'immagine
% oppure l'immagine intera
useSubset = true; 
if useSubset
    fprintf('Seleziona un''area di foresta densa disegnando un rettangolo sull''immagine...\n');
    
    figCrop = figure('Name', 'Seleziona l''area di test', 'WindowState', 'maximized');
    imshow(imgDouble);
    title('traccia un rettangolo su un''area di foresta');
    
    rect = getrect(figCrop); 
    close(figCrop);
    
    % Applico il ritaglio a entrambe le matrici (rgb e chm) usando lo stesso rettangolo
    img = imcrop(img, rect);       
    imgDouble = imcrop(imgDouble, rect); 
    chmCrop = imcrop(chmDouble_Resized, rect); 
    %Ritaglio anche la banda ibrida
    imgForGrowing = imcrop(imgForGrowing, rect);

    comp_a_crop = imcrop(comp_a_filtered_norm, rect);

    figure('Name', 'Diagnostica Spettrale Subset', 'Position', [100 100 1300 450]);
    
    subplot(1,3,1);
    imshow(imgDouble); 
    title('Subset Originale (RGB)');
    
    subplot(1,3,2);
    % Uso jet o parula per evidenziare i gradienti biochimici della componente 'a'
    imshow(comp_a_crop, []); 
    colormap(gca, 'jet'); 
    title('Componente a (ICELab)');
    colorbar;
    
    subplot(1,3,3);
    % Questa è l'immagine effettiva usata per la crescita della regione
    imshow(imgForGrowing, []); 
    title('Banda Ibrida (Input Algoritmo)');
    
    drawnow; % Forza il rendering della figura prima di procedere
    fprintf('  Visualizzazione diagnostica generata. Subset: %d x %d pixel.\n', size(img,1), size(img,2));
    
    fprintf('  Crop applicato: %d x %d pixel.\n', size(img,1), size(img,2));
else
    fprintf('Elaborazione dell''intera immagine in corso...\n');
    chmCrop = chmDouble_Resized;
    cMin = 1; 
    rMin = 1; 
    % utilizzo la larghezza e la l'altezza originali dell'immagine
    cMax = origW;  
    rMax = origH;  
end

% Aggiorno le dimensioni dell'immagine dopo l'eventuale ritaglio
[rows, cols, ~] = size(imgDouble);
% Applico una leggera sfocatura per aiutare i seed a cadere al centro della chioma
% (fonde leggermente insieme i micro-picchi di luce e le micro-ombre
imgDouble = imgaussfilt(imgDouble, 0.5);

%% 2. PARAMETRI
% Soglia di arresto: l'espansione si ferma se la similarità scende sotto il 5%
ALPHA = 0.10;       
% Deviazione standard radiometrica: modella la tolleranza alle variazioni di luminosità (ombre/foglie) nella chioma
% SIGMA1 = 0.36; 
SIGMA1 = 0.22;
% Filtro morfologico (in pixel): scarta i cluster troppo piccoli (es. < 0.3 mq) considerandoli rumore o arbusti
MIN_AREA = 30;      
% Limite spaziale (in pixel): raggio massimo di espansione (es. 4 metri). Ottimizza i tempi e impone limiti biologici
MAX_BOX_RADIUS = 40; 

%% 3. INDIVIDUAZIONE SEED POINTS (PICCHI LUMINOSI)
fprintf('--- Finding Seed Points ---\n');

% Maschera vegetazione (Otsu)
% seleziono la banda del NIR
% nir = imgDouble(:,:,1);
target_per_seeds = imgForGrowing;
level = graythresh(target_per_seeds);
% Rilasso la soglia del 10% per includere i bordi meno luminosi della chioma
nirThresh = level * 0.9;
binaryTree = target_per_seeds > nirThresh;
% Apertura (erosione + dilatazione): rimuove micro-rumore isolato
binaryTree = imopen(binaryTree, strel('disk', 1));


% MASCHERA DI ALTEZZA DAL CHM
% considero solamente le zone che hanno un'altitudine maggiore di 2 metri
ALTEZZA_MINIMA = 2.0; 
maskAltezza = chmCrop >= ALTEZZA_MINIMA;

% FUSIONE DELLE MASCHERE: Deve essere verde e alto almeno 2 metri
binaryTree = binaryTree & maskAltezza; 

% Applico un filtro gaussiano forte, al fine di estrarre un unico picco max
% per ogni albero
sigma_seeds = 1.8; 
hybrid_smoothed = imgaussfilt(target_per_seeds, sigma_seeds);

% Trova i picchi (cime degli alberi) utilizzando la trasformata H-Maxima
% con soglia 0.06
seedsBinary = imextendedmax(hybrid_smoothed, 0.06); 
% combina con la maschera precedente 
seedsBinary = seedsBinary & binaryTree; 

% estraggo le coordinate (riga e colonna) dei seed validi
[seedRows, seedCols] = find(seedsBinary);
if ~isempty(seedRows)
    % converto le coord 2d in indici lineari 1d
    seedVals = sub2ind(size(target_per_seeds), seedRows, seedCols);
    seedVals = target_per_seeds(seedVals); 
    % ordino in maniera decrescente per dare priorità agli alberi più
    % chiari e alti
    [~, sortIdx] = sort(seedVals, 'descend');
    seedList = [seedRows(sortIdx), seedCols(sortIdx)];
else
    seedList = [];
    warning('Nessun seed trovato. Prova ad abbassare il valore in imextendedmax.');
end
fprintf('  Seeds trovati: %d\n', length(seedList));

%FILTRO SPAZIALE: DISTANZA MINIMA TRA I SEMI
% imposto la distanza a 20 pixel (2 metri reali)
MIN_SEED_DIST = 20; 
filteredSeedList = [];
fprintf('  Applicazione filtro spaziale (Distanza minima: %d px)...\n', MIN_SEED_DIST);

for i = 1:size(seedList, 1)
    pt = seedList(i, :);
    if isempty(filteredSeedList)
        % il primo seme entra nella lista
        filteredSeedList = [filteredSeedList; pt];
    else
        % calcola la distanza euclidea tra il candidato e tutti i semi già
        % inseriti nella lista
        dists = sqrt((filteredSeedList(:,1) - pt(1)).^2 + (filteredSeedList(:,2) - pt(2)).^2);
        % se la distanza minima tra il candidato e i semi già salvati è maggiore o uguale a 20, lo salva.
        % Altrimenti, lo scarta
        if min(dists) >= MIN_SEED_DIST
            filteredSeedList = [filteredSeedList; pt];
        end
    end
end
seedList = filteredSeedList;
nSeeds = size(seedList, 1);
fprintf('  Seeds rimasti dopo il filtro spaziale: %d\n', nSeeds);

%% 4. LOOP PRINCIPALE 
% mappa finale (ogni albero avrà un id numerico)
finalLabelMap = zeros(rows, cols);
% contatore degli alberi validati
currentLabel = 0;
% maschera globale per evitare che gli alberi si sovrappongano
processedMask = false(rows, cols); 

hWait = waitbar(0, 'Segmentazione in corso...');
nSeeds = size(seedList, 1);

for k = 1:nSeeds
    if mod(k, 20) == 0
        waitbar(k/nSeeds, hWait, sprintf('Albero %d / %d', k, nSeeds));
    end
    
    startPt = seedList(k, :);
    
    % se il seed corrente è stato già inglobato dall'espansione di un
    % albero precedente, lo salto
    if processedMask(startPt(1), startPt(2))
        continue;
    end
    % fase di esplorazione: sigma2 fisso per capire le dimensioni
    sigma2_init = 4; 
    regionInit = growRegionFast(imgForGrowing, startPt, SIGMA1, sigma2_init, ALPHA, processedMask, MAX_BOX_RADIUS);
    % ignoro se è un falso positivo molto piccolo
    if sum(regionInit(:)) < 10
        continue; 
    end
    
    % calcolo dinamico di sigma2 (adatto il vincolo spaziale alle
    % dimensioni reali appena calcolate)
    props = regionprops(regionInit, 'EquivDiameter', 'Centroid');
    if isempty(props), continue; end
    % raggio stimato
    rho = props(1).EquivDiameter / 2;
    
    % Formula del paper di Erikson: sigma2 = 2*rho / sqrt(-2*ln(alpha))
    denom = sqrt(-2 * log(ALPHA));
    sigma2_est = (2 * rho) / denom;
    % Impongo un tetto massimo di sicurezza
    sigma2_est = min(sigma2_est * 1.3, 40); 
    
    % ottimizzazione del seed: uso il baricentro per cercare il centro più
    % luminoso
    cent = props(1).Centroid; 
    searchRad = max(2, round(rho * 0.5));
    
    [cGrid, rGrid] = meshgrid( -searchRad:2:searchRad, -searchRad:2:searchRad );
    candC = round(cent(1)) + cGrid(:);
    candR = round(cent(2)) + rGrid(:);
    
    valid = candR>0 & candR<=rows & candC>0 & candC<=cols;
    candR = candR(valid);
    candC = candC(valid);

    % Monte Carlo sampling: testo al max 15 candidati per velocizzare
    MAX_ATTEMPTS = 15;
    if length(candR) > MAX_ATTEMPTS
        idx = randperm(length(candR), MAX_ATTEMPTS);
        candR = candR(idx);
        candC = candC(idx);
    end
    
    bestMetric = Inf;
    bestRegion =[];
    foundCandidate = false;
    
    % espansione finale
    for i = 1:length(candR)
        candPt = [candR(i), candC(i)];
        
        if processedMask(candPt(1), candPt(2))
            continue;
        end
        
        % espando la regione candidata usando il sigma2 dinamico
        candRegion = growRegionFast(imgForGrowing, candPt, SIGMA1, sigma2_est, ALPHA, processedMask, MAX_BOX_RADIUS);
        
        area = sum(candRegion(:));
        if area < MIN_AREA
            continue; 
        end
        
        % check validità della regione: deve essere composta per l'85% da
        % pixel verdi e con altezza < 2 m
        intersectTree = candRegion & binaryTree;
        ratioTree = sum(intersectTree(:)) / area;
        if ratioTree < 0.85 
            continue; 
        end
        
        % cerco il contorno più scuro
        perim = bwperim(candRegion);
        meanContourVal = mean(target_per_seeds(perim));
        
        if meanContourVal < bestMetric
            bestMetric = meanContourVal;
            bestRegion = candRegion;
            foundCandidate = true;
        end
    end
    
    % salvataggio dell'albero
    if foundCandidate
        currentLabel = currentLabel + 1;
        finalLabelMap(bestRegion) = currentLabel;
        % setto a true nella maschera la sezione appena processata per evitare
        % sovrapposizioni successive
        processedMask(bestRegion) = true;
    else
        processedMask(startPt(1), startPt(2)) = true;
    end
end
close(hWait);
fprintf('  Segmentazione completata. Alberi trovati: %d\n', currentLabel);
fprintf('  Tempo totale: %.2f s\n', toc(tTotal));

%% 5. VISUALIZZAZIONE
figure('Name', 'Segmentation Results', 'Position',[100 100 1200 500]);
subplot(1,3,1);
imshow(img); title('Original Image');
imwrite(img, 'og_image_icelab_filtered_dense.png');

subplot(1,3,2);
rgb = label2rgb(finalLabelMap, 'jet', 'k', 'shuffle');
imshow(rgb); title(['Segmentazione (N=', num2str(currentLabel), ')']);
imwrite(finalLabelMap, 'segmentation_icelab_filtered_dense.png');

subplot(1,3,3);
base_contours = bwperim(finalLabelMap > 0);
% applico la dilatazione per aumentare lo spessore dei contorni
thickness = strel('disk', 1); 
thick_contours = imdilate(base_contours, thickness);
color = [1 0 0];
B = imoverlay_custom(img, thick_contours, color);
imshow(B); title('Contorni');
imwrite(B, 'contorni_icelab_filtered_dense.png');

%% --- PREPARAZIONE DATI PER LA VALIDAZIONE ---
fprintf('Preparazione dei dati per la validazione...\n');

% Se ho fatto il crop, utilizzo rect. Altrimenti, prendo l'intera immagine.
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

imgPath = 'ABBY_rgb/2019_ABBY_3_555000_5067000_image.tif';
shpPath = 'ABBY_labels/2019_ABBY_3_555000_5067000_image.shp';
[~, R] = readgeoraster(imgPath);
truthShapes = shaperead(shpPath);

save('dati_validazione.mat', 'finalLabelMap', 'cMin', 'rMin', 'cMax', 'rMax', 'R', 'truthShapes', 'chmCrop');
fprintf('Dati salvati con successo in dati_validazione.mat!\n');

%  FUNZIONI HELPER
% -- IMPLEMENTAZIONE DELL'ALGORITMO DI "FUZZY REGION GROWING" -- 
% calcola la probabilità che un pixel appartenga alla chioma in base a 
% colore e distanza
function mask = growRegionFast(img, seed, s1, s2, alpha, globalMask, maxRad)
    [H, W, channels] = size(img);
    r0 = seed(1); c0 = seed(2);
    
    % ottimizzazione: non lavora sull'immagine intera ma 
    % ritaglia un quadrato di 40 px attorno al seed 
    rMin = max(1, r0 - maxRad); rMax = min(H, r0 + maxRad);
    cMin = max(1, c0 - maxRad); cMax = min(W, c0 + maxRad);
    
    imgCrop = img(rMin:rMax, cMin:cMax, :);
    
    lr0 = r0 - rMin + 1;
    lc0 = c0 - cMin + 1;
    
    [colsGrid, rowsGrid] = meshgrid(1:size(imgCrop,2), 1:size(imgCrop,1));
    % vincolo spaziale: più ci allontaniamo dal seed, più il valore di
    % probabilità spaziale mu2 va verso lo zero, seguendo una curva
    % gaussiana controllata da sigma2
    distSq = (rowsGrid - lr0).^2 + (colsGrid - lc0).^2;
    mu2 = exp(-0.5 * distSq / (s2^2));
    
    seedColor = img(r0, c0, :);
    % vincolo radiometrico: prendo il colore del seed e calcolo quanto ogni
    % altro pixel nel crop è diverso da esso. Se un pixel è molto più
    % scuro, la sua probabilità radiometrica va verso lo zero, controllata
    % da sigma1
    diffCol = imgCrop - seedColor;
    if channels > 1
        colDistSq = sum(diffCol.^2, 3);
    else
        colDistSq = diffCol.^2;
    end

    mu1 = exp(-0.5 * colDistSq / (s1^2));
    
    % la probabilità totale di appartenenza è il prodotto delle due
    % probabilità
    muTotal = mu1 .* mu2;
    % se il prodotto supera la soglia minima alpha, viene annesso alla
    % regione
    localMask = muTotal > alpha;
    % mantengo nella regione solo i pixel adiacenti al seed centrale
    % (connettività a 8)
    localMask = bwselect(localMask, lc0, lr0, 8);
    % riempio i buchi neri all'interno della maschera
    localMask = imfill(localMask, 'holes');
    
    globalCrop = globalMask(rMin:rMax, cMin:cMax);
    % sottraggo i pixel che appartengono già ad altri alberi 
    localMask = localMask & ~globalCrop;
    mask = false(H, W);
    % inserisco il quadrato considerato nella mappa dell'immagine originale
    mask(rMin:rMax, cMin:cMax) = localMask;
end
% -- FUNZIONE DI VISUALIZZAZIONE --
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