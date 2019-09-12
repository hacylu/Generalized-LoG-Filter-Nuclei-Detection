addpath(genpath([pwd '\GeneralLoG']));

IMData4SymmetricVoting_1={'8913.tif'};
% IMData4SymmetricVoting_1={'IU_007_1.png'};
% IMData4SymmetricVoting_1={'117_H&E_07_7.tif'};
%IMData4SymmetricVoting_1={'117_Feulgen_1_7.tif'};

for i=1:length(IMData4SymmetricVoting_1)
    curIMName=IMData4SymmetricVoting_1{i};
    curIM=imread(curIMName);
    curIMsize=size(curIM);
    [curIM_norm] = normalizeStaining(curIM);
    curIM_normRed=curIM_norm(:,:,1);
    %% using general LoG,
    %%% initial segmentation
    R=curIM_normRed; I=curIM;
    ac=5;    % remove isolated pixels(non-nuclei pixels)
    [Nmask,cs,rs,A3]=XNucleiSeg_GL_Thresholding(R,ac,I,1);      %% thresholding based binarization
%     show(Nmask); LshowBWonIM(Nmask,I);
    %%% gLoG seeds detection
    largeSigma=16;smallSigma=12; % Sigma value
    ns=XNucleiCenD_Clustering(R,Nmask,largeSigma,smallSigma);  %% To detect nuclei clumps
    
    figure(1),imshow(I);
    hold on,plot(ns(2,:),ns(1,:),'y+');
    hold on,plot(cs,rs,'g*');  
end