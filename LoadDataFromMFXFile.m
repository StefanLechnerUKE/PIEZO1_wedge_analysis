% myfile = "240422-174335_minflux.mat";
% 

function [traces]=LoadDataFromMFXfile(myfile)
load(myfile);


%% isolate localization and all other parameter arrays from the iterations

% localizations
itr_loc = itr.loc;
% counts at offset
itr_eco = itr.eco;
% counts at center position
itr_ecc = itr.ecc;
% center frequency ratio
itr_cfr = itr.cfr;
% frequency at offset
itr_efo = itr.efo;
% frequency at center
itr_efc = itr.efc;
% diffChannel center frequency ratio
itr_dcr = itr.dcr;
% status localization process (0=ok)
itr_sta = itr.sta;
% fluorescence background measured during MINFLUX acquisition (in count frequency)
itr_fbg = itr.fbg;
% diameter of patern itr
itr_ext = itr.ext;

%subsequent parameters reflect changes of microscope hardware
% galvanometer pos x
itr_gvx = itr.gvx;
% galvanometer pos y
itr_gvy = itr.gvy;
% electro-optotically deflector pos x
itr_eox = itr.eox;
% electro-optotically deflector pos y
itr_eoy = itr.eoy;
% deformable mirror
itr_dmz = itr.dmz;
% lcy
itr_lcy = itr.lcy;
% lcx
itr_lcx = itr.lcx;
% lcz
itr_lcz = itr.lcz;


% Convert traceID (tid), time (tim), valid traces or not (vld), act 405nm ON or not (act), stickiness (sky = Amount of tries to localize an emitter before terminating the localizaton, set to 2 here) into double for import
tid = double(tid)';
tim = double(tim)';
vld = double(vld)';
act = double(act)';
sky = double(sky)';

%Extract time (last recorded trace, in s)
recordingTime = max(tim);

%Convert ecc and eco param (int32 format) into numbers
eco_num = double(itr_eco);
ecc_num = double (itr_ecc);


%For 3D, 2D and valid traces only. Change 10 to 5 for 2D mfx
% define x,y,z coordinates for each localization and link the trace ID

    if size(itr_loc,2) == 10;
        x = itr_loc(:,10,1);
        y = itr_loc(:,10,2);
        z = itr_loc(:,10,3);
        % x=x-min(x);
        % y=y-min(y);
        z=z-min(z);
        xyz = [x y z];
        traces = [x(:,1),y(:,1),z(:,1),tid(:,1),itr_efo(:,7),itr_cfr(:,7),itr_fbg(:,7),tim(:,1)];

    else size(itr_loc,2) == 5;
        x = itr_loc(:,5,1);
        y = itr_loc(:,5,2);
        z = itr_loc(:,5,3);
        xyz = [x y z];
        traces = [x(:,1),y(:,1),z(:,1),tid(:,1),itr_efo(:,5),itr_cfr(:,5),itr_fbg(:,5),tim(:,1)];
    end


%Add remaining parameters to the table, depending for 3D or 2D
    if size(itr_loc,2) == 10;
        myIdx=10;
        myIdx1=7;
    else size(itr_loc,2) == 5;
        myIdx=5;
        myIdx1=5;
    end
    eco = eco_num(:,myIdx1,1);
    ecc = ecc_num(:,myIdx1,1);
    cfr = itr_cfr(:,myIdx1,1);
    dcr = itr_dcr(:,myIdx,1);
    ext = itr_ext(:,myIdx,1);
    sta = itr_sta(:,myIdx,1);
    efo = itr_efo(:,myIdx1,1);
    efc = itr_efc(:,myIdx1,1);
    fbg = itr_fbg(:,myIdx,1);
    gvx = itr_gvx(:,myIdx,1);
    gvy = itr_gvy(:,myIdx,1);
    eox = itr_eox(:,myIdx,1);
    eoy = itr_eoy(:,myIdx,1);
    dmz = itr_dmz(:,myIdx,1);
    lcy = itr_lcy(:,myIdx,1);
    lcx = itr_lcx(:,myIdx,1);
    lcz = itr_lcz(:,myIdx,1);

end
