
% The dvFileInputs function lets the user select a Deltavision 
% multi-channel .dv image file, split the original file into the individual
% image planes, visualize them, and save them as separate .ome.tiff files.

% The indexing of the new .ome.tiff files follows the indexing of the .dv
% file - when imaging the sample channel 1 = bright field, 2 = FITC, 3 =
% TRITC.

% dvfile contains both image and metadata information from the Deltavision
% file.  Used bfopen function to extract data, image data assigned to R,
% where image data is then split out into R1.

function [dvfile,R1,n,R] = dvFileInputs;
% Matlab's imread function does not work with .dv files.  Therefore
% importing the bio-formats toolbox was required - see documentation here:
% https://docs.openmicroscopy.org/bio-formats/6.2.1/users/matlab/
% https://docs.openmicroscopy.org/bio-formats/6.2.1/developers/matlab-dev.html
% https://github.com/mcianfrocco/Matlab/blob/master/u-track-2.1.1/software/bioformats/bfopen.m

[dvfile, dvpath] = uigetfile('*.*');
R = bfopen(dvfile);   % bfopen is a function from the bio-formats toolbox

%%
% Bfopen worked great to open the .dv file.  In order to open the .dv 
% file, we need to unwrap the image planes.  The documentation above 
% details how this works.

RCount = size(R,1);
n = size(R{1,1},1);
R1 = R{1,1};
metadataList = R{1,2};

%%
% The variable n is the index number and represents the number of image
% planes in the .dv file.  The var dvfile is a string which is assigned the
% file name via the uigetfile function.  I subtract 3 chars to remove the
% .dv extension, and then add the index number in its place.

for n=1:n;
    R1plane = R1{n,1};
    R1label = R1{n,2};
    figure('Name',R1label(:,:));
    imshow(R1plane,[]);
    pause(1);
    dvfilename = [dvfile(1:end - 3) '_' num2str(n,'%d') '.ome.tiff'];
    bfsave(R1plane,dvfilename);
end

clearvars -except dvfile R1 n R


end

