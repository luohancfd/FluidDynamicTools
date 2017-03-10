Fname = 'sphere0';
% Read Eilmer3 exported plot3D data
[Name]  = ReadName([FName,'.nam']);
[ Flow,~,~ ] = ReadFlowEilmer3( [FName,'.f'],0 );
[Grid,~] = ReadCellGrid( [FName,'.g'],0);
[PGrid,~, ~] = ReadGrid([FName,'.grd']);

% concate blocks
map = reshape(1:40,8,[])';
[ NewGrid,NewFlow ] = ConcateBlocks( map,Grid,'Flow',Flow,'Skip',1 );
NewPGrid = ConcateBlocks( map,PGrid,'Skip',0 );

% expord correct plot3d format data
[nblock] = ExportXYZ( NewGrid,'sphere.g');   %cell centered grid
[nblock,Nvar] = ExportFUNCTION( NewFlow,'sphere.f' );
[nblock] = ExportXYZ( NewPGrid,'sphere.grd'); %vertex-centered grid
copyfile([Fname,'.nam'],'sphere.nam');