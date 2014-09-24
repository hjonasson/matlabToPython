'''
This is a translation of MarginsViaPointPaths.m into Python
'''
import numpy as np
import shapefile

desamp 						= 	2
unstretchedCrustThickness	=	38
defZoneShapeFile 			= 	'DefZones/SWPAC_DefZones_20140220.shp'
rotFile 					=	'../RotationFiles/Wright_901_20140220.rot'
crustThickFile 				=	'DefZones/crust1_CrustalThickness_Padded.nc'

movPlate 					=	901
fixPlate 					= 	804
outFileStem 				=	'PACMBL-PointPathTest-20140220'
outDir						=	'DefZones/'

finitePoleTimes				=	[73.600, 79.000, 90.000]

fixBasinID 					=	89000 + fixPlate
movBasinID 					=	89000 + movPlate

t0 							=	89 #startTime = 89;
tn 							=	74 #endTime = 74;
dT							=	1
startOfRifting 				=	89
endOfRifting 				=	74


fixPlateIDRoot 				=	1000
movPlateIDRoot 				=	2000

### Translated above, to be done below ###
'''
[CatCOBF,RestoredPointsF] = MarginsViaPointPaths_function(Params);

Params.MovPlate = FIXPLATE;
Params.FixPlate = MOVPLATE;
Params.FixPlateIDRoot = 2000;
Params.MovPlateIDRoot = 1000;

[CatCOBM,RestoredPointsM] = MarginsViaPointPaths_function(Params);
'''
### Translate MarginsViaPointPaths_function.m ###
# to do
# getMargin





def marginsViaPointPathsFunction(defZoneShapeFile,fixPlate,movPlate,desamp):

	# Get margin geometries
	fixMargin 				=	getMargin(defZoneShapeFile,fixPlate,movPlate,desamp)
	MovMargin 				=	getMargin(defZoneShapeFile,movPlate,fixPlate,desamp)

	# Assign unique IDs to each point



	return catCOB,restoredPoints

def getMargin(shapeFile,plate1,plate2,desamp):
	thisMargin = True
	return thisMargin

def mexShape(shapeFile,*args):
	
	filterFields			=	args[0::2]
	filterValues			=	args[1::2]
	sf 						= 	shapefile.Reader(shapeFile)
	sAttributes 			= 	sf.records()
	sFieldNames 			= 	sf.fields
	selectionIndex 			=	np.zeros(len(sAttributes))

	for i in range(len(filterFields)):

		filterFieldIndex 	= [j-1 for j in range(len(sFieldNames)) if sFieldNames[j][0] == filterFields[i]]
		filterFieldList 	= [j[filterFieldIndex[0]] for j in sf.iterRecords()]
		for k in range(len(filterFieldList)):
			if filterFieldList[k] == filterValues[i]:
				selectionIndex[k] = selectionIndex[k] + 1

	sIndex 					= [l for l in range(len(selectionIndex)) if selectionIndex[l] == len(filterFields)]
	data = [sf.shape(m) for m in sIndex]
	return data





### To be translated ###

'''
Get margin geometries, assign unique IDs to each point
FixMargin.Seedpoints.ID = FixMargin.Seedpoints.ID+Params.FixPlateIDRoot;
MovMargin.Seedpoints.ID = MovMargin.Seedpoints.ID+Params.MovPlateIDRoot;

Create point paths from each point on the moving margin, for the time
range within which we think the overlaps will occur
ReconTimes = [Params.EndTime:Params.deltaT:Params.StartTime];
%[FinitePoles] = GetPolesUsingRotfnd(Params.RotFile,Params.FixPlate,Params.MovPlate,ReconTimes)
[FinitePoles] = GPlatesGetRotationPole(ReconTimes,Params.FixPlate,Params.MovPlate,Params.RotFile)
for i=1:size(FinitePoles,1)
    [PointPaths.X(:,i),PointPaths.Y(:,i)] = ...
        EulerRotation(FinitePoles(i,4),FinitePoles(i,3),FinitePoles(i,5),MovMargin.Seedpoints.X,MovMargin.Seedpoints.Y);
    PointPaths.T(i) = FinitePoles(i,2);
end

figure
plot(MovMargin.Seedpoints.X,MovMargin.Seedpoints.Y,'-ro')
hold on
plot(FixMargin.Seedpoints.X,FixMargin.Seedpoints.Y,'-ro')
plot(PointPaths.X',PointPaths.Y','.b-')
axis equal xy


Fix Margin COB --> 3D cartesian 
[FixSdptxyz] = llp2uv3d(FixMargin.Seedpoints.X,FixMargin.Seedpoints.Y);

an array of distance along line
Xi = [0;cumsum(LL2Distance(FixMargin.Seedpoints.X,FixMargin.Seedpoints.Y))];

Loop over each point path
for i=1:size(PointPaths.X,1);
    
This point path to 3d cartesian
    [PPxyz] = llp2uv3d(PointPaths.X(i,:),PointPaths.Y(i,:));
    
    try
Intersection between this point path and COB
        [FixSdptIntxyz,ind.PPxFixSdpt,ind.FixSdptxPP] = LineIntersectionOnSphere(PPxyz,FixSdptxyz);
        
Intersection Point back to Long/Lat
        [FixSdptInt.Lon(i),FixSdptInt.Lat(i)] = uv3d2llp(FixSdptIntxyz);
        
find the distance along line of the intersection, so we can correctly
order all the points at a later stage
        intX(i) = interp1(FixMargin.Seedpoints.X(ind.FixSdptxPP:ind.FixSdptxPP+1), ...
            Xi(ind.FixSdptxPP:ind.FixSdptxPP+1), ...
            FixSdptInt.Lon(i),'linear');
        
NEED TO ADD SOME LINES TO DETERMINE THE TIME OR BREAKUP FOR EACH POINT
BASED ON WHERE IT INTERSECTS THE POINTPATH LINE
        intTIME(i) = interp1(PointPaths.X(i,ind.PPxFixSdpt:ind.PPxFixSdpt+1), ...
            ReconTimes(ind.PPxFixSdpt:ind.PPxFixSdpt+1), ...
            FixSdptInt.Lon(i),'linear');
        
        
    catch
        disp('No intersection found')
        FixSdptInt.Lon(i) = nan;
        FixSdptInt.Lat(i) = nan;
        intX(i) = nan;
        intTIME(i) = nan;
    end
    
end

plot(FixSdptInt.Lon,FixSdptInt.Lat,'bs')

Use the along-line distance of each point
CatX = [Xi;intX'];
[Y,I] = sort(CatX);
CatLon = [FixMargin.Seedpoints.X;FixSdptInt.Lon'];
CatLat = [FixMargin.Seedpoints.Y;FixSdptInt.Lat'];
CatID = [FixMargin.Seedpoints.ID;MovMargin.Seedpoints.ID];
CatTIME = [zeros(size(FixMargin.Seedpoints.ID));intTIME'];
CatCOB.X = CatLon(I);
CatCOB.Y = CatLat(I);
CatCOB.ID = CatID(I);
CatCOB.OverlapAge = CatTIME(I);

ind = ~isnan(CatCOB.X);
CatCOB.X = CatCOB.X(ind);
CatCOB.Y = CatCOB.Y(ind);
CatCOB.ID = CatCOB.ID(ind);
CatCOB.OverlapAge = CatCOB.OverlapAge(ind);

Copied from GridTestBigLoop
[FinitePoles] = GetPolesUsingRotfnd(Params.RotFile,Params.FixPlate,Params.MovPlate,[Params.EndOfRifting Params.StartOfRifting]);
ReconTimes = [Params.EndOfRifting:Params.deltaT:Params.StartOfRifting];
[FinitePoles] = GPlatesGetRotationPole(ReconTimes,Params.FixPlate,Params.MovPlate,Params.RotFile)
Poles.End.Long = FinitePoles(1,4);
Poles.End.Lat = FinitePoles(1,3);
Poles.End.Angle = FinitePoles(1,5);
Poles.Start.Long = FinitePoles(end,4);
Poles.Start.Lat = FinitePoles(end,3);
Poles.Start.Angle = FinitePoles(end,5);

[FLArray] = MakeSmallCircles(CatCOB,FixMargin.BoundingPolygon,Poles);


[Grid] = LoadCrustalThicknessGrid(Params.CrustThickFile,Params.UnstretchedCrustThickness);


Restore Points based on the calculated small circles and the crustal
thickness
[RestoredPoints] = RestoreSmallCircles(FLArray,Grid,Params.UnstretchedCrustThickness);
[FinitePoles] = GPlatesGetRotationPole(ReconTimes,Params.FixPlate,Params.MovPlate,Params.RotFile)



Added 2014/01/11
Get instantaneous velocity at 
[dX1,dY1] = EulerRotation(FinitePoles(end,4),FinitePoles(end,3),FinitePoles(end,5),CatCOB.X,CatCOB.Y);
[dX2,dY2] = EulerRotation(FinitePoles(end-1,4),FinitePoles(end-1,3),FinitePoles(end-1,5),CatCOB.X,CatCOB.Y);
Vel = LL2Distance([dX1,dX2],[dY1,dY2]);

CatCOB.InitVel = Vel;
keyboard
'''