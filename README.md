# FluidDynamicTools
This repository contains most of the post processing scripts I wrote in my work related to CFD, DSMC, MD and other numerical simulation methods of fluid dynamics. I hope it could make the life of my colleagues easier. Most of the scripts are written in Matlab and python. Feel free to use them and distribute them FREELY. But I don't take any responsibility for the scripts.

## Plot3D_Tools
These are some matlab scripts to preprocess plot3d format data. They were developed when I was trying to postprocess results calculated from [Eilmer 3](http://cfcfd.mechmining.uq.edu.au/eilmer3.html). So use it carefully and you may need to change some codes. A good reference of plot3d file is [Plot3D Manual](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19900013774.pdf). The scripts include the following ones:
 - **ReadCellGrid.m:** Read Eilmer3 generated cell-centered grid. (You may not need it in your work)
 - **ReadGrid.m:** Read normal plot3d vertex-centered multiblock grid file (ext: g, xyz)
 - **ReadFlowEilmer3.m:** Read plot3d multiblock function file from Eilmer3 (ext: f, fmt)
 - **ReadFlow.m:** Read plot3d multiblock function file
 - **ReadName.m:** Read name file (ext: nam, name)
 - **ConcateBlocks.m:** Merge multiblock grid file and function file (Only test for my 2D cases which have uniformed divided multi-block grid)
 - **ExportXYZ.m/ExportFUNCTION.m:** Export correct plot3d format
 - **PLOT3D_post.m:** script to involke the post processing of previous sample
 - _sphere0.**_: Example input of PLOT3D_post script
 - _sphere.**_: Desired output from PLOT3D_post script
Use these scripts with care. In addition, I strongly recommend use [mat2tecplot](https://www.cfd-online.com/Forums/tecplot/103860-mat2tecplot.html) to export the data to tecplot.

## Tecplot Tools
 - **mat2tecplot**.m: This is an archive of [mat2tecplo](https://www.cfd-online.com/Forums/tecplot/103860-mat2tecplot.html). Instructions are included in the file. It builds a connection between Matlab to Tecplot by outputing binary tecplot data file.
 - **mat2dat.m**: This is my own script to export data from Matlab to Tecplot. Instead of exporting binary file, it exports ASCII file, which is human readable. Arguments of the function is similar to mat2tecplot. Currently it only supports X-Y line or 2D data.
 - **tec2mat.m**: It serves the opposite function compared to *mat2dat.m*. I usually used it when I need to use Matlab to postprocess results from other software.
 - **tec2py.py**: Python version of *tec2mat.m*.
 - **ExtractStreamLine.py**: Extract staganation streamline starting at I=0 from a 2D contour and export the results to a dat file.
 - **tecplot.mplstyle**: Fake tecplot style in matplotlib of python. Used by call `matplotlib.style.use('https://raw.githubusercontent.com/luohancfd/FluidDynamicTools/master/Tecplot_Tools/tecplot.mplstyle')`
 -
## Publication Tools
 - **ExportEPS.py**: A simple python script to batch export *eps* files from Tecplot layout. To use it, you must install Tecplot (ver > 2017) with pytecplot support and have a valid license for pytecplot.
 - **cleanbib.py**: Clean bibitem for AIAA conference. It adds "AIAA Paper XXX" to the "series" field. To use it, you must have valid doi for the bibitems.

## Eilmer3
These are some modifications I made for [Eilmer 3](http://cfcfd.mechmining.uq.edu.au/eilmer3.html).
 - **VT relaxation model**: A new model called PolyFit, which can be used to for polynominal fit VT relaxation time or Eq. 5 in [paper](http://aip.scitation.org/doi/abs/10.1063/1.4813070). Check *energy-exchange-relaxation-time.cxx* for detail. A example input *o2-VT-QCT.lua* is included.
 - **poshax2tec.py**: Convert poshax3 output to tecplot format
 - **compileXX**: one-click script to compile the code
Eilmer3 is licsensed with GPL3.0. I greatly appreciate them to make a compressible CFD solver open source.

## Thermo_Chemical_Properties
 - Scripts to load thermochemical database and calculate properties
 -
# Licenses:
All the programs are licensed with [GPL 3.0](https://www.gnu.org/licenses/gpl-3.0.html). You are welcome to distribute the codes and you are not required to cite anything. But if you make any modifications, you are required to open source. I also hope you can submit a pull request if possible.
