# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 16:13:40 2018

@author: Han
"""
import logging
logging.basicConfig(level=logging.INFO)
import tecplot as tp
from tecplot.constant import *
import os,sys


def ExportStreamLine (workdir,file1,zonename):
    # workdir = "E:\\Research\\CFD\\RAMCII\\GridConvergence\\91x100"
    datafile = os.path.join(workdir, file1)
    file2 = 'streamline.dat'
    dataset = tp.data.load_tecplot(datafile)
    frame = tp.active_frame()
    frame.plot_type = tp.constant.PlotType.Cartesian2D
    
    # Setup up vectors and background contour
    plot = frame.plot()
    plot.vector.u_variable = dataset.variable('u')
    plot.vector.v_variable = dataset.variable('v')
    plot.contour(0).variable = dataset.variable('t')
    plot.contour(1).variable = dataset.variable('tv')
    plot.contour(2).variable = dataset.variable('p')
    
    # get starting point of stagnation streamline
    x = dataset.zone(0).values('x')
    y = dataset.zone(0).values('y')
    
    x0 = x[0]
    y0 = y[0] 
    
    plot.show_streamtraces = True
    plot.show_contour = True
    plot.fieldmap(1).contour.show = True
    
    # extract the data from streamline
    streamtraces = plot.streamtraces
    streamtraces.add(seed_point=[x0, y0],stream_type=Streamtrace.TwoDLine)
    tp.macro.execute_command('$!ExtractStreamtraces')
    dataset.zone(1).name = zonename   
      
    tp.data.save_tecplot_ascii(os.path.join(workdir,file2),
        zones=[1],
        include_text=False,
        precision=9,
        include_geom=False,
        use_point_format=True)

#tp.export.save_png('streamtrace_2D.png', 600, supersample=3)
#
#tecplot.save_layout('test.lay')
## Add streamtraces and set streamtrace style
#streamtraces = plot.streamtraces

if __name__ == '__main__':
    ExportStreamLine (sys.argv[1],sys.argv[2],sys.argv[3])