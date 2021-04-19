from os import listdir
from os.path import isfile, join
import numpy as np
import re
from . import disptools
from scipy.interpolate import CubicSpline, interp1d
from argparse import ArgumentParser


# change from metal to Si units
global metal2si, regex

metal2si = 1.60217662e-19*1e20*6.0221409e23*1000
regex = '[+-]?[0-9]' # negative or positive decimal regex

def find_files(wdir:str, pattern:str) -> list:
    """Return ordered list of log file names from a directory

    Use a filename pattern to find the corresponding log files

    :param wdir: working directory
    :type wdir: str
    :param pattern: filename pattern with wildcard
    :type pattern: str
    :return: list of file names in order
    :rtype: list
    """
    mypath = wdir+"/"
    fpattern = re.compile(pattern)
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))] # select only files from listdir
    onlylogfiles = [f for f in onlyfiles if fpattern.match(f)] # select files that only matech pattern
    onlylogfiles.sort(key=lambda f: int(''.join(re.findall(regex, f)))) # sort file names according to regex
    
    return [join(mypath,f) for f in onlylogfiles], get_slider_vals(onlylogfiles),  # return absolute path

def get_slider_vals(files:list) -> list:
    """Return values for the slider from file names

    The file pattern name is used to extract the slider values instead of just using file names

    :param files: file list
    :type files: list
    :return: values for the slider used in plotting
    :rtype: list
    """
    pvals = [''.join(re.findall(regex, f)) for f in files]
    return pvals

def process_log_files(files:list, nintp=10) -> tuple[np.array, list, list, dict]:
    """Return dispersion data for files with interpolation

    :param files: list of file names to process
    :type files: list
    :param SPs: symmetry points
    :type SPs: list
    :return: the x-values, dispersion points, interpolation, and coordinates of special points using cubic splines
    :rtype: tuple[np.array, list, list, list]

    .. todo:: Add custom SP feature
    .. todo:: put qr inside disp_data like disp_inter
    """
    nfiles = len(files)

    # Initializa lists for q,eig for files
    eig_full = []
    q_full = []
    
    # Get dispersion data from each file
    for ifile in range(nfiles):
        eig, q = disptools.get_data(files[ifile])
        eig = np.sqrt(eig*metal2si)/1e12/2/np.pi  # change to THz
        eig_full.append(eig)
        q_full.append(q)

    # Number of eigenvalues for each q
    neigvals = eig_full[0].shape[1]

    # Create data plot
    from .symmetry import FCC
    # qx qy qz qr w1 w2 w3 ....
    SPs = [FCC.SIGMA, FCC.X, FCC.W, FCC.K, FCC.SIGMA, FCC.L, FCC.U, FCC.W, FCC.L]
    SPsNames = ['Γ', 'X', 'W', 'K', 'Γ', 'L', 'U', 'W', 'L']
    disp_data_full = [] # dispersion x,y values
    disp_inter_full = [] # interpolated data

    for ifiles in range(nfiles):
        qr, disp_data, SPqr = disptools.build_curve(
        SPs, eig_full[ifiles], q_full[ifiles])
    
        # interpolated dispersion data
        disp_inter = np.zeros((nintp*(len(SPs)-1), neigvals+1))

        # Interpolate the data
        # iqestart = 0 # unused
        intpx = np.linspace(qr[0], qr[-1], nintp*(len(SPs)-1))
        for ieig in range(neigvals):
            ydat = disp_data[:, 3+ieig]
            f = CubicSpline(qr, ydat, bc_type='clamped')
            disp_inter[:, ieig+1] = f(intpx)

        disp_inter[:, 0] = intpx

        disp_data_full.append(disp_data)
        disp_inter_full.append(disp_inter)
    
    return qr, disp_data_full, disp_inter_full, dict(zip(SPqr, SPsNames))

def make_interactive_plot(qr:np.array, pvals, disp_data_full:list, disp_inter_full:list, labels:dict) -> None:
    """ Create the interactive Bokeh plot as a standalone html web page

    :param qr: x-values for the points
    :type qr: np.array
    :param pvals: slider values
    :type pvals: list
    :param disp_data_full: list of arrays containing q vectors and eigenvalues for points
    :type disp_data_full: list
    :param disp_inter_full: list of arrays containing interpolated x and eigen values
    :type disp_inter_full: list
    :param labels: labels for the x-axis indicating special points
    :type labels: dict
    :return: None

    .. todo:: add output file name option
    """
    
    # Plotting
    from bokeh.plotting import figure,output_file
    from bokeh.io import show, output_notebook, push_notebook
    from bokeh.models import ColumnDataSource, Panel, CustomJS, Circle, Line, Div, HoverTool
    from bokeh.models.widgets import CheckboxGroup, Slider, RangeSlider
    from bokeh.layouts import column, row
    from bokeh.application import Application
    
    nfiles = len(disp_data_full)
    neigvals = disp_inter_full[0].shape[1]-1 # first column is x-values 
    
    # Output file name (web page)
    output_file("disp_app.html")

    # Create the figure
    p = figure(title="Dispersion curves for FCC",
        y_range=(0, 2.25),
        x_range=(0, qr[-1]+0.001),
        x_axis_label="q", y_axis_label="ω (THz)",
        plot_width=800)

    # Glyphs for points and lines (intepolation)
    raw_glyphs = []
    inter_glyphs = []

    # Make dictionaries for Bokeh ColumnDataSource
    draw = {"xraw":qr}
    dint = {"xint":disp_inter_full[0][:,0]}
    for ip in range(nfiles):
        for ie in range(neigvals):
            ei = "e"+str(ip)+"_"+str(ie)
            draw[ei] = disp_data_full[ip][:,ie+3]
            dint[ei] = disp_inter_full[ip][:,ie+1]

    srcraw = ColumnDataSource(draw)
    srcint = ColumnDataSource(dint)
    
    # Build glyphs for each pressure point
    for ip in range(nfiles):
        lalpha = 0
        for ie in range(neigvals):
            ei = "e"+str(ip)+"_"+str(ie)
            if ip==0: lalpha = 1
            raw_glyphs.append(Circle(x="xraw", y=ei, line_alpha=lalpha,fill_alpha=lalpha))
            inter_glyphs.append(Line(x="xint", y=ei, line_alpha=lalpha))
            p.add_glyph(srcraw,raw_glyphs[-1])
            p.add_glyph(srcint,inter_glyphs[-1])


    # Styling
    p.xaxis.ticker = list(labels.keys()) # SPqr
    p.xaxis.major_label_overrides = labels


    # Widgets
    EIGEN_LBLS = ["Eigen "+str(i+1) for i in range(neigvals)]
    eigen_chkbx_grp = CheckboxGroup(labels=EIGEN_LBLS, active=[i for i in range(neigvals)])
    slider = Slider(start=1, end=nfiles, value=1, step=1, title="Pressure series")
    #pvals = [''.join(re.findall(regex, f)) for f in onlylogfiles ]
    div = Div(text="Pressure = "+pvals[0]+" bar", name="pval")

    # Callbacks for widgets
    ## Slider Callback
    slider.js_on_change("value", CustomJS(args=dict(glp=raw_glyphs,glq=inter_glyphs,neig=neigvals,lbl=div, vals=pvals, chkbx=eigen_chkbx_grp), code="""
        var i;
        var eigs = chkbx.active
        for (i=0; i<glp.length; i=i+1){
            if (i>=(this.value-1)*neig && i < this.value*neig){
                if( eigs.includes(i%neig) ){
                    glp[i]["fill_alpha"] = 1;
                    glp[i]["line_alpha"] = 1;
                    glq[i]["line_alpha"] = 1;
                } else {
                    glp[i]["fill_alpha"] = 0;
                    glp[i]["line_alpha"] = 0;
                    glq[i]["line_alpha"] = 0;
                }
                console.log(i,'YES')
            }else{
                glp[i]["fill_alpha"] = 0;
                glp[i]["line_alpha"] = 0;
                glq[i]["line_alpha"] = 0;
                console.log(i,'NO')
            }
        }
        lbl["text"] = "Pressure = " + vals[this.value-1] + " bar"
    """))

    ## Checkbox callback
    eigen_chkbx_grp.js_on_click(CustomJS(args=dict(glp=raw_glyphs,glq=inter_glyphs,neig=neigvals,sldr=slider), code="""
        var actv = this.active
        var i
        console.log("sldr ",sldr.value)
        for (i=(sldr.value-1)*neig; i<sldr.value*neig; i=i+1){
            console.log('i, i%neig res = ', i, i%neig, actv.includes(i%neig))
            if (actv.includes(i%neig)){
                glp[i]["fill_alpha"]=1;
                glp[i]["line_alpha"]=1;
                glq[i]["line_alpha"]=1;
            } else {
                glp[i]["fill_alpha"]=0;
                glp[i]["line_alpha"]=0;
                glq[i]["line_alpha"]=0;
            }
        }
    """))

    ## Add hover widget
    hover = HoverTool(tooltips=[("ω (THz)", "$y")])
    p.add_tools(hover)
    
    ## Make the layout
    layout = column(eigen_chkbx_grp,p,row(slider,div))
    
    ## Make the webpage
    show(layout)