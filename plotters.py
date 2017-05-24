import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import analysis_tools as tools
import scipy.spatial as spatial
from scipy.stats import norm, rayleigh
from itertools import cycle
from scipy import stats
import numpy as np
from pylab import *
import matplotlib


def plot_bax_oligo(data, fit=True, verbose=False):
    color = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    colorcycle = cycle(color)
    for id, timeseries in data.items():
        acolor = colorcycle.next()
        x, y = tools.get_xy(timeseries)
        plt. plot(x, y, color = acolor, label= "mito" + ' ' + str(id), linewidth=2)
        if not fit:
            continue  
        if fit:
          xfit, yfit, _, _ = tools.bax_fit(timeseries, id)
          plt.plot(xfit, yfit, color = acolor, linewidth=2, label= "mito " + str(id) + " fit")
    
    plt.xlabel('Time (minutes)')
    plt.ylabel('Percent Baseline')
    plt.axhline(y = 0, linewidth=3, color='black')
    plt.axvline(x=0, linewidth=3, color="black")
    plt.xticks(size = 14)
    plt.yticks(size = 14)

    plt.legend(loc=2, shadow=False, prop={'size':8}, bbox_transform=plt.gcf().transFigure)
    #plt.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure, prop={'size':2})	    
    plt.title(title)    
    plt.show()

def plot_scatter(x, y):
    plt.scatter(x, y)
    plt.axhline(y = 0, linewidth=3, color='black')
    plt.axvline(x=0, linewidth=3, color="black")
    plt.ylim(0, ymax = 3)
    plt.xlim(0, xmax = 1000)
    plt.xlabel("Cell")
    plt.ylabel("Time (min)")
    plt.show()

def plot_oligo_double_axis(x, x2, y, y2):
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  ax1.plot(x, x2, 'o', color = 'black')
  ax1.set_ylabel('Rate of Oligomerization')
  ax2 = ax1.twinx()
  ax2.plot(y, y2, 'o',  color = 'r')
  ax2.set_ylabel('Time of initiation', color='r')
  for tl in ax2.get_yticklabels():
    tl.set_color('r')
  plt.xlim(0, max(x) +1)  
  plt.ylim(0, max(y2) +20)
  plt.show()


def plot_cytc(data, fit=True, verbose=False):
  color = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
  colorcycle = cycle(color)
  for id, timeseries in data.items():
    acolor = colorcycle.next()
    x, y = tools.get_xy(timeseries)
    #print x, y, id
    if not fit:
      plt.plot(x, y, color = acolor, label=str(id), linewidth=2)
      continue  
    if fit:
      xfit, yfit, popt = tools.cytc_fit(timeseries, id)
      if xfit != None:
        #plt.plot(x, y, xfit, yfit, color = acolor, linewidth=2, label="mito " + str(id) + " fit")
        i = tools.sigmoid_inflection(xfit, popt[0], popt[1], popt[2], popt[3], popt[4])
        slope_d1 = tools.sigmoid_d1(xfit[i], *popt)
        if slope_d1 == abs(0) or math.isnan(slope_d1) == True:
            print slope_d1, "slope_d1"
            max_line = sum(yfit[:20]) / len(yfit[:20])
            min_line = sum(yfit[480:]) / len(yfit[480:])
            y_inflection = ((max_line - min_line) / 2) + min_line
            fl = tools.crosses(yfit, y_inflection)
            x_inflection = xfit[fl]
            print max_line, min_line, x_inflection
            b = y_inflection - (x_inflection * 0)
            x_line = range(len(xfit))
            for point in range(len(x_line)):
                y_line = (0* point) + b   
            plt.plot(xfit, yfit, color= 'green')
            plt.axvline(x= x_inflection, color = 'black')
            plt.plot(x_inflection, y_inflection, 'o', color= 'g')
            plt.show()
        elif i != len(xfit) - 1 and i != 0:
            slope = tools.sigmoid_d1(xfit[i], *popt)
            yy = []
            yymin = []
            bb = yfit[i] - (xfit[i] * slope)
            xx = range(len(xfit))
            for vals in range(len(xx)):
              yys = (slope * vals) + bb
              yy.append(yys)
            xxmin = np.arange(-500, 500, 5)
            for vals in range(len(xxmin)):
              yymins = tools.sigmoid(vals, popt[0], popt[1], popt[2], popt[3], popt[4]) 
              yymin.append(yymins)
            cmin_line = (sum(yymin[180:]) / len(yymin[180:]))  
            cmax_line= cmin_line + (2 * (yfit[i] - cmin_line))   
            cytc_decay = tools.crosses(yy, cmax_line)
            plt.plot(xfit, yfit, color = 'blue', linewidth=2)
            plt.axvline(x= cytc_decay, linewidth=2, color= 'red')
            plt.plot(xfit[i], yfit[i], 'o', color='r')

            #cyt_decay = xfit[not_zero[0]]
            #cyt_sat = xfit[not_zero[-1]]
            print slope_d1, i, "=parangs"
            plt.xlabel('Time (minutes)')
            plt.ylabel('Percent Baseline')
            plt.axhline(y = 0, linewidth=3, color='black')
    #plt.axvline(x= cyt_decay, linewidth=3, color="black")
    # # plt.axvline(x= cytc_decay, linewidth=2, color= 'red')
    # plt.xticks(size = 14)
    # plt.yticks(size = 14)
    # plt.legend(loc=2, shadow=False, prop={'size':8}, bbox_transform=plt.gcf().transFigure)      #          continue              
            plt.show()

def plot_dist_rate(data, xyzs, center):
    ks = {}
    for id, timeseries in data.items():
        x, y = tools.get_xy(timeseries)
        try:
            _, _, k = tools.get_sigmoid_fit(x, y, verbose=True)
            ks[id] = k
        except RuntimeError:
            continue
    dists = tools.distance_from_center(xyzs, center)
    print dists
    x = []
    y = []
    for id, _ in ks.items():
        x.append(dists[id])
        y.append(ks[id])
    plt.scatter(x,y)
    plt.show()


def plot_compiled_dist_rate(xs, ys):
    slope, intercept, r_value, p_value, std_err = stats.linregress(xs, ys)
    line = [slope * x + intercept for x in xs]
    print "rval", r_value, "pval", p_value
    plt.xlabel('Distance from nucleus (um)')
    plt.ylabel('Rate of Oligomerization (RFU / min)')
    plt.plot(xs, line, 'r-', label = "r = 0.06", linewidth=2)
    plt.plot(xs, line, 'r-', label = "pval = 0.25", linewidth=2)
    plt.scatter(xs, ys, facecolor="black", linewidth=2)
    legend = plt.legend(loc=1)
    frame = legend.get_frame()
    frame.set_linewidth(2)
    plt.rcParams['axes.linewidth'] = 2
    plt.axhline(y = -0.5, linewidth=3, color='black')
    plt.axhline(y = 3, linewidth=3, color='black')
    plt.axvline(x=0, linewidth=3, color="black")
    plt.axvline(x=30, linewidth=3, color="black")
    plt.show()

def plot_agg_initiation(x):
    bins_bound = np.linspace(0, 1000, 100)
    #events, edges, patches = hist(slopes, bins)
    #plt.plot(bins)#num_bins = 100
    n, bins, patches = plt.hist(x, bins_bound, facecolor = 'white', linewidth = 0)
    norm_data = []
    # print bins, bins_bound
    for i in range(len(n)):
        norm = (n[i] / sum(n))
        print norm, "=norm"
        norm_data.append(norm)
    norm_data.append(0)
    print norm_data, sum(norm_data), len(norm_data)
    width = bins_bound[1]- bins_bound[0]
    print width
    plt.bar(bins_bound, norm_data, width, color= "gray")
    #plt.hist(x, bins=num_bins, color = 'black')
    plt.xlim(xmin = 0, xmax = 1000)
    plt.ylim(ymax = 0.5)
    plt.xlabel('Time (min)')
    plt.ylabel(' Percent of Frequency * 100  (of mitochondria)')
    plt.title('Time of Initiation After STS/Hdac3 Addition')
    plt.show()
    print len(x)

def plot_oligo_rate_histo(norm_data):   
    # norm_data = []
    bins_bound = np.linspace(0, 5, 36)
    # n, bins, patches = plt.hist(slopes, bins_bound, normed=True, stacked = True)
    # total = sum(n)
    # for i in range(len(n)):
    #     norm = (n[i] / total)
    #     print norm, "=norm"
    #     norm_data.append(norm)
    # print norm_data, sum(norm_data)
    plt.hist(norm_data, bins=bins_bound, color = 'black')
    plt.xlabel('Rate of BAX Oligomerization\n(RFU min^-1)' , size = 16)
    plt.ylabel('Frequency\n(number of mitochondria)', size = 16)
    plt.legend()
    #plt.title('Rate of Aggregation per Mitochondria of One Cell ')
    #plt.ylim(ymax = 2)
    plt.xlim(xmax = 1)
    #plt.axes(frameon = False)
    #ax1.get_xaxis().tick_bottom()
    #ax1.get_yaxis().tick_left()
    #plt.rcParams['axes.linewidth'] = 2
    plt.axhline(y = 0, linewidth=3, color='black')
    plt.axvline(x=0, linewidth=3, color="black")
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.minorticks_off()
    plt.tick_params('both', length=8, width=2, which='major')
    plt.show() 


def plot_cell_start(xs, ys, xfit, yfit, k, xinfl, yinfl, y0):
    plt.plot(xs, ys, linewidth=2)
    plt.plot(xfit, yfit, linewidth=2)
    line = [k*x + yinfl - k*xinfl for x in xfit]
    # plt.plot(xfit, line, label="Maximum Rate", linewidth=2)
    y0s = [y0 for x in xfit]
    # plt.plot(xfit, y0s, label="Minimum Boundary", linewidth=2, color= "black")
    plt.ylim([0, 5])
    plt.xlim([0, 100])
    plt.xlabel("Relative Time (min)")
    plt.ylabel("Relative Fluorescence")
    plt.axhline(y=0, linewidth = 2.5, color= 'black')
    plt.axvline(x=0, linewidth = 2.5, color= 'black')

    # legend = plt.legend(loc=2)
    # frame = legend.get_frame()
    # frame.set_linewidth(2)
    # plt.rcParams['axes.linewidth'] = 2
    # plt.show()
