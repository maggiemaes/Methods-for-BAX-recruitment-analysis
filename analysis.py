from plotters import plot_bax_oligo, plot_dist_rate, plot_compiled_dist_rate, plot_agg_initiation, plot_oligo_rate_histo, plot_cytc, plot_scatter, plot_oligo_double_axis
import analysis_tools as tools
from scipy.spatial import distance
import matplotlib.pyplot as plt
import numpy as np
import data
import math
from itertools import cycle
from scipy import stats
from scipy.stats import norm
from collections import defaultdict
import os
import pandas as pd



def bax_oligo_fit(): 
    prefix="./data/"
    suffix=".csv"
      
    data_files = data.bax_oligo_data #Import data files from data.py
    bax_files = data_files.keys()
    ids = data_files.values()
    all_data = []
    bax_data = {}
    for i, f in enumerate(bax_files):
      fname = prefix + f + suffix 
      print f
      rows = tools.get_rows(fname, 4)
      cell_data = tools.get_val_time(rows, ids[i])
      tools.filter_1(cell_data)
      tools.normalize_val_time(cell_data)
      bax_data[f] = cell_data
      all_data.append(cell_data)
      print cell_data

    nfiles = len(bax_files)
    color = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    colorcycle = cycle(color)

    for i in range(nfiles):
      bax_file = bax_files[i]
      print bax_files[i]
      ids_list = ids[i]
      bdata = bax_data[bax_file]
      print bdata
      plot_bax_oligo(bdata, fit=True, verbose=True)
 

def cytc_release_fit():
    data_ch2 = data.bax_cytc_datach2
    cytc_files = data_ch2.keys()
    ids = data_ch2.values()
    
    prefix = "../cyt_c/channel2/"
    suffix = ".csv"

    cytc_data = {}
    for i, f in enumerate(cytc_files):
      fname = prefix + f + suffix
      rows = tools.get_rows(fname, 4)
      cell_data = tools.get_val_time(rows, ids[i])
      tools.filter_1(cell_data)
      tools.normalize_val_time(cell_data)
      cytc_data[f] = cell_data

    nfiles = len(cytc_files)
    color = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    colorcycle = cycle(color)

    for i in range(nfiles):
      cytc_file = cytc_files[i]
      ids_list = ids[i]
      cdata = cytc_data[cytc_file]
      print cytc_file
      plot_cytc(cdata, fit=True, verbose=True)

def retro_rates():
    prefix="../bax_agg/"
    suffix=".csv"
      
    data_files = data.bax_oligo_data #Import data files from data.py
    
    names = []
    ids_list = []
    on_slopes = []
    off_slopes = []
    param_list =[]
    for name, ids in data_files.items(): #Create lists of file names and corresponding ids
      names.append(name)
      ids_list.append(ids)

    for i in range(len(names)):
      fname = prefix+names[i]+suffix 
      print names[i]
      rows = tools.get_rows(fname, 4)
      val_time = tools.get_val_time(rows, ids_list[i])
      tools.filter_1(val_time)
      tools.normalize_val_time(val_time)
      scatterxval = i + 1
      cell_rate = [] 
      for id, timeseries in val_time.items():
          x, y = tools.get_xy(timeseries)
          try:
            xfit, yfit, popt, pcov = tools.retro_fit(timeseries, id)
            if popt is not None and not isinstance(pcov,float):
              split = (int(popt[1]) - x[0]) / 5
              print split, int(popt[1]), x[0], '=x0'
              front_timeseries = timeseries[:split]
              back_timeseries = timeseries[split:]
              # print front_timeseries, back_timeseries
              if len(front_timeseries) > 5:
                try:
                  x1fit, y1fit, popt1, pcov1 = tools.bax_fit(front_timeseries, id)
                  if popt1 != None: 
                    i = tools.sigmoid_inflection(x1fit, popt1[0], popt1[1], popt1[2], popt1[3], popt1[4])
                    if popt[1] > 0 and i != len(xfit) - 1:
                      slope = tools.sigmoid_d1(x1fit[i], *popt1)
                      print slope, "= slope", id
                      log = math.log(slope, 2)
                      on_slopes.append(log)
                      x0, k, y0, a, b = popt1[0], popt1[1], popt1[2], popt1[3], popt1[4]
                      param_data = {'id': id, 'a': a, 'b': b, 'k': k, 'x0': x0, 'y0': y0}
                      param_list.append(param_data)
                      plt.plot(x1fit, y1fit)   
                except RuntimeError:
                  continue                  
              if len(back_timeseries) > 5:
                try:
                  x2fit, y2fit, popt2 = tools.cytc_fit(back_timeseries, id)
                  print popt2, 'popt2', fname, id
                  if popt2 is not None and x2fit != None:
                    i2 = tools.sigmoid_inflection(x2fit, popt2[0],popt2[1], popt2[2],popt2[3], popt2[4])
                    if popt2[1] > 0 and i2 != len(x2fit) - 1:
                      slope2 = tools.sigmoid_d1(x2fit[i], *popt2)
                      if slope2 != 0: 
                        print slope2, "= slope2", id
                        print abs(slope2)
                        log2 = math.log((abs(slope2)), 2)
                        off_slopes.append(log2)
                        # x0, k, y0, a, b = popt2[0], popt2[1], popt2[2], popt2[3], popt2[4]
                        # param_data = {'id': id, 'a': a, 'b': b, 'k': k, 'x0': x0, 'y0': y0}
                        # param_list.append(param_data)
                        plt.plot(x2fit, y2fit)      
                  else: 
                    continue
                except RuntimeError:
                  continue    
              plt.axvline(x= popt[1])
              plt.plot(x, y, xfit, yfit)
              # plt.show()   
            else: 
              continue
          except RuntimeError: 
            continue    
    print on_slopes, off_slopes
    df = pd.DataFrame(param_list)
    # df.to_csv('Retro_ON.csv')

def oligo_rates():
    prefix="../bax_agg/"
    suffix=".csv"
      
    data_files = data.bax_oligo_data #Import data files from data.py
    
    names = []
    ids_list = []
    all_slopes = []
    log_slopes = []
    indiv_rates = {}
    cell_stds = []
    scatterx = []
    scattery = []
    popt_k = []
    popt_a =[]
    popt_b = []
    param_list =[]
    for name, ids in data_files.items(): #Create lists of file names and corresponding ids
      names.append(name)
      ids_list.append(ids)

    for i in range(len(names)):
      fname = prefix+names[i]+suffix 
      print names[i]
      rows = tools.get_rows(fname, 4)
      val_time = tools.get_val_time(rows, ids_list[i])
      tools.filter_1(val_time)
      tools.normalize_val_time(val_time)
      scatterxval = i + 1
      cell_rate = [] 
      for id, timeseries in val_time.items():
          x, y = tools.get_xy(timeseries)
          try:
            xfit, yfit, popt, pcov = tools.bax_fit(timeseries, id)
            if popt != None and not isinstance(pcov,float):
              print yfit[-10], yfit[-1], popt
              i = tools.sigmoid_inflection(xfit, popt[0], popt[1], popt[2], popt[3], popt[4])
              # print i, popt[1], len(xfit)
              if popt[1] > 0 and i != len(xfit) - 1:
                slope = tools.sigmoid_d1(xfit[i], popt[0], popt[1], popt[2], popt[3], popt[4])
                print slope, "= slope", id
                if slope > 0 and tools.filter_2(xfit, yfit) == None:
                  print fname, id
                  # plt.ylim(0,ymax=5)
                  # plt.plot(x, y, xfit, yfit)
                  # plt.show()
                  cell_rate.append(slope)
                  all_slopes.append(slope)
                  scatterx.append(scatterxval)
                  scattery.append(slope)
                  popt_k.append(popt[1])
                  popt_a.append(popt[3])
                  popt_b.append(popt[4])
                  x0, k, y0, a, b = popt[0], popt[1], popt[2], popt[3], popt[4]
                  x0_cov, k_cov, y0_cov, a_cov, b_cov = pcov[0,0]**0.5, pcov[1,1]**0.5, pcov[2,2]**0.5, pcov[3,3]**0.5, pcov[4,4]**0.5
                  param_data = {'id': id, 
                    'a': a, 'b': b, 'k': k, 'x0': x0, 'y0': y0,
                    'x0_cov': x0_cov, 'k_cov': k_cov, 'y0_cov': y0_cov, 'a_cov': a_cov, 'b_cov': b_cov}
                  param_list.append(param_data)
          except RuntimeError:
            continue  
      if sum(cell_rate) != 0: 
        cell_avg = sum(cell_rate)/len(cell_rate)
        cell_std = np.std(cell_rate)
          # print len(cell_rate), cell_avg, "=cellavg"
        cell_stds.append(cell_std)  
        indiv_rates[fname] = [cell_std]          
    for i in range(len(all_slopes)):  # log transformation for comparison between cell types
      log = math.log(all_slopes[i], 2)
      print log, all_slopes[i]
      if log < (-5.2 + (2*29)) and log > (-5.2 - (2*29)):
        log_slopes.append(log) 
      else: 
       continue  

    avg_rate = sum(log_slopes)/len(log_slopes)
    std = np.std(log_slopes)
    print log_slopes
    print "Avg=", avg_rate, "Std=", std, "n=", len(log_slopes), len(cell_rate), len(all_slopes)
    norm_data = []
    bins_bound = np.linspace(-10, 10, 51)
    n, bins, patches = plt.hist(log_slopes, bins_bound, color= 'white', linewidth = 0)
    for i in range(len(n)):
        norm_log = (n[i] / sum(n))
        norm_data.append(norm_log)
    norm_data.append(0)
    print norm_data, sum(norm_data), len(norm_data)
    width = bins_bound[1]- bins_bound[0]
    plt.bar(bins_bound, norm_data, width, color= "gray")
    plt.ylim(ymax=0.4)
    plt.xlabel("Log(Rate of Oligomerization) ")
    plt.ylabel("Percent Frequency within bin *100")
    # plt.scatter(scatterx, scattery)
    plt.show()

    df = pd.DataFrame(param_list)
    # df.to_csv('RGC_ONC.csv')
    avg_k = sum(popt_k) / len(popt_k)
    avg_a = sum(popt_a) / len(popt_a)
    avg_b = sum(popt_b) / len(popt_b)
    print 'k=', avg_k, 'a=', avg_a, 'b=', avg_b
   #For yMax Values of sigmoid curves. 




    #plot_oligo_rate_histo(norm_data)
    # cell1 = []
    # cell2 = []
    # cell3 = []
    # for i in range(len(scatterx)):
    #   if scatterx[i] == 1: 
    #     cell1.append(scattery[i])
    #   if scatterx[i] == 2: 
    #     cell2.append(scattery[i])
    #   if scatterx[i] == 3: 
    #     cell3.append(scattery[i])
    # mean1 = np.mean(cell1)
    # mean2 = np.mean(cell2)
    # mean3 = np.mean(cell3)      
    # std1 = np.std(cell1)
    # std2 = np.std(cell2)
    # std3 = np.std(cell3)
    # print "cell1=", len(cell1), mean1, std1, "cell2=", len(cell2), mean2, std2, "cell3=", len(cell3), mean3, std3
    
     ##For Manhattan plots
    # HCTDKO = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0086206896551724137, 0.0, 0.025862068965517241, 0.0, 0.034482758620689655, 0.077586206896551727, 0.077586206896551727, 0.13793103448275862, 0.20689655172413793, 0.17241379310344829, 0.094827586206896547, 0.051724137931034482, 0.043103448275862072, 0.034482758620689655, 0.017241379310344827, 0.0086206896551724137, 0.0086206896551724137, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0]    
    # HCTBAXKO = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.012195121951219513, 0.0, 0.024390243902439025, 0.036585365853658534, 0.024390243902439025, 0.14634146341463414, 0.073170731707317069, 0.18292682926829268, 0.12195121951219512, 0.12195121951219512, 0.06097560975609756, 0.04878048780487805, 0.073170731707317069, 0.012195121951219513, 0.036585365853658534, 0.012195121951219513, 0.012195121951219513, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0]
    # rgc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01048951048951049, 0.034965034965034968, 0.048951048951048952, 0.080419580419580416, 0.11888111888111888, 0.12237762237762238, 0.15034965034965034, 0.13636363636363635, 0.097902097902097904, 0.062937062937062943, 0.1048951048951049, 0.031468531468531472, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0]
    # hct = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0029498525073746312, 0.017699115044247787, 0.041297935103244837, 0.041297935103244837, 0.067846607669616518, 0.091445427728613568, 0.15929203539823009, 0.16224188790560473, 0.15339233038348082, 0.13569321533923304, 0.076696165191740412, 0.038348082595870206, 0.011799410029498525, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0]
    # d407 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0012804097311139564, 0.0, 0.0012804097311139564, 0.0, 0.0025608194622279128, 0.0012804097311139564, 0.0025608194622279128, 0.0064020486555697821, 0.017925736235595392, 0.014084507042253521, 0.040973111395646605, 0.060179257362355951, 0.10755441741357234, 0.17029449423815621, 0.17285531370038412, 0.16133162612035851, 0.14084507042253522, 0.04353393085787452, 0.024327784891165175, 0.012804097311139564, 0.0051216389244558257, 0.0038412291933418692, 0.0038412291933418692, 0.0025608194622279128, 0.0025608194622279128, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0]
    # mCherry = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0027548209366391185, 0.0013774104683195593, 0.011019283746556474, 0.0082644628099173556, 0.022038567493112948, 0.033057851239669422, 0.078512396694214878, 0.14325068870523416, 0.14462809917355371, 0.16115702479338842, 0.13085399449035812, 0.0743801652892562, 0.071625344352617082, 0.048209366391184574, 0.033057851239669422, 0.01790633608815427, 0.0055096418732782371, 0.0041322314049586778, 0.0041322314049586778, 0.0027548209366391185, 0.0, 0.0013774104683195593, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0]
    # # samp = norm.rvs(loc=0,scale=1,size=len(HCTDKO)) 
    # # print log_slopes, norm_data
    # mean, std_n = norm.fit(log_slopes)
    # print mean, std_n
    # plt.plot(bins_bound, norm.pdf(bins_bound, loc=mean, scale=std_n))
    # # bins_bound = np.linspace(-10, 10, 51)
    # # pdf_fitted = norm.pdf(bins_bound,loc=avg_rate,scale=std)
    # # plt.plot(bins_bound, rgc, ls = 'steps', label = 'RGC', linewidth = 2, color = 'red')
    # plt.plot(bins_bound, HCTDKO, ls= 'steps', label = "HCT116 DKO", linewidth = 2, color = 'red')
    # # plt.plot(bins_bound, HCTBAXKO, ls = 'steps', label = 'HCT116 BAX KO', linewidth = 2, color = 'blue')
    # plt.plot(bins_bound, hct, ls = 'steps', label = 'HCT116 WT', linewidth = 2, color = 'green')
    # plt.axhline(y=0, linewidth=2, color= 'black')
    # plt.axvline(x=-10, linewidth=2, color= 'black')
    # plt.legend(loc=2, shadow=False, prop={'size':12}, bbox_transform=plt.gcf().transFigure)
    # plt.xlabel('Rate of oligomerization \n log(2)RFU/min')
    # plt.ylabel('Frequency \n Percent of Total')
    # plt.ylim(ymax= 0.5)
    # plt.title('Rate of BAX Oligomerization')
    # plt.show()     


def initiation_time():  # Per Mitochondria And Per Cell
    prefix="../bax_agg/"
    suffix=".csv"
    data_files = data.bax_oligo_data #Import data files from data.py    
    names = []
    ids_list = []
    time_oligo = [] #initiation time. 
    time_complete = [] #completion time.
    scattery = []
    log_y = []
    avg_ymax =[]
    cell_rates = defaultdict(list)   ## key: [vals] = fname: [rates of ea. mito]
    cell_inits = defaultdict(list)
    cell_compls = defaultdict(list)  ###fname: [total completion time of mito]
    log_inits = defaultdict(list)
    cell_pcntrates = defaultdict(list) ## fname: [%rate or total cellrates of ea. mito]
    cell_pcntinits = defaultdict(list)
    avg_cell_rates = []     #list of avg rates(mito)/cell for all cells
    avg_cell_inits = []
    avg_cell_completions = []
    std_pcntof_avgrates = []         #list of std as percent over average per cell
    std_pcntof_avginits = []         #list of std as percent over average per cell
    std_cell_rates = []     #list of stdevs(mito)/cell for all cells
    std_cell_inits = []
    cell_number = []
    cell_number2 = []
    duration = []
    plos_list = []
    for name, ids in data_files.items(): #Create lists of file names and corresponding ids
      names.append(name)
      ids_list.append(ids)
    for i in range(len(names)): 
      print (i, names[i])
      label = (names[i])
      rows1 = tools.get_rows("../bax_agg/Image Acquisition Times.csv", 4)
      lag_name, lag_time = tools.get_lagtime(rows1, names[i])      
      #print lag_name, lag_time, "da;lkfsdjf"
      fname = prefix+names[i]+suffix 
      rows = tools.get_rows(fname, 4)
      val_time = tools.get_val_time(rows, ids_list[i])
      tools.filter_1(val_time)
      tools.normalize_val_time(val_time)
      for lagger in range(len(lag_name)):
      #print lag_name[lagger], names[i]
        if names[i] == lag_name[lagger]:
          cell_lag = lag_time[lagger]
          print lag_name[lagger], cell_lag
        else:
          continue
      # cell_lag = 50
      for id, timeseries in val_time.items():
        x, y = tools.get_xy(timeseries)
        # plt.plot(x, y, color = 'blue')
        xfit, yfit, popt, pcov = tools.bax_fit(val_time[id], id)
        if xfit == None:
          print "No fit"  
        else:
          i = tools.sigmoid_inflection(xfit, *popt)     #i = xfit where crosses zero. xfit[i] = inflection point.
          # print "SFIts", yfit[0], yfit[-1]
          # plt.plot(xfit, yfit, x, y)
          # plt.show()
          if popt[1] > 0 and i != len(xfit) - 1 and  tools.close(yfit[0], yfit[-1]) == False:
            # plt.plot(xfit, yfit, x, y)
            # plt.show()
            slope = tools.sigmoid_d1(xfit[i], *popt)
            y = []
            ymin = []  
            b = yfit[i] - (xfit[i]*slope)  #creating separate curve to determine max/min lines.
            xx=range(len(xfit))
            for value in range(len(xx)): 
              ys = (slope*value) + b
              y.append(ys)
            xmin = np.arange(-500, 500, 5)
            for value in range(len(xmin)):  
              ymins = tools.sigmoid(value, popt[0], popt[1], popt[2], popt[3], popt[4])
              ymin.append(ymins)


            max_line = (sum(ymin[180:]) / len(ymin[180:])) 
            minimumval = sum(ymin[:20]) / len(ymin[:20]) #to calc max Y value
            delta_y = (max_line - minimumval) + 1  #change in y plus norm value 1
            min_line= max_line - (2 * (max_line - yfit[i]))
            init_ymxb = tools.crosses(y, min_line)
            sat = tools.crosses(y, max_line)
            print max_line, '=maxline', min_line, minimumval, delta_y, '=deltay'
            # changey = min_line + ((max_line- min_line) * 0.05)
            # init_chan = tools.crosses(yfit, changey)
            # init_chanfix = (init_chan / 500) * len(x) + (x[0] - 1)
            # print init_chan, '=initchangey'
            # max_thry = (max_line - 0.05 * max_line)   # max threshold of extended curve. y value. 
            total_time = init_ymxb + cell_lag
            oligo_dur = int(sat-init_ymxb)
            if slope > 0 and math.isnan(slope) == False and init_ymxb != 499 and tools.close(slope, 0) == False and tools.close(oligo_dur, 0) == False:
              # print oligo_dur, sat, init_ymxb
              # log_dur = math.log(oligo_dur, 2)
              time_oligo.append(total_time)  # for population data
              cell_inits[label].append(total_time)  #for per cell data
              scattery.append(slope)        # for population data
              log_num = math.log(slope, 2)
              #log_init = math.log(total_time, 2)
              # if log_num < (-4 + (2*15)) and log_num > (-4 - (2*15)):
              cell_rates[label].append(log_num)       # for per cell data
              plos_data1 = {'Cell': label, 'id': id, "Mito_Rate": log_num, 'Mito_Initiation': total_time}
              plos_list.append(plos_data1) 
              # else: 
              #   continue
              time_complete.append(sat + cell_lag)  #population data
              #log_inits[label].append(log_init)
              cell_compls[label].append(oligo_dur)
              # duration.append(log_dur)   #log2 of duration of oligomerization creates normal distribution 
              avg_ymax.append(delta_y)   # for ymax values avg next to plot compiled curves
              # print init_ymxb, oligo_dur, avg_ymax, 'avg ymax'
            else: 
              continue
            print slope, '=slope'
    #        # print 'id=', id,  'init_ymxb', init_ymxb, 'total_time=', total_time, cell_lag, 'completion==', sat 
            # plt.plot(xfit, yfit, linewidth=2, color='blue') #x, y gives slope line
            # plt.plot(xmin, ymin, linewidth=2, color= 'green')
            # plt.scatter(xfit[i], yfit[i], color='red')
            # plt.axhline(y = max_line, color='black', linewidth=1)
            # plt.axhline(y = minimumval, color = 'black', linewidth=1)
            # # plt.axvline(x= init_chanfix, color= 'red', linewidth=1)
            # plt.axvline(x = init_ymxb, color  = 'black', linewidth=1)
            # plt.show() 
    avg_oligo = sum(time_oligo) / len(time_oligo)
    std_oligo = np.std(time_oligo)
    # print 'intiation time=', time_oligo, avg_oligo, std_oligo, len(time_oligo), cell_lag   
        # cell_rate_std = np.std(cell_rates)
        # cell_init_std = np.std(cell_inits)    
      # print names[i], "rate std", cell_rate_std, "initiation std", cell_init_std
    for item in range(len(names)):    # gets values for individual cells
      # print item, "fname, len of mitos, rates", names[item], len(cell_rates[names[item]])
      if len(cell_rates[names[item]]) >1 and len(cell_inits[names[item]]) >1 and len(cell_compls[names[item]]) >1:
        avg_ratepc = sum(cell_rates[names[item]]) / len(cell_rates[names[item]])  # avg of mitochondria in one cell
        std_ratepc = np.std(cell_rates[names[item]])                              #std of mito in one cell
        # print avg_ratepc, '=avg_key', names[item], std_ratepc
        avg_initpc = sum(cell_inits[names[item]]) / len(cell_inits[names[item]])
        std_initpc = np.std(cell_inits[names[item]])
        # print avg_initpc, std_initpc, names[item], #cell_rates[names[item]]
        avg_complpc = sum(cell_compls[names[item]]) / len(cell_compls[names[item]])
        std_complpc = np.std(cell_compls[names[item]])
        # print avg_complpc, std_complpc, "compeltions per cell"
        # plos_data2 = {'Cell': names[item], 'Cell_rate': avg_ratepc, 'Std_rate': std_ratepc, 'Cell_initiation': avg_initpc, 'Std_initiation': std_initpc}
        # plos_list.append(plos_data2) 
        avg_cell_rates.append(avg_ratepc) # list of avg rate within a single cell
        avg_cell_inits.append(avg_initpc)       
        avg_cell_completions.append(avg_complpc)
        std_cell_rates.append(std_ratepc)  #list of std of rates within a single cell 
        std_cell_inits.append(std_initpc)
        qs3, qs4, qs5, qs6 = avg_ratepc, std_ratepc, avg_initpc, std_initpc
        std_pcntof_avgrates.append(abs((std_ratepc) / (avg_ratepc)))  # COEFF OF VARIANCE =stdevs as a percentage of the cell average rate.
        std_pcntof_avginits.append(abs((std_initpc) / (avg_initpc)))
        # print std_pcntof_avgrates, fname, avg_cell_rates
    avg_pcrates = sum(avg_cell_rates)/ len(avg_cell_rates)  #total avg of rates/cells for avg line in graph
    avg_pcinits = sum(avg_cell_inits)/ len(avg_cell_inits)
    std_pcrates = np.std(avg_cell_rates)                    #std of total avgs. dont really need this
    std_pcinits = np.std(avg_cell_inits)
    avg_std_pcrates = sum(std_cell_rates) / len(std_cell_rates)  #avg of stdevs 
    avg_std_pcinits = sum(std_cell_inits) / len(std_cell_inits)
    ss_rates = np.std(std_cell_rates)                            #std of stdevs.
    ss_inits = np.std(std_cell_inits)
    avg_pcnt_pcrates = sum(std_pcntof_avgrates) / len(std_pcntof_avgrates)       #total avg of Variance(%cell rates/cell)
    avg_pcnt_pcinits = sum(std_pcntof_avginits) / len(std_pcntof_avginits)
    std_pcnt_pcrates = np.std(std_pcntof_avgrates)                            #std of stdevs.
    std_pcnt_pcinits = np.std(std_pcntof_avginits)
    cell_number.append(range(len(std_pcntof_avgrates)))
    # print "list of completions", avg_cell_completions, (sum(avg_cell_completions)/len(avg_cell_completions)), np.std(avg_cell_completions)
    print len(cell_number), len(std_pcntof_avgrates), std_pcntof_avgrates, avg_pcnt_pcrates
    plt.scatter(cell_number, std_pcntof_avgrates, label= "Rate Std/Avg (as percent)", color = 'blue', marker = 'o')
    plt.scatter(cell_number, std_pcntof_avginits, label= "Initiation Std/Avg (as percent)", color = 'red', marker = 'o')
    plt.axhline(y=avg_pcnt_pcrates, label="Avg cell rate", color = 'blue')
    plt.axhline(y=avg_pcnt_pcinits, label= "avg cell Initiation", color = 'red')
    plt.title("Rate/Init as Percent std/avg per cell")
    plt.legend(loc=1, shadow=False, prop={'size':8}, bbox_transform=plt.gcf().transFigure)
    plt.xlim(-5, 75)
    plt.ylim(-1, 2)
    plt.show()
    avg_ymaxs = sum(avg_ymax)/len(avg_ymax)
    std_ymax = np.std(avg_ymax)
    print "values needed", avg_pcrates, std_pcrates, avg_pcinits, std_pcinits, "avg,std of varianece",avg_pcnt_pcrates, std_pcnt_pcrates, avg_pcnt_pcinits, std_pcnt_pcinits 
    # print "per cell info", avg_cell_rates, avg_cell_inits, std_cell_inits, avg_ymaxs, std_ymax
    for it in range(len(scattery)):  # log transformation for comparison between cell types 
        log = math.log(scattery[it], 2)
        # if log < (-5 + 2*29) and log > (-5 - 2*29):
        log_y.append(log)
        # else:
        #   time_oligo.remove(time_oligo[it])
    plt.plot(time_oligo, log_y, 'o', color= 'blue')
    # plt.axhline(y = 0, linewidth=3, color='black')
    # plt.axvline(x=0, linewidth=3, color="black")
    plt.ylim(-10, ymax = 5)
    plt.xlim(0, xmax = 1000)
    plt.xlabel("Initiation Time (min)")
    plt.ylabel("Rate of Oligomerization (RFU/min)")
    plt.show()
    plot_agg_initiation(time_oligo)
    # outliers_rem = []          #creating list of completions without the outlier values
    # for t in range(len(avg_cell_completions)):
    #   print avg_cell_completions[t]
    #   if 3.6< avg_cell_completions[t] < 17.4:
    #     outliers_rem.append(avg_cell_completions[t])
    #     print 'add'
    #   else: 
    #     continue
    # new_avg = (sum(outliers_rem))/(len(outliers_rem))
    # new_std = np.std(outliers_rem)
    # print new_avg, new_std, len(outliers_rem), outliers_rem
    # dur_avg = (sum(duration))/(len(duration))
    # dur_std = np.std(duration)
    # print 'a', dur_avg, dur_std, len(duration)
    # plt.hist(duration)
    # plt.show()

    df = pd.DataFrame(plos_list)
    print plos_list
    df.to_csv('HCTKO_all.csv')

import plotters
def cell_time_plot():
  ids_list = [
    [7, 87, 119, 160, 225, 228, 235, 262, 284, 289, 301, 304],
    ]
  prefix="../bax_agg/"
  suffix=".csv"
  files = ["022014 D407 Bax cell1_Detailed",]

  for i in range(len(files)):
    print (i, files[i])
    fname = prefix+files[i]+suffix 
    rows = tools.get_rows(fname, 4)
    val_time = tools.get_val_time(rows, ids_list[i])
    tools.normalize_val_time(val_time)
    for id, timeseries in val_time.items():
        x, y = tools.get_xy(timeseries)
        try:
            xfit, yfit, popt, _ = tools.get_sigmoid_fit(x, y)
            i = tools.sigmoid_inflection(xfit, *popt)
            if popt[1] > 0 and i != len(xfit) - 1:
                slope = tools.sigmoid_d1(xfit[i], *popt)
                plotters.plot_cell_start(x, y, xfit, yfit, slope, xfit[i], yfit[i], popt[2])
        except RuntimeError:
          continue    

def init_diff(bdata, cdata):
  diff = {}
  ids = bdata.keys()
  bvals = bdata.values()
  cvals = cdata.values()
  for i in ids: 
    print(i)
    bx, by, bpopt, bpcov = tools.bax_fit(bdata[i], i)
    cx, cy, cpopt = tools.cytc_fit(cdata[i], i)
    _, sat = tools.start_time_fivepercent(cx, cpopt[3], cpopt[0], cpopt[1], cpopt[2])  
    thr, _ = tools.start_time_fivepercent(bx, bpopt[3], bpopt[0], bpopt[1], bpopt[2])  
    diff[id] = sat - thr   
  return diff

def bax_cytc_orig(): 
    data_ch1 = data.bax_cytc_datach1 #Import data files from data.py
    data_ch2 = data.bax_cytc_datach2

    bax_files = data_ch1.keys()
    cytc_files = data_ch2.keys()
    ids = data_ch1.values()
    ids2 = data_ch2.values()

    ch1_prefix="../cyt_c/channel1/"
    ch2_prefix="../cyt_c/channel2/"
    suffix=".csv"
    bax_data = {}
    cytc_data = {}
    filt_ids = {}  #ids filtered out for >5 data points.
    time_re_bax_initiation = []
    release_rate = []
    bax_rate  = []
    inf_slopes = []
    slope_d1s = []
    for i, f in enumerate(bax_files): 
      fname = ch1_prefix + f + suffix
      rows = tools.get_rows(fname, 4)
      cell_data = tools.get_val_time(rows, ids[i])
      tools.normalize_val_time(cell_data)
      tools.filter_1(cell_data)
      filt_ids[f] = cell_data.keys()
      bax_data[f] = cell_data     
    for i, f in enumerate(cytc_files):
      fname = ch2_prefix + f + suffix
      rows = tools.get_rows(fname, 4)
      cell_data = tools.get_val_time(rows, ids2[i])
      tools.normalize_val_time(cell_data)
      tools.filter_1(cell_data)
      for id, timeseries in cell_data.items():
        cell_data[id] = [[entry[0] / timeseries[-1][0], entry[1]] for entry in timeseries]
      cytc_data[f] = cell_data 
    
    #print "aa", filt_ids, "old", ids      
    for i in range(len(bax_files)):
      bax_file = bax_files[i]
      cytc_file = bax_file.replace("ch1", "ch2")
      print "file=", bax_file, cytc_file
      ids_list = filt_ids[bax_file]
      bdata = bax_data[bax_file]
      cdata = cytc_data[cytc_file]
      bax_initiate = float('NaN') 
      #diff = init_diff(bdata, cdata)
      #print(diff)
      for id in ids_list: 
        bx, by = tools.get_xy(bdata[id])
        cx, cy = tools.get_xy(cdata[id])
        bgood = True
        cgood= True 
        if bgood:
          bxfit, byfit, popt, pcov = tools.bax_fit(bdata[id], id)
          if bxfit == None:
            bgood = False
            print "bgood = False"  
          else:
            i = tools.sigmoid_inflection(bxfit, *popt)  #inflection point of fit_curve
            print "len_fit", len(bxfit), "len x/y", len(bx), len(by), "fit_inf=", i
            bslope_d1 = tools.sigmoid_d1(bxfit[i], *popt)
            print bslope_d1, "bslope_d1"
            if abs(bslope_d1) == 0 or math.isnan(bslope_d1) == True or tools.close(bslope_d1, 0):
              max_line = sum(byfit[:20]) / len(byfit[:20])
              min_line = sum(byfit[480:]) / len(byfit[480:])
              y_inflection = ((max_line - min_line) / 2) + min_line
              fl = tools.crosses(byfit, y_inflection)    #point where fit curve crosses infleciton point. 
              x_inflection = bxfit[fl]                  #defines x value of inflection point 
              print max_line, min_line, x_inflection 
              bax_initiation = x_inflection
              fast_rate = ((max_line - min_line) / (1))
              bax_rate.append(fast_rate)   
            else:
              y = []
              ymin = []  
              b = byfit[i] - (bxfit[i] * bslope_d1)  #creating separate curve to determine max/min lines.
              #print xinf, yinf, bxfit[i], byfit[i], 'here'
              x=range(len(bxfit))
              for value in range(len(x)): 
                ys = (bslope_d1*value) + b
                y.append(ys)
              xmin = np.arange(-500, 500, 5)
              for value in range(len(xmin)):  
                ymins = tools.sigmoid(value, popt[0], popt[1], popt[2], popt[3], popt[4])
                ymin.append(ymins)
              max_line2 = (sum(ymin[180:]) / len(ymin[180:]))  
              min_line2= max_line2 - (2 * (max_line2 - byfit[i]))
              print min_line2, max_line2, "max2"
              bax_initiation = tools.crosses(y, min_line2)
              bax_rate.append(bslope_d1)
          print "bax initiation", bax_initiation, bax_rate
              # print "fit_inf=", i, bxfit[i], popt[1], yinf, "bslope=", bslope, thr, "baxinitiate", bax_initiate, init_ymxb, "=init_ymxb"
              # plt.plot(bxfit, byfit, linewidth = 2)
              # #plt.scatter(xinf, yinf, color = 'red')
              # plt.axvline(x = init_ymxb, color = 'black')
              # plt.axvline(x = bax_initiate, color = 'green')
              # plt.axhline(y= max_line, color = 'black')
              # plt.axhline(y = min_line, color = 'black')
              # #plt.plot(x, y, color= 'red')
              # plt.show()

          if cgood:
            cxfit, cyfit, cpopt = tools.cytc_fit(cdata[id], id)
            if cxfit == None: 
              cgood = False
              print "No Cytc Fit"
              bax_rate.remove(bax_rate[-1])
            else: 
              cinf = tools.sigmoid_inflection(cxfit, *cpopt)    #determine inflection index
              # if cinf != len(cxfit) - 1 and cinf != 0:
              slope_d1 = tools.sigmoid_d1(cxfit[cinf], *cpopt)  #first derivative of inflection index gives slope of line at inflection point
              max_line = sum(cyfit[:20]) / len(cyfit[:20])
              min_line = sum(cyfit[480:]) / len(cyfit[480:])
              y_inflection = ((max_line - min_line) / 2) + min_line
              fl = tools.crosses(cyfit, y_inflection)    #point where fit curve crosses infleciton point. 
              x_inflection = cxfit[fl]                  #defines x value of inflection point 
              slope_decay = (min_line - max_line) / (1) #1 is max x distance for dataset
              print slope_decay, "=decay", slope_d1, "=d1"
              if slope_d1 > 0 and slope_decay > 0:
                cgood = False
              elif cinf == (len(cxfit) -1) or slope_d1 > 0 or slope_decay > 0:
                cgood = False  
              elif abs(slope_d1) == 0 or math.isnan(slope_d1) == True or tools.close(slope_d1, 0) == True or slope_d1 < slope_decay:
                release_rate.append(slope_decay)
                inf_slopes.append(slope_decay)
                  # print slope_decay, "=decay", slope_d1, "=d1"
                print max_line, min_line, x_inflection
                b = y_inflection - (x_inflection * 0)     #using x,y of inflection point, determine b of y=mx+b
                x_line = range(len(cxfit))
                for point in range(len(x_line)):
                  y_line = (0* point) + b   
                cytc_decay = x_inflection       #because slope is zero, x_inflection point denotes decay start
                # print slope_d1, "slope_d1", cytc_decay
                # plt.plot(cxfit, cyfit, color= 'green')
                # plt.axvline(x= x_inflection, color = 'black')
                # plt.plot(x_inflection, y_inflection, 'o', color= 'g')
                # plt.show()
              else:
                release_rate.append(slope_d1)
                slope_d1s.append(slope_d1)
                yy = []
                yymin = []
                bb = cyfit[cinf] - (cxfit[cinf] * slope_d1)
                xx = range(len(cxfit))
                for vals in range(len(xx)):
                  yys = (slope_d1 * vals) + bb
                  yy.append(yys)
                xxmin = np.arange(-500, 500, 5)
                for vals in range(len(xxmin)):
                  yymins = tools.sigmoid(vals, cpopt[0], cpopt[1], cpopt[2], cpopt[3], cpopt[4]) 
                  yymin.append(yymins)
                cmin_line = (sum(yymin[180:]) / len(yymin[180:]))  
                cmax_line= cmin_line + (2 * (cyfit[cinf] - cmin_line))
                cmax_line2 = sum(yymin[:20]) / len(yymin[:20])
                cytc_decay = tools.crosses(yy, cmax_line)
                # print "cslope", cslope, cytc_decay, 'cinf', cinf
                # plt.plot(xx, yy, color= 'red')
                # plt.xlim(xmax=65)
                plt.ylim(0, 10)
                # plt.plot(cxfit, cyfit, linewidth = 2)
                # plt.plot(cxfit[cinf], cyfit[cinf], 'o', color= 'red')
                # # #plt.scatter(xinf, yinf, color = 'red')
                # plt.axvline(x = c_init_ymxb, color = 'black')
                # #plt.axvline(x = bax_initiate, color = 'green')
                # plt.axhline(y= cmax_line, color = 'black')
                # plt.axhline(y = cmin_line, color = 'green', linewidth = 2)
                # plt.show()
            # print release_rate
            #print bax_file, id, "List", time_re_bax_initiation
        if bgood and cgood:             
          release_time = cytc_decay - bax_initiation
          print cytc_decay, bax_initiation, "values"
          time_re_bax_initiation.append(release_time)
          #print cytc_decay, bax_initiation, 'release intervals'
          plt.plot(cx, cy, cxfit, cyfit, color = 'green', linewidth = 2, label= "cytochrome c-GFP") 
          plt.plot(bx, by, bxfit, byfit, color = 'blue', linewidth=2, label = "mcherry-BAX")
          # plt.axvline(x = cytc_decay, color = 'green', linewidth= 1)  #Bax initiation
          # plt.axvline(x = bax_initiation, color = 'blue', linewidth=1)  # cytc Decay
          # plt.axhline(y=min(byfit), color = 'red')
          # plt.scatter(xinfl, yinfl, xfit_infl, yfit_infl, color= 'red')
          plt.axhline(y = 0, linewidth=3, color='black')
          plt.xlabel("Time (min)", size = 18)
          plt.ylabel("Fluorescence Intensity \n (Percent of Baseline)", size = 18 )
          plt.legend(loc=2, shadow=False, prop={'size':8}, bbox_transform=plt.gcf().transFigure)
          plt.xticks(size = 14)
          plt.yticks(size = 14)
          plt.ylim(ymax=5, ymin = -1)
          plt.xlim(xmax=100)
          plt.show()
          print bax_file, cytc_file, id, len(release_rate), len(bax_rate)
        print time_re_bax_initiation, 
        # print len(release_rate), len(bax_rate), len(inf_slopes), len(slope_d1s)
    # print 'regular=', min(slope_d1s), max(slope_d1s), 'infinity', len(inf_slopes), #min(inf_slopes), max(inf_slopes)
    # r = np.corrcoef(bax_rate, release_rate)[0, 1]
    avg_release = sum(release_rate) / len(release_rate)
    std = np.std(release_rate)
    # plt.hist(release_rate)
    print avg_release, std
    plt.plot(release_rate, bax_rate, 'o', color = 'green')
    plt.ylim(-5, 5)
    plt.show()
    slope, intercept, r_value, p_value, std_err = stats.linregress(release_rate, bax_rate)
    print r_value    

def plot(): #102114 protocol124
  baxxx = [26, 3, 24, 13, 3, 22, 25, 15, -7, 3, 14, 27, -6, -2, 0, 11, 17, 25, 29, 8, 9, 22, 24, 6, 30, 5, 23, 23, 6, 2, -3, 32, 15, -6, -2, 0, 5, 8, 8, 7, 18, 13, -3, 18, -2, 13, 0]
  bax = []
  bax1 = []
  for i in range(len(baxxx)):
    bax1.append(round(float(i)))
  for i in baxxx: 
    if i < 26 and i > -50:
      bax.append(i)

  avg = sum(bax) / len(bax)  
  std = np.std(bax)  
  print bax, "avg=", avg, "std", std, "n=", len(bax), len(baxxx), bax
  #param = stats.rayleigh.fit(bax, loc=0, scale= 150)
  bins = np.linspace(-25, 25, 20)
  #plt.plot(bins)
  events, edges, patches = plt.hist(bax, bins)
  #plt.hist(bax)
  #plt.ylim(ymax=75)
  #plt.xlim(xmin = -20, xmax=20)
  plt.xlabel ("Time (min)")
  plt.ylabel ("Number of Mitochondria")
  plt.axvline(x=0, color = 'black', linewidth=2)
  plt.axhline(y=0, color = 'black', linewidth=2)
  plt.show()

def box_whisker():    #baxKO= RGC on, DKO = RGC off
  data = []
  means = []
  # Hela = [0.23683173890718495, -2.24961206788583, -0.6263491716741517, -0.3339705376257015, 0.5077777288984477, -0.06601342769578361, -1.1257700261015249, 0.055206748464744715, -1.39506085515432, -2.807460567017833, -1.404623929276128, 1.0184621703139483, -1.6036862996871495, -1.404946094473967, -1.867038595265459, -0.638678197948608, -1.170418175827071, -0.5212585688822824, 0.05795480479428125, 1.2753269896593424, -0.4980504369105432, -1.9081241435765508, -0.07293106837113664, -1.5648547541691782, -0.4253450558249343, -1.4066791760862851, -1.9506887852754686, -1.8270113667594947, -2.5648257086151847, -1.8818940013407823, -3.2247937494661136, -2.9898218870025453, 0.47894108551453357, -1.5427435686895785, -2.1893485988930927, -1.3032414300378263, -0.26319088266799523, -1.411196869947902, -0.7212784161971041, -1.8315177074770894, -1.6972931715344806, -1.9516544566756953, -1.3199640189566304, -2.008068229259603, -0.015723558213525925, -2.6924454211114095, -0.4897828518325609, -1.3716859220355526, -0.807537309301328, -0.9199529346998983, -2.1304298536577275, -3.043155523561334, -1.995438131595078, -1.0275470761206407, -2.777297176235049, -1.809814937228777, -0.940730076774966, -0.6109274244531914, -1.7634253663411705, -1.6880586031394957, -0.08423581744628568, -0.8051132690425359, -1.1173509897468217, -3.636791528040705, -2.1709603553885315, 0.8647201725780853, -3.5548780594497953, -2.269710761004578, -1.0637952059577145, -2.0591825813444973, -0.09720536880384094, -1.8144524212044988, 0.2913193203157286, -0.9864259793949235, 0.6722275935395385, 0.29037382602724343, -0.9191276607768248, -1.1037373620235549, -2.167628103946756, -1.0753964497371407, -0.9046453833465229, -2.891181006966092, -1.2451320486487618, -3.5862317127780527, -0.73331056618498, -1.7048668330244148, 0.8827190095595048, 0.277477646550091, -2.9322529477511003, -0.5468602017241194, 0.38756535484523996, -1.3745896645309836, -0.8529625759469651, -0.9478777351272822, 1.3947166512724938, -0.05078306385590771, -1.120290749301348, 0.2084561021449008, -2.8403028566948407, -0.5012873524552545, -0.9047382187869704, -2.2642629670383334, -0.7260987345155563, -1.1835586040638055, -1.5830970351438507, -2.7455546347333155, -0.4993118933103362, -0.5405052528795764, -0.48637018317802583, -0.4831952399325376, -1.5794288021760505, -2.1839168486468368, 0.9017597066138929, -1.9613063630909544, -1.8991381737471835, 0.7933928778965106, -1.5699430551029134, -2.1258806660828755, -1.8164594241587655, -1.933178493850109, -1.4515734903662432, -1.9791998466594283, -2.587006812629882, 0.5538024454146823, -0.17204719784069863, 0.363908914745718, -0.8000094927405433, -1.048145251276714, -0.5688433808267676, -1.1310921635249862, 0.2776289956824792, -1.3762349695679983, -1.3545576491384321, -0.6952584570902056, 1.5850559642868625, -0.491546321881025, 0.49240960501489417, -0.05584823868111684, -1.616542724759839, 0.3951615122101898, -0.31170613168543554, 1.4509642874744175, 0.15948360269190434, 0.968930913311598, -0.7579997332879936, 0.8426916399949007, -2.08363642732748, -0.5110056378953983, -1.235898285697439, -1.1496938915886268, 0.8659235264084818, 0.008285175053814204, -0.43449097145065835, -0.8649247505991688, 1.5003533313825028, 1.475701175536639]
  # HCT = [-3.8188214040212536, -1.547372206680421, -2.381532291841083, -1.7233537158606513, -2.8643909641043215, -1.4461575131443414, -1.8461324517391915, -0.7583358949967546, -0.4568927653912538, -0.5649094982631293, -0.947271316499226, -1.8585138994035164, -1.4323855572696829, -1.8379118180733547, -2.1653654019310626, -1.218434894239759, -2.083174383315372, -1.5390578667286199, -0.8313347839995927, -2.506758315380836, -2.432400507621569, -0.800630831223478, -1.2513267299981063, -2.0978821430769794, -1.2733334704435335, -3.2066412888529205, -1.216995805647787, -1.3689944075121172, -0.9670134924181994, -2.9745308909919634, 1.4156871272632034, -1.2206501757370938, -0.6132176191788534, -1.911596319513689, -1.040256621954693, -1.1954864874782607, -1.449595454669271, -0.6300497182845524, -0.6404879358005556, -2.9706520506795426, -3.1543614242147653, 0.035734019154148026, -1.0213098797322848, 1.1375907768511204, -1.3314005777888713, -1.8233230155978164, -0.733714050338844, -2.2090804619438207, -2.6291833077170015, -1.6605660234301267, -1.8916828972856765, -0.6814230036044339, -1.3892721604938965, -1.1371228555489266, -0.6562956111052348, -1.7139565834492345, -2.102220196264516, -1.5391153142490257, -1.6848206339543201, 1.3678694662659665, -2.426603316514786, -2.2048953137582297, -2.400143795639514, 0.12941812071687264, -2.21210253455925, 1.7137236413050267, -1.4073024534410998, -0.8327129074058964, -1.6716766990009275, -0.5060025233708889, -1.344638006669225, -1.757816908992167, -0.9232569067934784, -1.899288001807817, -1.816550795923769, -0.3275460415315171, 0.3577167110880496, -4.665869553265852, -3.238331451943585, -3.0759351092046967, -0.7795279026699292, -0.9789417732937526, -1.980262091488971, 0.3186001263654846, -1.2403916772278565, -4.046922821383569, -3.589377347666131, -1.788977041895824, -0.8746030657696056, -1.0421364824028296, -0.949185588036039, -1.079686028882805, -1.7273387051085978, -1.7882180368625264, -1.5711564169644405, -0.7468739168372928, -0.9609281794710356, -1.122998463135572, -1.4118994057251988, -2.107093079016382, -1.4851967328303315, -1.7781726267597644, -1.5767611221610924, -0.9659214607749491, -1.7028310222789087, -2.5143947461319978, -1.1201072434429529, -2.533425624251427, -2.7780690216768424, -0.07060091091560634, -2.178302291768201, -0.8929172852142984, -2.2558236995513252, -1.9147237772241297, -3.6652434965205956, -2.3202275312438667, -2.9342626380276395, -0.7641160085001464, -2.352593072467158, -3.406968931185189, -2.6315466928519666, -2.74400971646063, -2.873467599808732, -3.174117995406509, -2.2495479928798323, -1.5894295612886935, -1.5157660737456347, -1.8289778811742585, -1.2773525608616991, -1.7418842434040147, -2.644989731263061, -0.6201125428054944, -0.9461274576828015, -1.262166242155328, -1.7572122161524526, -0.4351166146937627, -1.6659493349550731, -2.083483079016584, -2.3102239503468263, -0.9558995258583628, -2.455751285312725, -2.3933825423998125, -1.4840685820472292, -1.1167194709130825, -2.6647379107024656, -1.4424167980687836, -1.293317556716752, -2.561645204411296, -2.0349506055975395, -2.0911045528010535, -3.719777761199421, -1.8947454890650355, -1.621039862330471, -0.7116355943085179, -2.1118499422844312, -1.1201418623830135, -1.6138005637295594, -3.523249671977575, -0.9439799047414946, -0.789937080756306, 0.14906879817736757, -1.7166085413468077, -0.26798146476844226, 2.0550309792061627, -1.0417453425248484, -1.1393335727231468, -2.724843498820101, -3.695171462430091, -2.2507718623573907, -1.1265898418533506, 1.738275841380853, -1.3071701393006876, -5.355936628151209, -2.103689796339194, -2.523238089904965, -2.4767902969685034, -1.8808681088189727, -2.7690276733265033, -1.5006917551341228, -2.589862096534869, -2.720856606713678, -2.012993410516537, -1.9218589285919856, -3.1051557322190164, -2.5918419816863247, -1.5622767401041855, -1.9684136091803186, -2.4837788559199687, -1.6227012814860242, -0.12971836409155987, -2.9464942562410177, -0.037065876939440995, -2.149690464678456, -1.387118798620698, 0.3825107623351216, -2.7064161302417435, -0.7179204173727424, -1.1584859610320264, -0.7845227293085663, -1.5291113630580981, -3.3296430709733498, -3.926555122780028, -3.2782873848097505, -2.9403303134580994, -4.762076390116016, -3.7349158247169374, -3.0758372691695923, -3.032254032528961, -3.265160979552494, -0.6355540513899167, -0.720189642291861, -3.3398236370570658, 0.3664829053523651, -1.6229328891968497, -0.9809582459681639, -0.08325969288966044, -2.0222510448373354, -1.0348320443453978, -2.0733844480804073, -1.0146170091548927, -2.898921451438145, -2.632308494249264, -4.037510048369951, -2.2062750236350355, -1.229750390000674, -2.669512333228992, -2.7983124716591794, -0.958435518566647, -1.1296764994864776, 1.1895794117300167, -1.467476703724112, -1.223470379171596, -3.251101009228476, -2.2033843408291705, -1.0789491179377917, -3.7302318923834665, -1.9637033137300601, -1.8762547876337614, -1.7944341941248845, -1.695584380983178, -2.793974759855914, -1.9423158953539612, -0.06070587752748938, 0.5806863250181018, -2.8647873134849267, -1.7765537837591732, -2.1249557542938557, -2.5197231464621503, -1.5405516022456693, -3.065320366135996, -3.155939391899201, -2.3554044631798, -2.160632941638052, -2.123398800758548, -2.036590830949676, -3.0140447552890874, -4.130681695601219, -2.5111886077013126, -1.6192242588945722, -1.4336391758361215, -2.9766854550200295, -3.253923823189617, -2.224254003252688, -2.0510047500756725, -2.8506313930593983, -1.0480049804446503, -1.545497097096861, -3.1245073748143404, -1.7333934394630013, -2.16657695046465, -1.23497183392266, -1.1829207748659794, -2.092576069781622, -1.9545359809460512, -1.7928990151905577, -2.425259554136227, -1.1530066216429307, -0.950540824440458, -2.1644906674906172, -1.0220837619340226, 0.29387640664716735, -2.709593573964717, -0.29878780743343425, -0.9082210338511036, -1.2737792680715503, -3.671659803841352, -3.089391493895403, -1.2854967014540057, -1.9867326854290903, -0.7677994873434327, -1.3186422080253908, -2.5335226970077755, -1.8060814552979518, -2.3060558874815547, -1.7905879668829343, -2.7396939303819634, -2.912297298972527, -2.700565623596505, -3.3621093260923676, -2.7757351091914972, -3.870549979176692, -4.2316811860917225, -1.923911664700599, -1.73600503591019, -2.1962136449790144, -1.7277994526693716, -1.5374729357233166, -0.2891771195762975, 0.19626036154392706, -1.8558783426332495, -2.807944061542066, -2.119631931179658, -1.0814337687233166, -1.3481194779775685, -0.5670198863367263, -1.4042479415814195, -2.1783399473504192, 1.0740489534997812, -2.7876888477577118, -2.2266478147525275, -0.2116382824575959, -0.9594657894133763, -2.1592889959367265, -1.8169681081557274, -2.528771223422528, -4.247586888424918, -0.954090218320778, -0.9140225682614451, -1.501455140458419, -0.7660381542967035, -2.3508345015464953, -0.5783603055767442, -1.1645317313794277, -1.6542158664977062, -1.2354191116235034, -1.8555420923386186]
  # D407 = [-2.3368050326723107, -0.9604395429836351, -0.8159789358364736, -2.69001230565714, 0.042687351821541014, 0.5451164350184617, -1.5565223806230468, -1.309484004412145, -1.511683179065075, -2.2321542190153543, -2.221928303415502, -1.2704531088030135, -1.3254593184941053, -1.2592482742925946, -3.3023306659812564, -1.8986065392866247, -1.9544731772054236, -1.9559374916282262, -2.6078403425476204, -2.3192958606000618, -1.3141656055346387, -2.281937243160781, -1.9844824997578607, -2.7195645206604833, -2.9076967110323673, -2.375213621171911, -1.9104698887486165, -1.7682449013196286, -2.0830437244396847, -2.6356908932450103, -2.346741480006085, -2.5871669184650434, -2.807983699707957, -2.0641745333312755, -1.8363093813190035, -2.73832334743631, -2.156002576924819, -2.2985599283345213, -1.8519491154847683, -2.1286932512520957, -1.0478159463043033, -1.691033340219063, -1.0706207753701424, -1.2067596516965533, -1.394083605240578, -1.660184225578098, -1.7290581094304187, -2.505688073372543, -2.579427755776998, -1.7946145151544333, -1.5692099165747444, -2.24406248165822, -2.3907547978344517, -2.7898322507776654, -2.041021210385909, -1.4599702961842167, -3.1123166452529394, -2.108444844944421, -1.8527493367602643, -2.3700107567112463, -1.2352522347022565, -1.8734561096963827, -2.436954354453154, -2.151552060765871, -2.5395060688212245, -1.758148719717742, -1.9970129404987742, -2.31935634119232, -2.3239877900861767, -4.693389036695131, -1.9577408782594252, -0.9154351647793323, -1.180469674525513, -1.6593524215395448, -0.8612582336892513, -2.0887080255224206, -0.9954895746088241, -2.1402228992939447, -3.035043834732276, -2.2341475979131076, -1.9260394448799951, -0.5376449610649977, -1.3781699262388625, -0.9742578546124331, 2.0782499174863385, -0.588545653491503, -0.5480526151941648, 0.2137995870980992, -2.2082958503160564, -2.4352991319265618, -0.5960204339086504, -1.768723545640771, -1.4053698858685508, -1.3378131707475405, -2.426833519233202, -2.13468368733802, -2.830082622596404, -2.858176948665572, -0.8892086483943171, -0.055859088993665276, -1.5216545615613992, -1.5354431730041118, 1.2156110772223125, -0.9345778783402289, 0.19842482076642137, -0.6063884179453781, -1.4169309360793403, -1.2963443420640257, -0.9456962087260545, 2.324110184720176, -0.029730624627283417, -1.30276079691913, -0.9352039781720877, -3.2357145706008303, -1.338839685970891, -2.068916810519003, -0.8278083941650792, -1.9426879927098988, -2.981920606734751, -0.6880412693033467, -1.1015255052071844, -2.1041208495940027, -0.8622599986295313, -1.7250666892973874, -0.22284689765707866, -1.0240420233363776, -0.6574953114999447, -2.213812831619469, -2.200382148888963, -1.1510741144561571, 1.3214743468000463, -1.3759379022537301, -3.1604371939716893, -1.7491945425739317, -0.660057703416759, -1.8536639123683514, -0.9361948682480309, -0.8149582022823418, -1.2398435335783982, -1.8319919055754035, -2.6754712330397283, -0.8380030237531861, -2.4775765151773808, -2.0487131533422964, -1.493717415729145, -1.950355854109715, -1.6814958741804449, -2.102242799879415, -1.524149149423406, -0.9783479768261188, -3.542552037415868, -1.097327070680054, -1.3452904591015418, -1.5889376776303983, -1.8992730752490472, -1.6281126551294334, -1.8858511013348493, -2.05014656030429, -2.9336482960489527, -0.9410961285224363, -1.8722406072322413, -1.3864344072486663, 0.13795924351142058, -2.659058510186328, -1.318637137643352, -1.5736965794468343, -2.0322775988516284, -0.003494779335152044, -1.1998941360254272, -1.1055514961352757, -0.3526692190970694, -2.95645761378954, -1.2214534943756612, -1.8477509644029624, -1.0912327257736214, -2.3924019448536034, -1.110508570355993, -1.7815776339820553, -2.2752455784274304, -0.7359679978227096, -1.4633887995518795, -1.5053768660501667, -1.5615997286517334, -2.6596896611182976, -1.4630313864619136, -0.9125438315573098, 0.7013579900275129, -1.4592156145221886, -2.014315522795021, -1.7824918102126275, -0.8137789415666952, -1.6517835047470435, -1.2958746027495476, -1.498158957820845, -2.370383831116162, -1.8810289813839738, -1.5933165820915725, -1.1545113457240719, -1.7164476163881628, -1.0864150513769026, -1.1169662635954707, -0.8919841090497719, -0.9820899777671802, -1.5744348025296073, -1.540050713186303, -2.0697083135155614, -1.8035826140714928, -1.32767844741257, -1.0560818820658215, -3.306481138500573, -2.5311174557228617, -2.197785503625951, -3.5838120399066256, -2.386496495758965, -3.4478642873407983, -1.1277078140577483, -1.5951222561793688, -4.72697939072751, -4.075215465659392, -4.684724695767821, -4.893358890027522, -2.4770209103152565, -4.557541768602269, -2.4294311613128388, -2.775137887159885, -1.3791539111867788, -1.5855465266225874, -3.2786542482516383, -4.014568411831406, -2.431782219945827, -2.762031026526355, -4.006944481243173, -1.1619693770346842, -3.0537734055289993, -2.5514257473584747, -3.317213784915006, -2.5602358752785324, -1.263417179279086, -2.4949527600491797, -3.1424740526216417, -2.974911961651783, -2.498689666179835, -1.4395825842140544, -1.0986499845345947, -2.2965515285516633, -1.8473078130947937, -1.563842247108899, -1.758453835660446, -1.7888814259090078, -1.9805307182605532, -2.6264517040003796, -1.5975903718578577, -1.0507547888327555, -1.8476136867816424, -1.2436556611707303, -2.5655191516863702, -1.4283933491054577, -1.839780062091717, 1.6848876687248124, -1.561315054899277, -3.326608413091606, -1.6153760475387344, -2.053898818918296, -2.121715022496005, -2.478817874494937, -0.9301751052324592, -0.7752417111853597, -0.3092631380264327, -0.8896042970702827, -1.2082024640919902, -1.1623290990113329, -0.6218219395473572, -2.039190490760908, -2.5468681728210814, -0.8156521959141166, -1.0296058198878213, -1.9523883452926731, -1.1927326287764568, -1.2935056260325088, -1.5792866656258502, -0.8616138270244287, -1.4518528800904367, -1.121566018312378, -2.387006370110005, -1.3772883080315408, -1.5785341891989373, -2.261878536717288, -2.0292112794644845, -3.261688572325039, 0.31788215481798476, -1.512984388741324, -1.1010149818114594, -2.1183739908180588, -0.44439888689885765, -1.8266097801472514, -1.946285014497391, -2.6860985406678997, -0.723718242549591, 1.9080324181711317, -1.1197540623617719, -0.9796266635592886, -1.9631613293080237, -1.6575956551078175, -1.8287333983964893, -2.613409984624098, -1.8146121978220415, -2.484122851827109, -1.975060317246072, -1.6312112001012822, -2.7935939537000896, -1.4231162905853203, -1.8888368166804737, -3.6846079413838346, -2.5077470081645563, -3.1060506806358847, -1.5835635727645518, -2.6074694154941547, -2.346385864814118, -1.9814689929325746, -3.267304512528861, -1.3353843809856867, -1.8955811808941543, -1.490968029940252, -2.0654431042927333, -1.0673171550163947, -2.4375710570644302, -2.1367465571031232, -3.1166184961116916, -1.5153890938514478, -2.1854977351954723, -2.0794550241431793, -1.6675380218496105, -1.8718244027708437, -2.668972708154182, -1.8201735591042294, -1.3975414678221678, -1.9471288863954659, -1.5291730797570966, -1.969028885792285, -2.7419048446705605, -0.7701421041579677, -1.7834627855364256, -1.8656221471134329, -2.072386792784951, -1.4043352438573629, -0.8845689406489898, -1.716414992186671, -1.214784624408845, -2.734585067777144, 0.08298400045000646, -2.6884316793088345, -3.336936406800368, 0.2220718564499679, -1.6957136616373885, 1.580372255924936, -1.0183018582139691, -1.1974467848133779, -0.4569619713855802, -1.7375233722637895, -1.90961447484648, -0.8868774124172544, -2.270346136314556, -1.755658412159703, -1.1821170060859083, -1.314865899728904, -2.162369849372663, -1.2372456037788928, -1.1835889446880106, -1.7477619458292375, -1.4969456453128933, -0.8281734060833286, -3.0246171923331033, -0.8244719875094539, -1.7116863499559616, -2.1301648738042824, -1.3783563911010754, -1.3407424628349969, -0.7136730536197325, -1.9205348281350811, -3.526432799057221, -1.1720749928098755, -1.9754401658386174, -1.5865041911600917, -2.045563014316179, -2.456863305753631, -1.4407849723421853, -1.1184024940526534, -3.196279854984314, -1.692810117131118, -2.2800158721077235, -2.1536438820134274, -1.6437586974052874, -2.1248989032603594, -2.310854769039136, -0.764133357190666, -0.3954112564777015, -0.6236718076976913, -1.454348019009177, 0.14999606802665122, -0.02505937996241781, -0.11999764214924226, 0.8085177038486395, -1.5693803917280575, -1.121777704261443, -1.5777599385423127, -0.7313430649209386, -0.5266879351778891, -0.18349218972736092, -1.2526227754230634, -1.2121082062486968, -1.0724565275207443, -0.9330855383894446, -0.8988118502474275, -1.2229502165372923, -1.0319871926712945, -1.2571927730699841, -1.1244574478449074, -0.9029008841657841, -0.3404206710749027, -1.5552641306548862, -1.4444599324049952, -2.690933294990251, -0.18251039124792343, -3.366727698712963, -2.6310661368942605, -3.083223924917869, -2.7740478385159117, -3.5319942065454435, -4.766965101407633, -4.004488114031976, -2.3159321996830653, -3.071634495698617, -2.2667679227999833, -3.5362135810650326, -3.26379779492183, -4.128481412519842, -3.2146374708479346, -2.399212298842897, -3.329629323108388, -0.6755123419947352, -2.4701672727781148, -2.3411661320420474, -2.914831085947031, -2.6694909851506723, -3.370078020442042, -3.6049251823962933, -2.2122086557686917, -2.32445685528492, -2.1494279627130726, -3.4892041442336588, -2.3436784936210304, -1.4777892114858278, -3.7577571702828623, -2.3831213440428685, -1.0055602160589612, -2.156800518073221, -2.077265094804175, -2.207419028106878, -1.190907705799435, -1.6322950570875838, -1.620427989757875, -1.8438618508740834, -1.3873475748076232, -1.304192889942662, -2.044382400174537, -1.9931070028834998, -2.5037929183312544, -2.1243914973870477, -1.9018440165126878, -1.8891583478779919, -2.6147518697192287, -0.6848728023651155, -1.1575377387461372, -2.434635873501319, -1.9111781499302827, -2.4828813795254794, -1.0914484613104145, -2.8893375873980505, -2.295199562634246, -2.190855423485137, -1.7469262323110466, -1.134335170726995, -2.3280960380972773, -2.3109107417628163, -2.4320530723502123, -1.272000218663389, -2.0007102375112544, -1.6906025487199698, -2.640879894833283, -1.9967812079525138, -1.886578274760967, -0.9894604161501269, -2.8977990239317735, -2.3243736340982455, -2.193761026220385, -4.821803820822976, -4.712338890730119, -5.81502569906952, -4.1837524539332485, -3.3321130549547564, -2.431191086562477, -2.29964011393281, 0.8466138184011132, -1.675018652296381, -3.875303000567302, -2.213759076259226, -1.178394013285726, -4.083708867970206, -1.114104520625193, -3.9028017816180665, -2.0802109429558233, -4.060210749089278, -1.4473593777888454, -1.2604878066768606, -1.726151980593447, -1.3285181813433953, -1.6969888650046787, -2.0362431689050755, -0.9720628258211361, -0.343447207247471, -1.6643124570938184, -1.5568419270545537, -1.9313901554630863, -0.16181615107403896, -1.362526571025361, -0.8050555660326764, -2.161810588224371, -1.8347400261886375, -2.215368751527702, -1.0182971588885694, -0.8669327483770292, -2.343654617383461, -1.0235297268347776, -0.9320393982786485, -1.2021675760450503, -0.7566983479210081, -1.337232331688012, -0.73339145039397, -2.020867073668722, -2.701125906100277, -0.8881987404499442, -2.1784417758935897, -2.9012917841835772, -2.3280340093323515, -2.492172758811285, -2.1854156714198796, -2.2781096579318647, -3.4558192008918187, -2.9557265434039888, -2.6656698649893826, -2.6797060205323264, -1.8291860093458538, -0.9088026536919671, 0.032384791152583706, -1.128085698537346, -0.1417322421572041, -0.7546524547008262, -0.1216256806815664, -1.4012007492438072, -1.6651680641189337, -1.9136929450424842, 0.8568244170048627, -1.1488170967705849, 0.5482316085804425, -1.2187044764710608, -0.17733408698525832, -1.0161937213440737, -1.1269430126330942, -0.9453642579883785, 0.007148415481918298, -0.6454646297352322, -0.9947884776111233, -1.313697427610297, -1.8840768997668331, -1.1612045425905957, -0.7514468123833974, -1.6031509812535276, -1.9960704880488473, -1.3065350594922465, -1.6908970256323743, -2.066589588604974, -0.8669777713460479, -1.9086673834130594, -2.1202166412223655, -1.129498066675864, -1.432460192786826, -1.9152453782527992, -0.770952907553733, -1.8395135369586249, -1.9212975767650422, -1.1785842239540958, -1.8571033308509899, -2.026518370536327, -1.4692206081878019, -2.1964753642589665, -6.60979475339655, -3.9017068460108173, -3.260267299440124, -5.631552643184072, -1.596409029240322, -3.782531301693265, -1.7351279796081285, -2.664328539606227, -1.8363123830781474, -1.589885228043347, -1.6889788857700088, -2.131432509056274, -1.5675563219025532, -2.2540579356385897, -1.1688551819344752, -1.9089917440518833, -2.1564245411942706, -1.8019086024204594, -2.192025137912445, -0.9106942922095846, -0.5621359279540812, -2.938001230250796, -0.8069164203085253, -1.0271847476536116, -2.646691485846517, 0.45847113818843305, -2.103790777084908, -1.9351775642633078, -3.1340438825212384, -2.961683175145331, -1.5075487907917116, -2.3005241176445703, -1.6084248993298345, -4.014842423465852, -1.8276119558749746, -4.318125669965453, -3.089617608501897, -1.5570184254614783, -3.3590819163119097, -3.205930793853663, -2.897105930469105, -2.514963973902852, -2.627346112689678, -2.409345877816076, -2.7700379525501364, -2.1461903815442382, -2.0023201080903075, -2.392921831384185, -2.6660608785702777, -2.375736134340219, -2.5127823530189963, -3.2078863368566335, -1.518630936578613, -1.7648867762269416, -0.1558516447631984, -1.5075269430214913, -2.8507738629576176, -0.9716228648515516, -0.7201988576147417, -2.591954523802858, -3.847738104203156, -3.071960968560534, -2.2855726464507358, -1.8219998130446151, -2.781127919346007, -1.2974671568132716, -2.8510859464798366, -1.89742395937838, -1.3615812782077512, -1.84425951212813, -3.164638192636406, -0.9957617904107476, -2.4144863080765533, -1.2001801102381489, -1.9144789913840543, -2.3760915589944838, -2.0142396517938366, -1.3500388027493133, -2.9795188065731635, -2.443347105561319, -2.6467808949211813, -2.827553060210785, -1.1714390057888258, -1.8060247131601788, -2.525374156997649]
  # W = [-1.9686652855913873, -1.8154789917843022, -1.98454995488554, -1.629107295296061, -2.1288247989255424, -2.2377187545704627, -1.9463281704502948, -1.383322153161686, -1.4444332317651636, -1.860966512913813, -2.3225817641929876, -1.5431814822761678, -2.842297174953922, -2.698023616457717, -1.608587715824514, -3.023832743714741, -2.4256907453372123, -2.3713726064786145, -2.1369548165281813, -0.5619229984160562, -3.0605023824215376, -1.5811260159432632, -2.0836163932755154, -1.5239233543528532, -3.778715652149526, -2.954933856299178, -2.2167023722264605, -2.844716711645708, -2.0906724751223273, -1.9236142543879688, -3.389486965498626, 0.6199349459975027, -0.9509171224972204, -1.665730596302534, -0.4660946133788351, -0.012808342586384766, -0.9143080538385099, -0.6096826796447983, -1.123801303144753, -0.7252089604486273, -1.252299914007617, -0.21425475501727811, -1.368992193113575, -1.3952504515205233, -0.7631610269837155, -4.157735447813404, -4.235042643101414, -3.932714174758111, -4.203045547794048, -4.612770722063855, -5.088473706377181, -3.6734961606078, -3.5625863835352534, -5.709134011143131, -4.079460211049163, -4.9745775921562645, -4.772780030339673, -5.363374269354892, -6.193208250625197, -5.796002299094382, -4.54677797463596, -4.226502258485209, -3.530194030701992, -2.7246996778127826, -1.8913527868052458, -1.5075688868265575, -2.75535593449339, -3.3944283386005094, -3.6344623701606724, -3.371006147181204, -3.1493116329512865, -2.325503800234303, -1.0664378611117384, -2.61716976880512, -2.724234319387649, -1.7383894432076268, -1.7144068812218605, -1.6345948111785151, -2.459181560017098, -0.7820081770719883, -1.538271913630224, -1.3991509780604754, -2.4421214908206026, -2.1892001837023654, -1.7863266703147858, -2.392004733124951, -2.029576575896748, -1.2853691506380418, -4.962108779578446, -1.8733278808362213, -2.1364693370425116, -0.7940664805340873, -3.0522824833406785, -2.058385390471308, -3.2199173933388474, -0.23114480148071806, -1.3480301440805582, -2.8361814199254987, -2.2400901147807715, -3.9069496097795753, -2.502911520867621, -1.1655158794386011, -2.8648126762141235, -3.012997832022107, -2.5338888914343873, -0.6273003359716145, -2.5966438224535513, -4.032019944944817, -4.073180011035177, -2.4999582436965793, -0.55152290159562, -3.3565280954645096, -1.849913325136527, -0.37226403461094054, -3.8006804922385555, -1.1169683095900946, -0.6448640005629469, -1.6524411043765737, -3.4146267689376884, -0.7266485396873601, -2.207845996981209, -2.2026820687552333, -4.016550650810945, -1.0570385951424064, -2.905489585873088, -1.899427871845703, -3.00490523468512, -2.464086445045941, -2.258655980132902, -1.9027287500723848, -1.47384063523314, -4.020292451531605, -1.8539326730562797, -1.5959502432562973, -2.297131100091666, -2.1000786093981336, -1.7136918392713703, -0.9896812510328588, -2.2157596302015614, -0.4956644581949658, -0.6212780265221205, -1.302805497645281, 2.1671012572403923, -2.6535860000680427, -0.7168858713413603, -1.4621857150134008, -1.9485902687991508, -1.3773751334663404, -3.892976070480381, -0.18218483318707882, -2.6641117042024685, -0.771960539396183, -1.9048311311259674, -2.3769937559028946, -1.179486934022786, -2.487278841627734, -2.18149792910476, -0.6587179809247005, -0.0727194210105588, -2.1773560537290533, -0.9081554828012146, -0.8563575604190512, -1.7968357974538722, -3.333149516406519, -3.1132136205272634, -2.6783773068540007, -3.3033196961186913, -2.481217455924971, -2.539276330711286, -2.720277738295391, -2.6482715942597306, -3.0831730501853514, -2.501701560262148, -2.1029316930532076, -2.959713810636916, -2.509959496756643, -0.6280085244280542, -1.5881600191953573, -2.030281355065566, 1.2364824176157634, 0.024230765952289506, -0.8266046909567718, -3.4203820881050566, -3.9020034752505057, -2.948867627446577, -1.0858654321208432, -0.766647681031759, -1.2819815615053578, -2.8954540058347726, -2.3791250588422095, -2.83888022258717, -0.74921082510868, -2.4912932829409886, -3.1834936405724115, -2.710282729217103, -3.8256410579241718, -3.5515744514304135, -2.4059276293050944, -1.7065881302424444, -2.1504582354859156, -3.3456870904756157, -3.883577702316912, -1.7219891369828662, -1.7165786753501309, -1.8376919227335502, -3.3331793412555664, -2.7058564485541843, -2.103792196505259, -1.6863000071867555, -3.074061054511263, -3.5019172995348247, -2.9599305464264956, -0.7766098284483277, -3.3203427171793396, -2.800241093053808, -0.29445822995919646, -1.101558207366532, -0.6612876251520773, -2.8448058396742266, -2.6262440289157856, -3.4014573167113977, -0.7329851889240324, -1.043240636104507, -2.2300263720201716, -1.713062725292956, 0.29185304506513, -3.6521507895458005, 0.7679331450681481, -1.6140067967301903, -2.6944818544999056, -3.222242991237385, -1.9403108387901087, -2.091097923086954, -1.733215119844524, -2.75535593449339, -1.5075688868265575, -1.8189549899069049, -3.567943919967316, -3.371006147181204, -3.6344623701606724, -2.6352441799424193, -2.03949320171187, -1.1903126465976495, -2.932245489570365, -3.1227934189432838, -3.14533270467227, -1.6671111359312063, -1.6517080169781757, -0.08972334131964954, -1.8088163407622664, -2.192223030361638, -0.6151418501826251, -3.7958914311707352, -3.2732114369611005, -3.011961969062057, -3.208905275828847, -2.5388861504275684, -4.294199783641142, -3.7087428208274265, -2.51471937208703, -4.592782621820012, -2.171101320987485, -3.600432889150504, -3.036188465066602, -4.612980227593731, -3.586443824831545, -4.307881459105266, -2.609377327105497, -4.42300121946692, -2.1785868107437074, -4.225978572581011, -3.3863696047914105, -4.286990090319613, -2.802327215322535, -2.227374666867734, -1.879550906943192, -2.998110925930662, -2.9226215660967156, -1.392981604589859, -3.1445151253585895, -0.5963394892198449, -2.212400115533137, -3.172144423576613, -4.2549614868256045, -2.422914449637171, -4.209818708309015, -2.875215818367452, -3.4245949393257007, -2.4854155797618027, -2.0646458935769023, -3.41407347793412, -3.864818252472089, -3.4461607491885236, -2.848219839845086, -4.15281702016437, -3.5491423831614273, -3.140916275810168, -3.0150753574420324, -3.0367975456834553, -3.9446188367240267, -2.6419709638700377, -0.4779643899376654, -2.3788713197558407, -1.2788626625076005, -2.8557510704712348, -2.785952013207263, -3.5894931185756773, -2.986898023320157, -2.114927044060712, -2.3849658087284094, -4.575717472106449, -1.656089835319415, -2.352235822254303, -2.4751892452993443, -1.9202196854369513, -2.752371898004765, -2.5708309294138054, -1.4088517875999313, -1.0397606633918748, -3.0704904854012054, -3.093531537547857, -3.226864854029559, -2.3632220522903338, -2.525614643236145, -2.634517183992183, -3.1861790175298905, -4.276128117185845, -2.223270395274599, -1.0133264559403812, -1.6239033764804747, -1.3863787051177776, -1.6105366677274628, -3.8019467071097353, -1.8034582995937336, -3.1372912309631857, -3.047254152565993, -1.8899395747502201, 0.4374433247197172, -0.43013966991082636, -0.6130286706300428, -2.137333331333708, -0.4731630885619336, -1.7336351041921882, -3.897415514232658, -0.6504093449453295, -0.9640788772569088, -0.2663012059035544, -2.5045401395008815, -2.846942233351066]
  ONCRGC = [-3.608823606956928, -6.238829667855235, -7.3622305706434945, -6.722902767979996, -6.9119423354409655, -5.985460214708538, -4.945076766351701, -7.4461787954655145, -2.7046505560917504, -6.300228251719708, -5.376988788153347, -5.79443094230511, -5.824100661076135, -3.014655519174094, -6.505338519920784, -7.420864056974152, -7.004240610368588, -6.631964551551678, -5.663863853577057, -0.10913177220943085, -6.02396377228499, -6.195459218891133, -5.245617609986049, -5.065283010606309, -8.313312149201472, -4.957423987662304, -4.718977171916304, -3.9942194195656207, -5.213046222325124, -6.169077974667619, -5.140813211553333, -5.068429390101577, -5.375684221105025, -5.572946475786283, -2.818190185630893, -5.3518670048921875, -4.28529565554304, -4.8835284268314005, -4.2864532577805, -5.190237877788709, -1.7822857652109196, -5.313257961259824, -4.829418357896248, -4.341090764562417, -4.144037148097446, -5.10266950705007, -5.989220769357143, -6.33333883716791, -6.249605640253469, -4.590988287768243, -6.55240140980232, -6.749666252772515, -6.333321133288076, -6.106836555123853, -5.133933123533699, -6.479759243881705, -6.485576374160978, -6.488432735591233, -6.729126136499057, -3.4174342985465787, -6.15442219976914, -6.111067923109414, -3.091157869735088, -2.0628161730609613, -7.633526929724818, -5.357888444019209, -6.64133183440137, -6.497768638379569, -4.683445206789783, -3.540721089930127, -4.472655626714897, -4.298327308478214, -3.8227920443279784, -4.7372985529025735, -3.287866590575283, -3.592852228722821, -4.86746237781185, -4.650612891881721, -4.392871580258158, -4.305040179047842, -5.37140835273469, -4.7900898826018645, -4.941911104268652, -5.102903170448788, -4.79981354844973, -3.9574373154695976, -3.9730521679942847, -1.987244050628747, -4.647464903962698, -3.7528299584447127, -4.579720040301317, -2.423441544989694, -4.594765525019566, -4.918116168604272, -4.649333390033458, -3.938989851441319, -4.552370838656892, -5.948745784610636, -4.878177238514214, -6.449292568715395, -7.134291915371182, -5.961263840557565, -5.859688395673487, -4.059260591694929, -4.00608428723516, -5.454211324257659, -4.7780929728476815, -2.007038252026353, -6.37676542196874, -6.783838645210516, -4.303977857522952, -0.5628474276281665, -4.1311275692665514, -6.090359604045767, -4.611028171505511, -1.42568972978699, -6.033598980240168, -5.00565549532391, -4.7494596229513695, -6.636714971913578, -3.2035538147580302, -6.0845156879337745, -2.8624524924592443, -4.4869339764452105, -4.302147357201511, -6.175516375433367, -4.8885118250088215, -3.3267696750185474, -4.650937505058863, -2.6760291760011903, -3.765043419199229, -2.7439512103482997, -5.282183726308741, -5.968509107237634, -4.020153169257999, -6.440154621962112, -2.9774910342088092, -5.805621385570058, -6.696275586222908, -6.603956299414696, -6.813595126857791, -6.882915551085488, -6.670280139833947, -6.405129421994238, -4.726978340476291, -4.086497291743163, -4.98499970501301, -5.733841688294804, -2.5662672013814976, -6.956713418779291, -6.352574956040612, -6.726953358201356, -5.038662710662867, -3.37086913990936, -6.848756406402631, -6.031461454154335, -4.9123953161982445, -6.9733509202608275]
  STSRGC = [-5.897307940998905, -7.802783251470472, -6.458016719336295, -6.565516897928632, -6.705406104022834, -6.189715419135893, -6.522771721310046, -6.072566919236705, -6.18194879754123, -7.437519664518906, -7.00285455980093, -5.232566190566097, -5.658354635735111, -5.693550914573179, -5.453097277604324, -2.486999191801452, -5.3746529097196785, -1.5135688099229858, -5.579991431230711, -2.2959677977340776, -5.0547288712224985, -6.060374079733725, -6.172274234377859, -5.25804563106207, -7.241844983891885, -5.848681872357697, -5.945091997962311, -5.989677547406103, -4.629954765310014, -5.132943452594273, -6.184538942644696, -5.639161575381423, -6.451897098151542, -5.849378006307023, -6.099071642343534, -4.836083685930779, -3.974506674236215, -5.295401896147439, -6.784603486508097, -5.712502008930156, -5.773510746567029, -6.463550793964399, -5.802019182642341, -6.1794506322048885, -7.903149036786915, -6.507166881605478, -7.564428731547107, -7.265381909877207, -7.040610943408335, -7.760034721758207, -8.084641076924681, -6.817712861703599]
  BAXKO = [-4.333671729503716, -4.702043574063155, -5.553689015940044, -5.3087196040590205, -5.18803813500619, -3.510956647609866, -5.343106129875339, -4.9378591071429785, -4.04425819850874, -3.6921688695438144, -5.491945997322772, -4.84353830917061, -4.900764128046027, -2.517423057270856, -3.322644368877201, -4.540189748538708, -6.076096941296875, -5.1352701025310035, -4.456482248515596, -4.729773414632914, -4.615023147563886, -5.420132403771758]
  DKO = [ -7.447465796882185, -10.333587402361893, -7.140045607608002, -7.082019078062345, -4.738237060474038, -8.067835720571628, -10.03297335851225, -7.722130853006652, -6.532298695226293, -9.938776550189884, -7.613016594722381, -8.684551326633716, -24.619533706975492, -7.544586034470632, -6.9985387635838565, -9.449531499932004]
      # data.append(D407) 
  # data.append(HCT) 
  # data.append(Hela)
  # data.append(W)
  data.append(ONCRGC)
  data.append(STSRGC)
  data.append(BAXKO)
  data.append(DKO)
  labels = '1==RGCONC, 2=RGCSTS, 3=Retro_ON, 4=Retro_off'
  # avg_cyt = sum(D407) / len(D407)
  # avg_smac = sum(HCT) / len(HCT)
  # avg_aif = sum(Hela) / len(Hela)
  # avg_drp1 = sum(W) / len(W)
  avg_onc = sum(ONCRGC) / len(ONCRGC)
  avg_sts = sum(STSRGC) / len(STSRGC)
  avg_bax = sum(BAXKO) / len(BAXKO)
  avg_DKO = sum(DKO) / len(DKO)
  # means.append(avg_cyt)
  # means.append(avg_smac)
  # means.append(avg_aif)
  # means.append(avg_drp1)
  means.append(avg_onc)
  means.append(avg_sts)
  means.append(avg_bax)
  means.append(avg_DKO)
  # std_cyt = np.std(cytc)
  # std_smac = np.std(smac)
  # std_aif = np.std(aif)
  plt.boxplot(data, notch=False, sym='o', vert=True, whis=1.5,
        positions=None, widths=None, patch_artist=False,
        bootstrap=None, usermedians=None, conf_intervals=None)
  # plt.plot(means, [1, 2, 3, 4], 'k |')
  plt.xlabel(labels)
  plt.ylabel('RFU/min')
  plt.title("Rate of RecruitmentRGC conditions")
  plt.ylim(ymin = -15, ymax= 5)
  plt.xlim(0, 5)
  plt.show()
 
  # print avg_cyt, std_cyt, len(cytc), avg_smac, std_smac, len(smac), avg_aif, std_aif, len(aif)

def sigmoid_test():
 #'RGC ONC': [50, 0.79, 1, 2.11, 6.59], 'RGC + STS': [55, 0.46, 1, 10.71, 19.51], 
#'661W':[50, 0.46, 1, 2.27, 1.26], 'HCT WT':[50, 0.46, 1, 3.13, 1.2], 'D407':[50, 0.42, 1, 3.52, 1.41], 'HeLA':[48, 0.82, 1, 3.32, 1.4], #'HCT DKO': [50, 1.15, 1, 6.15, 2.93], 'HCT BAXKO':[50, 0.98, 1, 4.71, 1.82], 
  p0 ={'HeLA DRP1KO':[50, 1.88, 1, 5.15, 4.12], 'DRP1KO + WTDRP1':[50, 1.32, 1, 5.75, 2.43], 'WT + R247E':[50, 1.68, 1, 3.92, 2.26],'HeLA':[50, 0.82, 1, 3.32, 1.4] }
  # p0 = {'RGC ONC': [50, 0.23, 1, 0.94, 1.93],  'Retro_ON': [50, 0.26, 1, 1.17, 3.32], 'Retro_OFF': [64, 0.15, 1, 0.69, 2.08],'RGC + STS': [64, 0.11, 1, 0.81, 1.86]}
  xs = (range(100))
  ymax = []
  for key, val in p0.items():
    i = tools.sigmoid_inflection(xs, *val)
    ys = []
    for x in xs:
      y = tools.sigmoid(x, *val)
      ys.append(y)
      slope = tools.sigmoid_d1(xs[i], *val)
    avg = sum(ys[90:]) / len(ys[90:])
    ymax.append(avg) 
    print ymax
    plotters.plot_cell_start(xs, ys, xs, ys, slope, xs[i], ys[i], val[2])
    plt.plot(xs, ys, label=key, linewidth=2)
  legend = plt.legend(loc=2)
  frame = legend.get_frame()
  frame.set_linewidth(2)
  plt.rcParams['axes.linewidth'] = 2
  plt.ylim(0, ymax=5)
  plt.show()

  

if __name__ == "__main__":
  initiation_time()    

