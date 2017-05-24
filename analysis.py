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

 

if __name__ == "__main__":
  initiation_time()    

