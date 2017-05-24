import numpy as np
from scipy.spatial import distance
from scipy.optimize import curve_fit
import csv 
import math

def get_rows(fname, nheader):
  counter = 0
  rows = []
  with open(fname, 'rU') as f:
    csvrows = csv.reader(f, delimiter=',',quotechar='|')
    for row in csvrows:
      if counter > nheader - 1:
        rows.append(row)
      counter += 1
  return rows

def get_columns(fname, nheader):  
  counter = 0
  rows = []
  columns = []
  with open(fname, 'rU') as f:
    csvrows = csv.reader(f, delimiter=',',quotechar='|')
    for row in csvrows:
      if counter > nheader - 1:
        rows.append(row)
      counter += 1  
  headers = rows[0]

  #column_1 = [row[0] for row in rows[1:]]
  #i = 0
  #for column in rows: 
  #  column[i] = [row[i] for row in rows[1:25]]
  #  i += 1
  #  columns.append[column[i]]
  #  if i + i[1] == i:


  #print column[0]  
  print rows[0] #columns    


def get_val_time(rows, ids):
  val_time = {}
  for id in ids: 
    val_time[id] = []
  for row in rows:
    if len(row[5]) == 0:
      continue
    id_str = row[5][-4:]
    id = int(id_str)
    if id in ids:
      val_time[id].append([float(row[0]), int(row[4])])
  return val_time    

def get_norm_const(timeseries):  #First value in timeseries
    norm_const = 0
    i = 0
    while norm_const <= 0:
      norm_const = timeseries[i][0]
      i += 1
    return norm_const

def log_val_time(val_time):
  for id, timeseries in val_time.items():
    val_time[id] = [[math.log10(float(entry[0])), entry[1]] for entry in timeseries]

def normalize_val_time(val_time):  #normalize data to 1 for comparison
  for id, timeseries in val_time.items():
    norm_const = get_norm_const(timeseries)
    val_time[id] = [[entry[0] / norm_const, entry[1]] for entry in timeseries]

def filter_1(val_time): #Number of data points necessary for curve fit function
  for id, timeseries in val_time.items():
    if len(val_time[id]) < 5:
      del val_time[id]
    else: 
      continue 
  return val_time      

def filter_2(xfit, yfit): #Remove data that lacks plateau at min/max thresholds  
  diff = yfit[-10] - yfit[-1]
  print diff, "difffff"
  if diff < 0.025:
    return None


def get_max_bax(timeseries):
    norm_const_max = 0
    i = 0
    while norm_const_max <= 0:
      norm_const_max = timeseries[-1][0]
      print norm_const_max
      i += 1
    return norm_const_max

def cytc_norm_baseline(timeseries):
    norm_const_cytc = 0
    i = 0
    while norm_const_cytc <= 0: 
      norm_const_cytc_time = timeseries[i][-1]
      i += 1
    return norm_const_cytc

def normalize_max_bax(val_time):
  for id, timeseries in val_time.items():
    norm_const_max = get_max_bax(timeseries)
    val_time[id] = [[entry[0] / norm_const_max, entry[1]] for entry in timeseries]

def sigmoid(x, x0, k, y0, a, b):
  y = a / (b + np.exp(-k*(x-x0))) + y0
  return y

def sigmoid_d1(x, x0, k, y0, a, b):
  y = (a * k*np.exp(k*(x+x0))) / (b*np.exp(k*x) + np.exp(k*x0))**2
  return y

def sigmoid_d2(x, x0, k, y0, a, b):
  y = (a * k**2 * np.exp(k* (x + x0)) * ((np.exp(k*x0))-(b*np.exp(k*x)))) / (b*np.exp(k*x) + np.exp(k*x0))**3
  return y  

def linear(x, m, b):
  y = (m*x) + b
  return y

def gaussian(x, a, b, c):
  y = a * np.exp(-(((x-b)**2) / (2 * c**2)))
  return y

def gauss_norm(x, a, b, c):
  y = (a / (c*np.sqrt(2*np.pi))) * np.exp(-(((x-b)**2) / (2 * c**2)))

def get_sigmoid_fit(x, y, verbose=False):
  p0 = [35, 0.1, 1.0, 5.0, 1.0]
  popt, pcov = curve_fit(sigmoid, x, y, p0)
  xfit = np.linspace(x[0], x[-1], 500)
  yfit = sigmoid(xfit, *popt)
  if verbose: 
    print "y = {3} / ({4} + exp(-{1}*(x-{0}))) + {2}".format(*popt)
  return xfit, yfit, popt, pcov 

def get_fit(x, y, p0, function):
  popt, pcov = curve_fit(function, x, y, p0)
  xfit = np.linspace(x[0], x[-1], 500)
  yfit = function(xfit, *popt)
  return xfit, yfit, popt, pcov   

def bax_fit(data, id):
  bx, by = get_xy(data)
  try:
    p0 = [30,  0.1, 0.5, 3.0, 1.0] ##For TC cells
    xfit, yfit, popt, pcov = get_fit(bx, by, p0, sigmoid)
    if close(yfit[0], yfit[-1]):
      try: 
        p0 = [60, 0.2, 1.0, 6.0, 1.0]
        xfit, yfit, popt, pcov = get_fit(bx, by, p0, sigmoid)
        if close(yfit[0], yfit[-1]):
            xfit, yfit, popt, pcov = None, None, None, None 
            print "No Fit" 
      except RuntimeError or RuntimeWarning:
        xfit, yfit, popt, pcov = None, None, None, None 
        print "No Fit"
  except RuntimeError or RuntimeWarning:
    xfit, yfit, popt, pcov = None, None, None, None
    print "No Fit"  
  if xfit == None: 
    try:     
      p0 = [50, 0.25, 1.0, 1.0, 4.0] ## For explants
      xfit, yfit, popt, pcov = get_fit(bx, by, p0, sigmoid)
      print "Explant"
      if close(yfit[0], yfit[-1]):
        xfit, yfit, popt, pcov = None, None, None, None 
        print "No Fit" 
      if xfit == None:
        p0 = [100, 0.25, 1.0, 1.0, 4.0]  
        xfit, yfit, popt, pcov = get_fit(bx, by, p0, sigmoid)
        print "Explant"
        if close(yfit[0], yfit[-1]):
          xfit, yfit, popt, pcov = None, None, None, None 
          print "No Fit" 
    except RuntimeError or RuntimeWarning:
      xfit, yfit, popt, pcov = None, None, None, None
      print "No Fit"  
  return xfit, yfit, popt, pcov

def cytc_fit(data, id):
  cx, cy = get_xy(data)
  try:  
    p0 = [60, 0.05, 1.0, 3.0, 1.0]
    cxfit, cyfit, cpopt, _ = get_fit(cx, cy, p0, sigmoid)
    if close(cyfit[0], cyfit[-1]):
      try: 
        p0 = [60, 0.5, 0.7, 13.0, 1.0]
        cxfit, cyfit, cpopt, _ = get_fit(cx, cy, p0, sigmoid)
        if close(cyfit[0], cyfit[-1]):
          try: 
            p0 = [30,  0.1, 0.5, 3.0, 1.0]
            cxfit, cyfit, cpopt, _ = get_fit(cx, cy, p0, sigmoid)
            if close(cyfit[0], cyfit[-1]):
              cxfit, cyfit, popt = None, None, None
              print "No cytc Fit"
          except RuntimeError or RuntimeWarning:
            cxfit, cyfit, cpopt = None, None, None
            print "No cytc Fit"   
      except RuntimeError or RuntimeWarning:
        cxfit, cyfit, cpopt = None, None, None
        print "No cytc Fit"         
  except RuntimeError or RuntimeWarning:
    cxfit, cyfit, cpopt = None, None, None
  if cxfit == None: 
    try: 
      p0= [100, -0.1, -1, 3, 1]
      cxfit, cyfit, cpopt, _ = get_fit(cx, cy, p0, sigmoid)
      if close(cyfit[0], cyfit[-1]):
        cxfit, cyfit, popt = None, None, None
        print "No cytc Fit"
    except RuntimeError or RuntimeWarning:
      cxfit, cyfit, cpopt = None, None, None
  return cxfit, cyfit, cpopt


def get_lagtime(rows, names):
  #fname = "../bax_agg/Image Acquisition Times.csv"
  #dtype = {'names': ('name', 'time'), 'formats': ('S40', 'i4')}
  #lagtimes = np.loadtxt(fname, dtype=dtype, skiprows=4, delimiter=',', usecols=(0, 2))
  lag_name = [] 
  lag_time = []
  for row in rows: 
    lag_name.append(str(row[0]))
    lag_time.append(int(row[2]))
  return lag_name, lag_time    

  #fname = "../bax_agg/Image Acquisition Times.csv"
  #dtype = {'names': ('name', 'time'), 'formats': ('S40', 'i4')}
  #lagtimes = np.loadtxt(fname, dtype=dtype, skiprows=4, delimiter=',', usecols=(0, 2))

# def get_image_lagtime(names):  #time that needs to be added to start of acquisition
#   fname = "../bax_agg/Image Acquisition Times.csv"
#   rows = get_rows(fname, 4)
#   lagtime = {}
#   for names in rows:  
#     print type(names[0]), type(names[2])
#     lagtime[names[0]] = float(names[2])
#   #print lagtime
#   return lagtime    


def sigmoid_inflection(xs, x0, k, y0, a, b):  #inflection point for xfits 
  ys = [sigmoid_d2(x, x0, k, y0, a, b) for x in xs]
  return crosses(ys,0)

def inflection_index(xs, function, parameters):
    ys = [function(x, *parameters) for x in xs]  
    return crosses(ys, 0)

def orig_inflection(xfitinfl, xs, timeseries): #convert xfit inflection into raw data inflection x value, then to corresponding timepoint.
    raw_xinfl = (float(xfitinfl) / 500) * (len(xs))
    for entry in timeseries.items():
      print "entry1", entry[1]
      if entry[1] == 0:
        time_infl = raw_xinfl
      if entry[1] > 0:
         time_infl = raw_xinfl + entry[0]
    return time_infl


def crosses(f, min_line):
  #print "inside", min_line
  for i in range(len(f)):
    if i == (len(f)-1): 
      return i
    y1 = f[i]
    y2 = f[i+1]
    #print y1, y2, min_line
    above = y1 > min_line and y2 > min_line
    below = y1 < min_line and y2 < min_line
    if not (above or below):
      #printprint y1, y2, i, min_line
      return i

def get_start_time(xs, k, xinfl, yinfl, min_line):
  ys = [(k*x + (yinfl - k*xinfl)) for x in xs]
  print "outside", min_line, crosses(ys, min_line)
  return crosses(ys, min_line)

def start_time_fivepercent(xs, inf, slope):
  thr = inf - (2.0 / slope)
  sat = inf + (2.0 / slope)
  return thr, sat

def close(x, y, eps=1e-7):
    return abs(x - y) < eps

