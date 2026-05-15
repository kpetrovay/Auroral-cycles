#!/usr/bin/env python 
# program aursim.py 
import sys
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings('ignore')   # suppress warnings - use only after development finished!

rng = np.random.default_rng()   # random number generator

N=1000   # simulate N cycles
#M=10    # mean event rate/year
#F=3    # active/quiet event rate ratio
F = float(sys.argv[1])    # active/quiet event rate ratio
M = float(sys.argv[2])    # mean event rate/year
datafilename='mcsim_results.dat'

sigma=1.15    # st.dev. of observed cycle lengths (1.22 if integer years)
steps=10   # no. of time steps /year
mincyclength=6   # cycles shorter than this are artefacts to be removed
      # [empirical extrema-resolution threshold, calibrated by MC tests]

# activity on/off at phase...
aon=0.18
aoff=0.81

#cycno=arange(N)
clengths=rng.normal(11.0,sigma,N)   # cycle lengths, mean 11 yr, stdev 1.15 yr
# starting points of the 4 unequal quadrants in each cycle:
qa=np.array([np.sum(clengths[0:j]) for j in range(N)])  # cycle minima
qb=qa+aon*clengths                          # active phase onset times
qc=qa+0.5*clengths                          # cycle midpoints
qd=qa+aoff*clengths                         # quiet phase onset times

# array of time steps:
t=np.arange(np.sum(clengths)*steps)/steps

# define masks (arrays of 1s and 0s) for active/quiet phase;
activemasks=np.array([np.logical_and((t > qb[j]) , (t < qd[j])) for j in range(N)]).astype(int)
activemask=np.max(activemasks,axis=0)
quietmask=(np.logical_not(activemask)).astype(int)

# now generate the number of auroral storms each time interval
# assuming Poisson distrib. with overall mean 
# C*[F*(aoff-aon)+(1-aoff+aon)] = C*[1+(F-1)*(aoff-aon)]
# For C=1 the arbitrary rate is 1 for quiet phase and F for active phase
# C= M/arbmeanrate  where arbmeanrate is the rate for C=1
arbmeanrate= 1+(F-1)*(aoff-aon)
Cfac= M/arbmeanrate
eventrate=Cfac*(quietmask+F*activemask)/steps  # rate/step
aurno_step=rng.poisson(eventrate) 

##print(np.sum(aurno)/(t[-1]-t[0]))  # to check mean event rate/year ~ M
##print(aurno)

# now re-discretize to calendar years:
tint=t.astype(int)
yearno=int(t[-1])+1
years=np.arange(yearno)   # array of calendar years

# yearly auroral event totals:
aurno = np.zeros_like(years)
for i, y in enumerate(years):
    aurno[i] = np.sum(aurno_step[tint == y])

############ ANALYSIS ##################


# Apply smoothing:

aurnossm=(np.roll(aurno,2)+np.roll(aurno,-2)+2*np.roll(aurno,1)+2*np.roll(aurno,-1)+2*aurno)/8
aurnossm[0] = (2*aurno[0] + 2*aurno[1] + aurno[2]) / 5
aurnossm[1] = (2*aurno[0] + 2*aurno[1] + 2*aurno[2] + aurno[3]) / 7
aurnossm[-2] = (aurno[-4] + 2*aurno[-3] + 2*aurno[-2] + 2*aurno[-1]) / 7
aurnossm[-1] = (aurno[-3] + 2*aurno[-2] + 2*aurno[-1]) / 5

qa_int=qa.astype(int)
qb_int=qb.astype(int)
qc_int=qc.astype(int)
qd_int=qd.astype(int)

xstudy=aurnossm  # evaluate their performance one by one

############## qa and qc ###########################

# array of local minima of xstudy:
qa_sim = years[1:-1][
    (xstudy[1:-1] <= xstudy[:-2]) &
    (xstudy[1:-1] < xstudy[2:])
]

# Combine minima closer than mincyclength into one:
qa_emp = []
#qc_remove = []
i = 0
while i < len(qa_sim):
    group = [qa_sim[i]]
    i += 1
    while i < len(qa_sim) and qa_sim[i] - group[-1] < mincyclength:
        # maxima between the two minima should be removed
        #qc_remove.extend(qc_sim[(qc_sim > group[-1]) & (qc_sim < qa_sim[i])])
        group.append(qa_sim[i])
        i += 1
    group = np.array(group)
    weights=xstudy[group]
    weights[ weights == 0 ] = 0.1
    weights = 1/weights                 # weight minima by their depth, or ...
    #weights=np.where(xstudy[group-1500] == np.min(xstudy[group-1500]), 1, 0)  # use deepest 
    qa_emp.append(int(np.round(np.sum(group * weights) / np.sum(weights))))
qa_emp = np.array(qa_emp)
# remove intervening maxima

mins=qa_emp
clengths=mins[1:]-mins[:-1]

#### Regularization: shift minim btw. pair of very short - very long cycles by 1 yr:
i = 1
while i < len(qa_emp)-2:
    if clengths[i-1] < 9 and clengths[i] > 13:
        mins[i]+=1
    if clengths[i-1] >13 and clengths[i] > 9:
        mins[i]-=1
    i += 1

### alternative regularization:
#i = 1
#while i < len(qa_emp)-2:
#    if clengths[i-1]-clengths[i] < -3 :
#        mins[i]+=1
#    if clengths[i-1]-clengths[i] > 3 :
#        mins[i]-=1
#    i += 1


qa_emp=mins
clengths=mins[1:]-mins[:-1]
#print("clengths mean and st.dev.:", round(np.mean(clengths), 2), round(np.std(clengths), 2))

qa=mins[:-1]
qb=qa+aon*clengths
qc=qa+0.5*clengths
qd=qa+aoff*clengths
qb = qb.astype(int)
qc = qc.astype(int)
qd = qd.astype(int)

N=len(qa)
eventsa=np.empty(N,dtype=int)
eventsb=np.empty(N,dtype=int)
eventsc=np.empty(N,dtype=int)
eventsd=np.empty(N,dtype=int)
ratesa=np.empty(N,dtype=int)
ratesb=np.empty(N,dtype=int)
ratesc=np.empty(N,dtype=int)
ratesd=np.empty(N,dtype=int)

for i in range(N):
  yearai=np.where(years == mins[i])[0][0]
  yearbi=np.where(years == qb[i])[0][0]
  yearci=np.where(years == qc[i])[0][0]
  yeardi=np.where(years == qd[i])[0][0]
  yearaii=np.where(years == mins[i+1])[0][0]
  ##print(yearai, aurno[yearai], "   ", yearbi, aurno[yearbi])
  eventsa[i]=np.sum(aurno[ yearai : yearbi ])
  eventsb[i]=np.sum(aurno[ yearbi : yearci ])
  eventsc[i]=np.sum(aurno[ yearci : yeardi ])
  eventsd[i]=np.sum(aurno[ yeardi : yearaii ])
  ratesa[i]=np.nan_to_num(np.sum(aurno[ yearai : yearbi ])/(yearbi-yearai))
  ratesb[i]=np.nan_to_num(np.sum(aurno[ yearbi : yearci ])/(yearci-yearbi))
  ratesc[i]=np.nan_to_num(np.sum(aurno[ yearci : yeardi ])/(yeardi-yearci))
  ratesd[i]=np.nan_to_num(np.sum(aurno[ yeardi : yearaii ])/(yearaii-yeardi))

Femp=round((np.sum(eventsb)+np.sum(eventsc))/(np.sum(eventsa)+np.sum(eventsd)), 2)
#with np.printoptions(precision=2):
    #print("F= ", Femp)

epct=np.sum(eventsa+eventsb+eventsc+eventsd)/100
qpcta=round(np.sum(eventsa)/epct, 2)
qpctb=round(np.sum(eventsb)/epct, 2)
qpctc=round(np.sum(eventsc)/epct, 2)
qpctd=round(np.sum(eventsd)/epct, 2)
#print("Cumulative percents per quadrant: ", qpcta, qpctb, qpctc, qpctd)

#print("Scatter in event rates: ", round(np.std((eventsa+eventsb+eventsc+eventsd)/clengths), 2))

#print("Scatter in total events: ", round(np.std(eventsa+eventsb+eventsc+eventsd), 2))

Farr=((eventsb)+(eventsc))/(eventsa+eventsd)
Farr=Farr[np.isfinite(Farr)]
#with np.printoptions(precision=1):
    #print("Scatter in F: ", round(np.std( Farr ), 2) )


qall = np.dstack((qa,qb,qc,qd)).flatten()
eventall=np.dstack((eventsa,eventsb,eventsc,eventsd)).flatten()
rateall=np.dstack((ratesa,ratesb,ratesc,ratesd)).flatten()





############ Evaluation ################

def nearest_distances(a, b):
    """For each element of a, return distance to closest element of b."""
    dists = np.zeros(len(a), dtype=int)
    for i, x in enumerate(a):
        dists[i] = np.min(np.abs(b - x))
    return dists

#print("Minima: ", len(qa_emp), len(qa_int))

# qa_emp -> qa_int
dist_emp_to_int = nearest_distances(qa_emp, qa_int)
#print("qa_emp compared #with qa_int:")
rate1=round(100 * np.mean(dist_emp_to_int <= 1),1)
#print("  % within <=1 year:", rate1)
rate2=round(100 * np.mean(dist_emp_to_int <= 2),1)
#print("  % within <=2 year:", rate2)

# qa_int -> qa_emp
dist_int_to_emp = nearest_distances(qa_int, qa_emp)
#print("qa_int compared #with qa_emp:")
#print("  % within <=1 year:", 100 * np.mean(dist_int_to_emp <= 1))
#print("  % within <=2 year:", 100 * np.mean(dist_int_to_emp <= 2))

#print("\nM, F, rate1, rate2, Femp:")
#with np.printoptions(precision=1):
   #print(M, F, rate1, rate2, Femp) 

#print sys.argv[1], sys.argv[2], sys.argv[3], revBpol, revWSOB, revdipmom, halfmaxnpolarb, halfmaxWSOB
printout =  sys.argv[1] + '\t' + sys.argv[2] + '\t' + str(rate1) + '\t' + str(rate2) + '\t' + str(Femp) + '\t'  + str(round(np.mean(clengths), 2)) + '\t' + str(round(np.std(clengths), 2)) + '\n'
with open(datafilename, 'a') as outtable:
    outtable.write(printout)
outtable.closed

#print(printout)

#sys.exit("Bye!")

