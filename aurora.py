#!/usr/bin/env python 
# program aurora.py 

import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# read xls into a pandas dataframe:
#aurdatadf = pd.read_excel('annual_aurorae_1500_1650.xlsx', dtype={"Year": int, "Observations": int})
aurdatadf = pd.read_excel('data/annual_aurorae_nokorean_1500_1650.xlsx', dtype={"Year": int, "Observations": int})
#aurdatadf = pd.read_excel('data/annual_aurorae_korean_only_1500_1650.xlsx', dtype={"Year": int, "Observations": int})
#aurdatadf = pd.read_excel('occidental_annual_aurorae_1500_1650.xlsx', dtype={"Year": int, "Observations": int})
#aurdatadf = pd.read_excel('aurorae_1500_1650_Krivsky1996.xlsx', dtype={"Year": int, "Observations": int})
#aurdatadf = pd.read_csv('aurorae_KrivskyPejml.csv',  skipinitialspace=True)

#print(aurdatadf.head())

# convert dataframe to numpy array, discarding header line:
#aurdata=pd.DataFrame(aurdatadf.iloc[1:]).to_numpy
aurdata=pd.DataFrame(aurdatadf).to_numpy()

years=aurdata[:,0]
aurno=aurdata[:,1]
M=np.sum(aurdata)/(years[-1]-years[0])    # mean event rate
M_EMAP=np.sum(aurno[61:141])/82 

# 12221 smoothing:
aurnossm=(np.roll(aurno,2)+np.roll(aurno,-2)+2*np.roll(aurno,1)+2*np.roll(aurno,-1)+2*aurno)/8
aurnossm[0] = (2*aurno[0] + 2*aurno[1] + aurno[2]) / 5
aurnossm[1] = (2*aurno[0] + 2*aurno[1] + 2*aurno[2] + aurno[3]) / 7
aurnossm[-2] = (aurno[-4] + 2*aurno[-3] + 2*aurno[-2] + 2*aurno[-1]) / 7
aurnossm[-1] = (aurno[-3] + 2*aurno[-2] + 2*aurno[-1]) / 5

xstudy=aurnossm  

############## qa and qc ###########################

# array of local minima of xstudy:
mins_sim = years[1:-1][
    (xstudy[1:-1] <= xstudy[:-2]) &
    (xstudy[1:-1] < xstudy[2:])
]

mincyclength=6   # cycles shorter than this are artefacts to be removed
      # [empirical extrema-resolution threshold, calibrated by MC tests]
# activity on/off at phase...
aon=0.18
aoff=0.79

# Combine minima closer than mincyclength into one:
mins = []
i = 0
while i < len(mins_sim):
    group = [mins_sim[i]]
    i += 1
    while i < len(mins_sim) and mins_sim[i] - group[-1] < mincyclength:
        group.append(mins_sim[i])
        i += 1
    group = np.array(group)
    weights=xstudy[group-1500]
    weights[ weights == 0 ] = 0.1
    weights = 1/weights                 # weight minima by their depth, or ...
    #weights=np.where(xstudy[group-1500] == np.min(xstudy[group-1500]), 1, 0)  # use deepest 
    mins.append(int(np.round(np.sum(group * weights) / np.sum(weights))))
mins = np.array(mins)
# remove intervening maxima

clengths=mins[1:]-mins[:-1]

#### Regularization: shift minima btw. pair of very short - very long cycles by 1 yr:
i = 1
while i < len(mins)-2:
    if clengths[i-1] < 10 and clengths[i] > 12:
        mins[i]+=1
    if clengths[i-1] >12 and clengths[i] > 10:
        mins[i]-=1
    i += 1

print("Minima: ", mins)

#reduce st.dev. in long-short pairs:
# 12+6-> 11+7, 14+6 -> 13+7, 14+9 -> 13+10 :
mins_full= np.array([ 1511, 1519, 1530, 1541, 1548, 1560, 1567, 1577, 1587, 1596, 1609, 1621, 1634, 1644 ])

mins_krivsky=np.array([1510, 1523, 1536, 1543, 1550, 1559, 1566, 1577, 1588, 1596, 1609, 1618, 1626, 1640])
mins_nokorean=np.array([1512, 1523, 1539, 1550, 1559, 1566, 1577, 1587, 1596, 1609, 1619, 1627, 1640])
#or with strong regularization:
mins_nokorean=np.array([1512, 1523, 1540, 1550, 1559, 1566, 1577, 1587, 1597, 1609, 1619, 1627, 1640])


#mins=mins_krivsky

clengths=mins[1:]-mins[:-1]
print("clengths: ", clengths, np.mean(clengths[4:12]), np.std(clengths[4:12]))

qa=mins[:-1]
qb=qa+aon*clengths+0.5
qc=qa+0.5*clengths+0.5
qd=qa+aoff*clengths+0.5
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
ratesbc=np.empty(N,dtype=int)
ratesad=np.empty(N,dtype=int)

for i in range(N):
  yearai=np.where(years == mins[i])[0][0]
  yearbi=np.where(years == qb[i])[0][0]
  yearci=np.where(years == qc[i])[0][0]
  yeardi=np.where(years == qd[i])[0][0]
  yearaii=np.where(years == mins[i+1])[0][0]
  #print(yearai, aurno[yearai], "   ", yearbi, aurno[yearbi])
  eventsa[i]=np.sum(aurno[ yearai : yearbi ])
  eventsb[i]=np.sum(aurno[ yearbi : yearci ])
  eventsc[i]=np.sum(aurno[ yearci : yeardi ])
  eventsd[i]=np.sum(aurno[ yeardi : yearaii ])
  ratesa[i]=np.nan_to_num(np.sum(aurno[ yearai : yearbi ])/(yearbi-yearai))
  ratesb[i]=np.nan_to_num(np.sum(aurno[ yearbi : yearci ])/(yearci-yearbi))
  ratesc[i]=np.nan_to_num(np.sum(aurno[ yearci : yeardi ])/(yeardi-yearci))
  ratesd[i]=np.nan_to_num(np.sum(aurno[ yeardi : yearaii ])/(yearaii-yeardi))
  ratesbc[i]=np.nan_to_num(np.sum(aurno[ yearbi : yeardi ])/(yeardi-yearbi))
  ratesad[i]=np.nan_to_num((np.sum(aurno[ yearai : yearbi ])+np.sum(aurno[ yeardi : yearaii ]))/((yearbi-yearai)+(yeardi-yearci)))

fempfac = 1/(aoff-aon) - 1
with np.printoptions(precision=2):
    print("F= ", fempfac*(np.sum(eventsb)+np.sum(eventsc))/(np.sum(eventsa)+np.sum(eventsd)))
    print("Contrast 1560-1640: F= ", fempfac*(np.sum(eventsb[4:12])+np.sum(eventsc[4:12]))/(np.sum(eventsa[4:12])+np.sum(eventsd[4:12])))
    print("Mean rate 1560-1640: M= ", M_EMAP)
    

print("Cumulative events per quadrant: ", np.sum(eventsa[4:12]), np.sum(eventsb[4:12]), np.sum(eventsc[4:12]), np.sum(eventsd[4:12]))
qnorm=np.sum(eventsa[4:12]) + np.sum(eventsb[4:12]) + np.sum(eventsc[4:12]) + np.sum(eventsd[4:12])
qnorm *= 0.01
with np.printoptions(precision=4):
    print("Cumulative percentage per quadrant: ", round(np.sum(eventsa[4:12])/qnorm,1), 
    round(np.sum(eventsb[4:12])/qnorm,1), round(np.sum(eventsc[4:12])/qnorm,1), round(np.sum(eventsd[4:12])/qnorm,1))

print("Odd: ", eventsb[6]+eventsb[8]+eventsb[10], eventsc[6]+eventsc[8]+eventsc[10])
print("Even: ", eventsb[5]+eventsb[7]+eventsb[9], eventsc[5]+eventsc[7]+eventsc[9])
print("Odd2: ", eventsb[4]+eventsb[6]+eventsb[8]+eventsb[10], eventsc[4]+eventsc[6]+eventsc[8]+eventsc[10])
print("Even2: ", eventsb[3]+eventsb[5]+eventsb[7]+eventsb[9], eventsc[3]+eventsc[5]+eventsc[7]+eventsc[9])

print("Characteristics of individual cycles: ")

with np.printoptions(precision=1):
    print("Event rates: ", (eventsa+eventsb+eventsc+eventsd)/clengths)

print("Total events: ", (eventsa+eventsb+eventsc+eventsd))

with np.printoptions(precision=1):
    print("Femp: ", ratesbc/ratesad)
    print("F= ", fempfac*(eventsb+eventsc)/(eventsa+eventsd))


qall = np.dstack((qa,qb,qc,qd)).flatten()
eventall=np.dstack((eventsa,eventsb,eventsc,eventsd)).flatten()
rateall=np.dstack((ratesa,ratesb,ratesc,ratesd)).flatten()


plt.ion()	 # set interactive mode, so fig.is redrawn every draw() commanfig = plt.figure(1,figsize=(10,5))
fig = plt.figure(1,figsize=(12,4))
##plt.suptitle('$\lambda_R=$ ' + str(lambdar) + r' $\quad\tau=$ ' + str(tau))

maxevent=15
plt.xlim([1510,1650])
plt.ylim([0,maxevent])
plt.xlabel('Year')
plt.ylabel('Events / year')

# Minor ticks: mid-centuries, unlabelled
minor_ticks = np.arange(1530, 1650, 10)
plt.xticks(minor_ticks, minor=True)

#ax.set_axisbelow(True)

plt.vlines(mins,0,maxevent,color='#CCC',linewidth=25, zorder=-1.0)

#plt.vlines(mins_krivsky+0.25,0,maxevent,linestyle='--',color='#DDD', zorder=-1)

mins_FS=np.array([1511, 1522, 1532, 1539, 1549, 1558, 1567, 1575])
plt.vlines(mins_FS-0.3,0,maxevent,linestyle='-.',color='#7B7', zorder=-1)

mins_U=np.array([1510, 1525, 1533, 1542, 1552, 1565, 1574, 1584, 1595, 1609, 1620, 1632])
plt.vlines(mins_U+0.35,0,maxevent,linestyle='--',color='#7B7', zorder=-1)

mins_MCB=np.array([1517, 1531, 1548, 1558, 1566, 1576, 1584, 1597, 1606, 1621, 1627, 1637])
plt.vlines(mins_MCB+0.25,0,maxevent,linestyle=':',color='#A66', zorder=-1)

mins_Be=np.array([ 1509, 1530, 1543, 1550,  1559, 1573, 1586, 1597, 1608, 1620, 1630 ]) 
plt.vlines(mins_Be+0.25,0,maxevent,linestyle='--',color='#A66', zorder=-1)

mins_GSN=np.array([1620, 1631, 1640])
plt.vlines(mins_GSN+0.15,0,maxevent,linestyle='--',color='#99D', zorder=-1)

plt.vlines(mins,0,maxevent,color='k')

barwidths=((np.roll(qall,-1)-qall).astype(float))
barwidths[-1]=mins[-1]-qd[-1]
plt.bar(qall.astype(float),rateall.astype(float),width=barwidths,align='edge',color='#DA9')
#plt.bar(qall,eventall)

#plt.text(1513.5,75,'T-10',size='large')
for i in range(2,N):
    plt.text((qa[i]+qd[i])/2-0.75,maxevent-1.2,('T'+str(i-9)),size='large')\

plt.tight_layout()


input("Press [enter] to terminate.")

plotfilename="out.png"
#plotfilename = 'algebraic_' + str(tau) + '_' + str(lambdar) + '.png'
plt.savefig(plotfilename)

sys.exit("Bye!")



