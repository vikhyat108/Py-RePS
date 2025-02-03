import numpy as np
import matplotlib.pyplot as plt
import Pk_library as PKL
import readgadget
from matplotlib import gridspec

fb = 0.0490/0.3175
fc = 0.2685/0.3175

####################		PLOT-1 T(k)

file = open("./Tk_99_30.txt", "r")		
x1 = np.genfromtxt("./Tk_99_30.txt", usecols=0)
tk_b = np.genfromtxt("./Tk_99_30.txt", usecols=2)
tk_c = np.genfromtxt("./Tk_99_30.txt", usecols=1)
tk_n = np.genfromtxt("./Tk_99_30.txt", usecols=5)


file = open("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_rescaled_transfer_z99.0000.txt", "r")		
x2 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_rescaled_transfer_z99.0000.txt", usecols=0)
db_0 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_rescaled_transfer_z99.0000.txt", usecols=2)
dc_0 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_rescaled_transfer_z99.0000.txt", usecols=1)
dn_0 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_rescaled_transfer_z99.0000.txt", usecols=5)

ratio_b = np.interp(x2,x1,tk_b)/db_0 - 1
ratio_c = np.interp(x2,x1,tk_c)/dc_0 - 1
ratio_n = np.interp(x2,x1,tk_n)/dn_0 - 1

plt.rcParams["axes.linewidth"] = 2.0
plt.rcParams["xtick.major.size"] = 10
plt.rcParams["xtick.minor.size"] = 5
plt.rcParams["ytick.major.size"] = 10
plt.rcParams["ytick.minor.size"] = 5
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["legend.frameon"] = 'False'
plt.rcParams['figure.figsize'] = [20, 7]
plt.rcParams['font.family']="serif"
plt.rc("text", usetex=True)
plt.rc("font", size=30)

plt.subplot(1, 3, 1)

plt.plot(x2,ratio_c,color='r',label = 'ratio_c')
plt.plot(x2,ratio_b,color='g',label = 'ratio_b')
plt.plot(x2,ratio_n,color='b',label = 'ratio_n')
plt.axhline(y=-0.001,ls = '--')
plt.axhline(y=0.001,ls = '--')
plt.xscale('log')
plt.ylim(-0.005,0.005)
plt.xlim(1e-5,10)
plt.xlabel('k (h/Mpc)',fontsize=25)
plt.ylabel('Ratio - 1',fontsize=30)
plt.title('T(k) for 2F(M_n=0.3 eV)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize = 20)




####################		PLOT-2 G.R.

file = open("./fb_30.txt", "r")		
x1 = np.genfromtxt("./fb_30.txt", usecols=0)
y1 = np.genfromtxt("./fb_30.txt", usecols=1)

file = open("./fc_30.txt", "r")		
y2 = np.genfromtxt("./fc_30.txt", usecols=1)

file = open("./fn_30.txt", "r")		
yn = np.genfromtxt("./fn_30.txt", usecols=1)

x2 = x1

#reps 00
file = open("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fb_z99.0000.txt", "r")		
x3 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fb_z99.0000.txt", usecols=0)
y3 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fb_z99.0000.txt", usecols=1)


#reps 99
file = open("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fc_z99.0000.txt", "r")		
x4 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fc_z99.0000.txt", usecols=0)
y4 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fc_z99.0000.txt", usecols=1)


#reps 99
file = open("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fn_z99.0000.txt", "r")		
x5 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fn_z99.0000.txt", usecols=0)
y5 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_fn_z99.0000.txt", usecols=1)


y6 = np.interp(x5,x1,yn)/y5 -1
y7 = np.interp(x3,x1,y1)/y3 -1
y8 = np.interp(x4,x1,y2)/y4 -1


plt.subplot(1, 3, 2)

plt.plot(x4,y8,color='g',label = 'ratio_c')
plt.plot(x3,y7,color='r',label = 'ratio_b')
plt.plot(x5,y6,color='b',label = 'ratio_n')
plt.ylim(-0.005,0.005)
plt.xscale('log')
plt.xlim(1e-5,10)
plt.axhline(y=-0.001,ls = '--')
plt.axhline(y=0.001,ls = '--')
plt.xlabel('k (h/Mpc)',fontsize=25)
#plt.ylabel('Ratio - 1',fontsize=20)
plt.title('Growth Rates for 2F(M_n=0.3 eV)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize = 20)





file = open("/home/vikhyat/Desktop/CLASS_Jan_2020/class_public/output/explanatory_0.3000_pk.dat", "r")
k = np.genfromtxt("/home/vikhyat/Desktop/CLASS_Jan_2020/class_public/output/explanatory_0.3000_pk.dat", usecols=0)
Pm_0 = np.genfromtxt("/home/vikhyat/Desktop/CLASS_Jan_2020/class_public/output/explanatory_0.3000_pk.dat", usecols=1)


file = open("/home/vikhyat/Desktop/CLASS_Jan_2020/class_public/output/explanatory_0.3001_pk.dat", "r")
k = np.genfromtxt("/home/vikhyat/Desktop/CLASS_Jan_2020/class_public/output/explanatory_0.3001_pk.dat", usecols=0)
Pm_99_cl = np.genfromtxt("/home/vikhyat/Desktop/CLASS_Jan_2020/class_public/output/explanatory_0.3001_pk.dat", usecols=1)


file = open("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_Pm_rescaled_z99.0000.txt", "r")		
x6 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_Pm_rescaled_z99.0000.txt", usecols=0)
y6 = np.genfromtxt("/home/vikhyat/Desktop/reps-jan17_check/output_030/mnu000_Pm_rescaled_z99.0000.txt", usecols=1)

file = open("./Pm_99_30.txt", "r")		
x1 = np.genfromtxt("./Pm_99_30.txt", usecols=0)
pm_99 = np.genfromtxt("./Pm_99_30.txt", usecols=1)


ratio_r = pm_99/np.interp(x1,x6,y6) - 1
ratio_cl = np.interp(x1,k,Pm_99_cl)/pm_99 - 1
ratio = np.interp(x1,k,Pm_99_cl)/np.interp(x1,x6,y6) - 1

plt.subplot(1, 3, 3)

plt.plot(x1,ratio_r,color='r',label = 'ratio_m_r')
#plt.plot(x1,ratio_cl,color='b',label = 'ratio_m_cl')
#plt.plot(x1,ratio,color='k',label = 'ratio_cl_r')
plt.ylim(-0.005,0.005)
plt.xlim(1e-5,10)
plt.axhline(y=0.001,ls = '--')
plt.axhline(y=-0.001,ls = '--')
plt.xscale('log')
plt.xlabel('k (h/Mpc)',fontsize=25)
#plt.ylabel('Ratio - 1',fontsize=20)
plt.title('P(k) for 2F(M_n=0.3 eV)',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize = 20)




plt.savefig('2F_chatgpt')
plt.show()
plt.close()






