import numpy as np
import math
import matplotlib.pyplot as plt
import timeit
import time
from pathlib import Path
from scipy.integrate import ode
from scipy.integrate import solve_ivp
import classy
from classy import Class
import subprocess
import os

##########################################################
##########	PATH TO CLASS FOLDER	##################
##########################################################
Class = Path("../CLASS")








start_time = timeit.default_timer()
h = 0.6711
w0 = -1.
wa = 0.
M_nu = 0.30
On0 = M_nu/(93.14*h*h)
Ob0 = 0.0490
Oc0 = 0.2613483
Ocb0 = Oc0 + Ob0
Or0 = ( 2.469e-5 )/(h*h)
Neff = 0.00641
OR0 = ( ( Neff*(7./8.)*pow(4./11.,4./3.) )+1.)*(Or0)

Tcmb_0 = 2.7255
N_nu = 3.0
kb = 8.617342e-5
wrong_nu = 1
Gamma_nu = 0.71611

file = open("./FF_GG/FF_GG.txt", "r")
y1 = np.genfromtxt("./FF_GG/FF_GG.txt", usecols=0)
FF1 = np.genfromtxt("./FF_GG/FF_GG.txt", usecols=1)
GG1 = np.genfromtxt("./FF_GG/FF_GG.txt", usecols=2)

if M_nu != 0.:
	if wrong_nu == 2:
		On0 = M_nu/(93.14*h*h)
	else:
		A = 1.0
		Gnu4 = Gamma_nu**4
		a4 = A**4
		pi4 = (math.pi)**4
		y = (M_nu*A)/(Gamma_nu*kb*N_nu*Tcmb_0)
		On1 = ((15.*Gnu4*N_nu*Or0/(a4*pi4))*np.interp(y,y1,FF1))
		

def ONE2(A):
	if wrong_nu == 0:
		if M_nu == 0.0:
			return 0.
		Gnu4 = Gamma_nu**4
		a4 = A*A*A*A
		pi4 = (math.pi)**4
		y = (M_nu*A)/(Gamma_nu*kb*N_nu*Tcmb_0)
		return ((15.*Gnu4*N_nu*Or0/(a4*pi4))*np.interp(y,y1,FF1))
	
	else:
		return On0/(A*A*A)

Om0 = Ocb0 + On0
fb = Ob0/(Ocb0)
fc = Oc0/(Ocb0)
fcb = Ocb0/Om0
fn = On0/(Om0)
O_lambda_0 = 1 - Om0 - OR0
print(Om0,O_lambda_0,OR0, fn + fcb)

def E2(A,ONE2_current):
	if wrong_nu == 1:
		if M_nu == 0.0: ONE2_current = 0.0
		else:
			Gnu4 = Gamma_nu**4
			a4 = A**4
			pi4 = (math.pi)**4
			y = (M_nu*A)/(Gamma_nu*kb*N_nu*Tcmb_0)
			ONE2_current = ((15.*Gnu4*N_nu*Or0/(a4*pi4))*np.interp(y,y1,FF1))
	
	else: ONE2_current = ONE2(A)
	
	return OR0/(A*A*A*A) + Ocb0/(A*A*A) + O_lambda_0*(pow(A,-3.0*(1.0+w0+wa))*np.exp(3.0*wa*(A-1.0))) + ONE2_current

def Ocb(A, E2):
	return ( (Ocb0)/(A*A*A) )/E2(A,ONE2(A))

def Or(A,E2):
	return ( (OR0)/(A*A*A*A) )/E2(A,ONE2(A))

def On(A,E2):
	return ONE2(A)/E2(A,ONE2(A))

def O_lambda(A, E2):
	return ( (O_lambda_0*(pow(A,-3.0*(1.0+w0+wa))*np.exp(3.0*wa*(A-1.0)))) )/E2(A,ONE2(A))


def func_1_3w(A):
	if wrong_nu!= 0:
		return 1.0
	else:
		if M_nu == 0.0:
			return 1.
		else:
			y = M_nu*A/(Gamma_nu*kb*N_nu*Tcmb_0)
			return y*(np.interp(y,y1,GG1)/np.interp(y,y1,FF1))



def A_func(Ocb, Or, On, O_lambda, A, E2, rk_1_plus_3w):	
	if wrong_nu == 1:
	        if M_nu == 0.0:
        	    ON = 0.0
        	else:
        	    Gnu4 = Gamma_nu**4
        	    a4 = A**4
        	    pi4 = (math.pi)**4
        	    y = (M_nu * A) / (Gamma_nu * kb * N_nu * Tcmb_0)
        	    ON = (15.0 * Gnu4 * N_nu * Or0) / (a4 * pi4) * np.interp(y,y1,FF1)
        	    ON /= E2(A,ONE2(A))
        	
        	if M_nu == 0.0:
        	    rk_1_plus_3w = 1.0
        	else:
        	    rk_1_plus_3w = 2.0 - y * (np.interp(y,y1,GG1)/np.interp(y,y1,FF1))
	
	else: ON = On(A, E2)
	return 0.5*(Ocb(A, E2) + 2.*Or(A, E2) + (1.+3.*(w0+wa*(1.-A)))*O_lambda(A, E2) + rk_1_plus_3w*ON - 2.)

def B(Ocb, On, A, E2):
	return -3.0*( Ocb(A,E2) + func_1_3w(A)*On(A,E2) )/2.0

def coupled_odes(t, y):
	db_i, xb_i, dc_i, xc_i, dn_i, xn_i = y
	
	A = np.exp(t)
	if M_nu == 0.0:
		k_fs = 1.0
	else:
		k_fs = -(E2(A,ONE2(A))/(1.34423*1.34423))*(A**4)*B(Ocb, On, A, E2)*(M_nu*M_nu/9.)
	
	rk_1_3w = func_1_3w(A)
	rk_1_plus_3w = 2.0 - rk_1_3w
	
	ddb_idt = -1.0*xb_i
	dxb_idt = B(Ocb, On, A, E2)*(fcb*((Ob0/(Ob0+Oc0))*db_i + (Oc0/(Ob0+Oc0))*dc_i) +fn*dn_i) + A_func(Ocb, Or, On, O_lambda, A, E2, rk_1_plus_3w)*xb_i  
	
	ddc_idt = -1.0*xc_i
	dxc_idt = B(Ocb, On, A, E2)*(fcb*((Ob0/(Ob0+Oc0))*db_i + (Oc0/(Ob0+Oc0))*dc_i) +fn*dn_i) + A_func(Ocb, Or, On, O_lambda, A, E2, rk_1_plus_3w)*xc_i
	
	ddn_idt = -1.0*xn_i
	if M_nu == 0.:
		dxn_idt = B(Ocb, On, A, E2)*(fcb*((Ob0/(Ob0+Oc0))*db_i + (Oc0/(Ob0+Oc0))*dc_i) +fn*dn_i) + A_func(Ocb, Or, On, O_lambda, A, E2, rk_1_plus_3w)*xn_i
	else:
		dxn_idt = B(Ocb, On, A, E2)*(fcb*((Ob0/(Ob0+Oc0))*db_i + (Oc0/(Ob0+Oc0))*dc_i) + (fn - (k_i*k_i/k_fs) )*dn_i) + A_func(Ocb, Or, On, O_lambda, A, E2, rk_1_plus_3w)*xn_i
	
	return [ddb_idt, dxb_idt, ddc_idt, dxc_idt, ddn_idt, dxn_idt]


T1 = timeit.default_timer() - start_time
print('T1 = ',T1)


# Define the working directory
working_dir = str(Class)

# Define the executable and parameter file
executable = "./class"
parameter_file = "parameters_2F_class.ini"

# Run the command in the CLASS directory
try:
    result = subprocess.run(
        [executable, parameter_file],
        cwd=working_dir,  # Change to the CLASS directory before executing
        capture_output=True,
        text=True
    )

    # Print standard output
    print("STDOUT:\n", result.stdout)

    # Print errors if any
    if result.stderr:
        print("STDERR:\n", result.stderr)

    # Check exit status
    if result.returncode != 0:
        print(f"Process finished with error code {result.returncode}.")
    else:
        print("Execution completed successfully.")

except Exception as e:
    print(f"Exception occurred: {e}")

print('DONE')



n = 20
z_initial = 99.0
z_final = 0.0
z_output = np.linspace(z_initial, z_final, num=n)



file = open(Class/"output_2F/_z1_tk.dat", "r")
k = np.genfromtxt(Class/"output_2F/_z1_tk.dat", usecols=0)
db = np.genfromtxt(Class/"output_2F/_z1_tk.dat", usecols=2)
dc = np.genfromtxt(Class/"output_2F/_z1_tk.dat", usecols=3)
dn = np.genfromtxt(Class/"output_2F/_z1_tk.dat", usecols=6)
knum = len(k)

file = open(Class/"output_2F/_z2_tk.dat", "r")
dbplus = np.genfromtxt(Class/"output_2F/_z2_tk.dat", usecols=2)
dcplus = np.genfromtxt(Class/"output_2F/_z2_tk.dat", usecols=3)
dnplus = np.genfromtxt(Class/"output_2F/_z2_tk.dat", usecols=6)

file = open(Class/"output_2F/_z3_tk.dat", "r")
dbminus = np.genfromtxt(Class/"output_2F/_z3_tk.dat", usecols=2)
dcminus = np.genfromtxt(Class/"output_2F/_z3_tk.dat", usecols=3)
dnminus = np.genfromtxt(Class/"output_2F/_z3_tk.dat", usecols=6)


beta_b = np.array(db)/np.array(dc)
beta_n = np.array(dn)/np.array(dc)

db_99 = 1.*beta_b
dc_99 = [1.]*knum
dn_99 = 1.*beta_n

def lin_interp_between(x, x0, x1, y0, y1):
	return y1 + ( (y1 - y0)/(x1 - x0) )* (x - x1)



t = 0
Db = [0]*knum
Dc = [0]*knum
Dn = [0]*knum
Xb = [0]*knum
Xc = [0]*knum
Xn = [0]*knum
factor_0_m = [0]*knum
factor_0_b = [0]*knum
factor_0_c = [0]*knum
xb_99 = [0]*knum
xc_99 = [0]*knum
xn_99 = [0]*knum




z_output = 99.0
zmin = z_output - 2.
zmax = z_output + 2.

xmin = np.log(1. + zmin)
xmax = np.log(1. + zmax)
bc_step = (xmax - xmin) / (50 - 1)

bc_zz = np.array( np.exp(xmin + np.arange(50) * bc_step) - 1 )

zminus = np.delete(bc_zz, 49)
zplus = np.delete(bc_zz, 0)

fB = np.zeros((knum, 49))
fC = np.zeros((knum, 49))
fN = np.zeros((knum, 49))



for j in range(1,50):
	file = open(Class/f"output_2F/_z{j+3}_tk.dat", "r")
	k = np.genfromtxt(Class/f"output_2F/_z{j+3}_tk.dat", usecols=0)
	dbminus = np.genfromtxt(Class/f"output_2F/_z{j+3}_tk.dat", usecols=2)
	dcminus = np.genfromtxt(Class/f"output_2F/_z{j+3}_tk.dat", usecols=3)
	dnminus = np.genfromtxt(Class/f"output_2F/_z{j+3}_tk.dat", usecols=6)
	
	file = open(Class/f"output_2F/_z{j+4}_tk.dat", "r")
	dbplus = np.genfromtxt(Class/f"output_2F/_z{j+4}_tk.dat", usecols=2)
	dcplus = np.genfromtxt(Class/f"output_2F/_z{j+4}_tk.dat", usecols=3)
	dnplus = np.genfromtxt(Class/f"output_2F/_z{j+4}_tk.dat", usecols=6)
	
	knum = len(k)
	
	for l in range(0,knum):
		fB[l,j-1] = (np.log(dbplus[l]/dbminus[l])/np.log((1.+zminus[j-1])/(1.+zplus[j-1])))
		fC[l,j-1] = (np.log(dcplus[l]/dcminus[l])/np.log((1.+zminus[j-1])/(1.+zplus[j-1])))
		if (N_nu!=0):
			fN[l,j-1] = (np.log(dnplus[l]/dnminus[l])/np.log((1.+zminus[j-1])/(1.+zplus[j-1])))
		else:
			fN[l,j-1] = 0.
		
	

y_b = [0]*knum
y_c = [0]*knum
y_n = [0]*knum
for l in range(0,knum):
	x = (zplus+zminus)/2.	
	y = fB[l,:]
	#print(x,y)
	# Perform linear least squares fitting (degree 1 means linear fit)
	a, b = np.polyfit(x, y, 1)
	
	# Generate the fitted line
	y_fit = a * x + b
	y_b[l] = np.array(a*99.0 + b, dtype=np.float64)



for l in range(0,knum):
	x = (zplus+zminus)/2.	
	y = fC[l,:]
	
	# Perform linear least squares fitting (degree 1 means linear fit)
	a, b = np.polyfit(x, y, 1)
	
	# Generate the fitted line
	y_fit = a * x + b
	y_c[l] = np.array(a*99.0 + b, dtype=np.float64)



for l in range(0,knum):
	x = (zplus+zminus)/2.	
	y = fN[l,:]
	
	# Perform linear least squares fitting (degree 1 means linear fit)
	a, b = np.polyfit(x, y, 1)
	
	# Generate the fitted line
	y_fit = a * x + b
	y_n[l] = np.array(a*99.0 + b, dtype=np.float64)

FB = y_b
FC = y_c
FN = y_n

T2 = timeit.default_timer() - start_time
print('Time taken to read CLASS code = ',T2)




file_1 = open("./2F_ICs/Pm_99_30.txt","w")
file_2 = open("./2F_ICs/Tk_99_30.txt","w")
file_3 = open("./2F_ICs/fb_30.txt","w")
file_4 = open("./2F_ICs/fc_30.txt","w")
file_5 = open("./2F_ICs/fn_30.txt","w")
file_6 = open("./2F_ICs/Pm_30_z=0.txt","w")

norm = [0.]*(knum+1)
delta_c = [0.]*(knum+1)
delta_b = [0.]*(knum+1)
growth_b = [0.]*(knum+1)
growth_c = [0.]*(knum+1)
growth_b_99 = [0.]*(knum+1)
growth_c_99 = [0.]*(knum+1)
delta_n = [0.]*(knum+1)
growth_n = [0.]*(knum+1)
growth_n_99 = [0.]*(knum+1)
delta_m = [0.]*(knum+1)
Dm = [0.]*(knum+1)

file = open(Class/"output_2F/_z54_pk.dat", "r")
pm_0 = np.genfromtxt(Class/"output_2F/_z54_pk.dat", usecols=1)


for i in range(0,knum):
	file_6.write('	'+str(k[i])+'	'+str(pm_0[i])+'\n')


a_f = np.log(1./(1.+z_final))
a_i = np.log(1./(1.+z_initial))
h1 = ( a_f - a_i )/(n-1)
a = a_i
zminus = z_initial - 0.1
zplus = z_initial + 0.1

for i in range(0,knum):
	a = a_i
	db_i = 1.*beta_b[i]
	dc_i = 1.
	dn_i = 1.*beta_n[i]
	#break
	xb_i = -db_i*FB[i]
	xc_i = -dc_i*FC[i]
	xn_i = -dn_i*FN[i]
	
	db_1 = db_i
	dc_1 = dc_i
	dn_1 = dn_i
	xb_1 = xb_i
	xc_1 = xc_i
	xn_1 = xn_i
	
	a_1 = a - h1
	if np.exp(a) > 0.01 :
        	x_0_exp = 1.0 / np.exp(a_1) - 1.0
        	x_1_exp = 1.0 / np.exp(a) - 1.0
        	db_99[i] = lin_interp_between(0., x_0_exp, x_1_exp, db_1, db_i)
        	dc_99[i] = lin_interp_between(0., x_0_exp, x_1_exp, dc_1, dc_i)
        	dn_99[i] = lin_interp_between(0., x_0_exp, x_1_exp, dn_1, dn_i)
        	xb_99[i] = lin_interp_between(0., x_0_exp, x_1_exp, xb_1, xb_i)
        	xc_99[i] = lin_interp_between(0., x_0_exp, x_1_exp, xc_1, xc_i)
        	xn_99[i] = lin_interp_between(0., x_0_exp, x_1_exp, xn_1, xn_i)
	
	k_i = k[i]
	
	initial_conditions = [db_i, xb_i, dc_i, xc_i, dn_i, xn_i]
	t_eval = np.linspace(a_i, a_f, n)
	
	solver = ode(coupled_odes)
	solver.set_integrator('dopri5', atol=1e-6, rtol=1e-4)
	solver.set_initial_value(initial_conditions, a_i)
	
	# Integrate the ODEs
	solution = np.zeros((n, len(initial_conditions) + 1))
	solution[0, 0] = a_i
	solution[0, 1:] = initial_conditions
	idx = 1
	
	while solver.successful() and solver.t < a_f:
		solver.integrate(t_eval[idx])
		solution[idx, 0] = solver.t
		solution[idx, 1:] = solver.y
		idx += 1
	
	db_i, xb_i, dc_i, xc_i, dn_i, xn_i = solution[n-1, 1], solution[n-1, 2], solution[n-1, 3], solution[n-1, 4], solution[n-1, 5], solution[n-1, 6]
	
	
	Db[i] = db_i
	Dc[i] = dc_i
	Dn[i] = dn_i
	Xb[i] = xb_i
	Xc[i] = xc_i
	Xn[i] = xn_i
	
	################## writng in output files ##########################
	
	norm[i] = 1./( fcb*( Oc0/(Ob0+Oc0)*Dc[i] + Ob0/(Ob0+Oc0)*Db[i] ) + fn*Dn[i] )
	delta_c[i] = Dc[i]*norm[i]
	delta_b[i] = Db[i]*norm[i]
	delta_n[i] = Dn[i]*norm[i]
	delta_m[i] = fcb*( Oc0/(Ob0+Oc0)*delta_c[i] + Ob0/(Ob0+Oc0)*delta_b[i]) + fn*delta_n[i]
	
	delta_c_99 = dc_99[i]*norm[i]
	delta_b_99 = db_99[i]*norm[i]
	delta_n_99 = dn_99[i]*norm[i]
	delta_m_99 = fcb*(Ob0/(Ob0+Oc0)*delta_b_99 + Oc0/(Ob0+Oc0)*delta_c_99) + fn*delta_n_99
	
	growth_b_99[i] = -xb_99[i]/db_99[i]
	growth_c_99[i] = -xc_99[i]/dc_99[i]
	growth_n_99[i] = -xn_99[i]/dn_99[i]
	growth_m_99 = -(   (1-fn)*(fc*xc_99[i] + fb*xb_99[i]) + fn*xn_99[i]   )/( (1-fn)*(fc*dc_99[i] + fb*db_99[i]) + fn*xn_99[i]  )
	
	Pm_99 = (delta_m_99)*(delta_m_99)*pm_0[i]
	
	tk_b = delta_b_99/(fn*delta_n_99 + (1-fn)*(fb*delta_b_99 + fc*delta_c_99) )
	tk_c = delta_c_99/(fn*delta_n_99 + (1-fn)*(fb*delta_b_99 + fc*delta_c_99) )
	tk_n = delta_n_99/(fn*delta_n_99 + (1-fn)*(fb*delta_b_99 + fc*delta_c_99) )
	
	file_1.write('	'+str(k[i])+'	'+str(Pm_99)+'\n')
	file_2.write('	'+str(k[i])+'	'+str(tk_c)+'	'+str(tk_b)+'	'+str(0.0)+'	'+str(0.0)+'	'+str(tk_n)+'	'+str(1.0)+'\n')
	file_3.write('	'+str(k[i])+'	'+str(growth_b_99[i])+'\n')
	file_4.write('	'+str(k[i])+'	'+str(growth_c_99[i])+'\n')
	file_5.write('	'+str(k[i])+'	'+str(growth_n_99[i])+'\n')
	
	

########	Hubble in format of ICs for Gadget-3
H = open("./2F_ICs/Hz.txt","w")
z_nstep = 1000
z_step = (np.log(1.0 + z_initial) - np.log(1.0 + z_final)) / (z_nstep - 1)
H0 = h*100.
for i in range(z_nstep-1,-1,-1):
	zz = np.exp(np.log(1.0 + z_final) + i * z_step)
	r12 = ONE2(1.0 / zz)
	e2 = E2(1.0 / zz, r12)
	r2 = 0.1 * np.sqrt(e2)
	H.write('	'+str(zz)+'	'+str(-1.0)+'	'+str(r2)+'\n')


########	Hubble in format of RePS
'''
H = open("./2F_ICs/Hz.txt","w")
z_nstep = 1000
z_step = (np.log(1.0 + z_initial) - np.log(1.0 + z_final)) / (z_nstep - 1)
H0 = h*100.
print(H0)
for i in range(z_nstep):
	zz = np.exp(np.log(1.0 + z_final) + i * z_step)
	e2 = E2(1.0 / zz, ONE2(1.0 / zz))
	r2 = 1. * H0 * np.sqrt(e2)
	H.write('	'+str(zz-1.)+'	'+str(r2)+'\n')
'''
T3 = timeit.default_timer() - start_time
print('Total time taken = ',T3)
