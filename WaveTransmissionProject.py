import cmath
import numpy as np
import matplotlib.pyplot as plt

def plot_V_and_I(z,*V):
    
    
    figure, axis = plt.subplots(1,2, figsize = (20,10))
    
    axis[0].set_xlabel('z[m]')
    axis[0].set_ylabel('V[v]')
    axis[0].set_title('V(z)')
    axis[0].grid()
    axis[0].plot(z,V[0][0:N])
    if(len(V)>1):
        axis[0].plot(z,V[1][0:N])
        axis[0].legend(['Nummeric Calculation', 'Theoretical Value'])
    
    axis[1].set_xlabel('z[m]')
    axis[1].set_ylabel('I[A]')
    axis[1].set_title('I(z)')
    axis[1].grid()
    axis[1].plot(z,V[0][N:2*N])
    if(len(V)>1):
        axis[1].plot(z,V[1][N:2*N])
        axis[1].legend(['Nummeric Calculation', 'Theoretical Value']) 

    plt.show()

#define Q1 constants
l = 0.1
N = 1000
delta = l/N
vp = 1.5*10**8
f = 10*(10**9)
lamb = vp/f
Beta = 2*(np.pi)/lamb
C = 8.89*(10**-11)
L = 5*(10**-7)
V = 1
w = 2*(np.pi)*f
RG = 50
Z0 = 75 
Zin = 50
ZL = 150
Gamma_G = ((RG-Z0)/(RG+Z0))
Gamma_L = ((ZL-Z0)/(ZL+Z0))
epsilon = 0.0001

# define I-V linkage matrix
M = np.zeros((2*N,2*N), dtype=complex)
solution_vector = np.zeros(2*N, dtype = complex)
solution_vector[0] = 1

#M content according and border conditions

#Condition no.1
M[0,0] = 1
M[0,N] = RG

#Condition no.2
M[1,0] = -1/delta
M[1,1] = 1/delta
M[1,N] = (w*L*delta/2)*1j
M[1,N+1] = (w*L*delta/2)*1j

#Main content part 1
for i in range(2,N):
    M[i,i-2] = -1/(2*delta)
    M[i,i] = 1/(2*delta)
    M[i,N+i-1] = (w*L)*1j

#Condition no.3
M[N,N-1] = 1
M[N,2*N-1] = -ZL

#Condition no.4
M[N+1,N-1] = 1/delta
M[N+1,N-2] = -1/delta
M[N+1, 2*N-1] = (w*L*delta/2)*1j
M[N+1, 2*N-2] = (w*L*delta/2)*1j

#Main content part 2
for i in range(N+2,2*N):
    M[i,i-2] = -1/(2*delta)
    M[i,i] = 1/(2*delta)
    M[i,i-N-1] = (w*C)*1j

#solve
variables_vector = np.linalg.solve(M,solution_vector)

z = np.linspace(0,l,N)
plot_V_and_I(z,variables_vector)

#Section 1.e

#According to recitation 3
z = np.linspace(-l,0,N)
Vo_plus = 0.6/(np.exp(1j*Beta*l)*(1-((RG-Z0)/(RG+Z0))*(Gamma_L)*np.exp(-2j*Beta*l)))

#According to recitation 3
Vo = np.real(((Vo_plus)*np.exp(-1j*Beta*z))*(1+(1/3)*(np.exp(2j*Beta*z))))
Io = np.real(((Vo_plus)/Z0)*np.exp(-1j*Beta*z)*(1-(1/3)*(np.exp(2j*Beta*z))))

z = np.linspace(0,l,N)

theo_values = np.hstack((Vo,Io))

#Lets compare the results graphically:

plot_V_and_I(z,variables_vector,theo_values)

#Calculate relative errors:
error_v = np.zeros(N)
error_i = np.zeros(N)

for i in range(N):
    error_v[i] = abs((theo_values[i] - np.real(variables_vector[i]))/theo_values[i])
    error_i[i] = abs((theo_values[N+i] - np.real(variables_vector[N+i]))/theo_values[N+i])

#Remove outliers that are more than 3 times larger than the mean value. Only 16/2000 elements were removed! Only 0.8% of the data.

error_v = error_v[error_v<(3*error_v.mean())]
error_i = error_i[error_i<(3*error_i.mean())]

print(f'Average relative error between voltage values was {np.round(error_v.mean()*100,2)}%')
print(f'Average relative error between voltage values was {np.round(error_i.mean()*100,2)}%')

# Question 2

# Now define L_match, C_match that were calculated analytically in qustion 2.1
L_match = 5.3033*(10**-8)
C_match = 4.714*(10**-12)
z1 = np.sqrt(L_match/C_match)

#And get to business with calculatin and plotting question 2.2 results

f_range = np.linspace(5*10**9,15*10**9,1000)
w_range = 2*np.pi*f_range

#According to a formula presented in the PDF
Zin = z1*(ZL+1j*z1*np.tan(w_range*np.sqrt(L_match*C_match)*(l/2)))/(z1 + 1j*ZL*np.tan(w_range*np.sqrt(L_match*C_match)*(l/2)))

gamma_in_matched = abs((Zin - Z0)/(Zin + Z0))

plt.plot(f_range/10**9, gamma_in_matched)

set_high = 0
set_low = 0

for k in range(N//2):
    if(gamma_in_matched[N//2-1-k]>0.1 and not set_low):
        lower_bound = f_range[N//2-1-k]
        set_low = 1
    if(gamma_in_matched[N//2+k]>0.1 and not set_high):
        higher_bound = f_range[N//2+k]
        set_high = 1

sym_range = min(2*abs(higher_bound - 10**10),2*abs(lower_bound - 10**10))

print(f'The symmetric bandwidth around 10GHz is {np.round(sym_range/10**9, 2)} GHz')

plt.xlabel('frequency(GHz)')
plt.ylabel('Gamma_in')
plt.title=('Gamma_in(f)')
plt.grid()
plt.show()




















