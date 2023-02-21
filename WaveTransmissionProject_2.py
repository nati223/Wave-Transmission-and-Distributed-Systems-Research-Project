import cmath
import numpy as np
import matplotlib.pyplot as plt

def calculate_gamma(L1,C1,alpha1,alpha2,w):
    
    M = np.zeros((2*N, 2*N), dtype=complex)
    solutions_v = np.zeros(2*N, dtype=complex)
    solutions_v[0] = 1
    
    #Condition no.1
    M[0,0] = 1
    M[0,N] = RG

    #Condition no.2
    M[1,0] = -1/delta
    M[1,1] = 1/delta
    M[1,N] = (w*L*delta/2)*1j
    M[1,N+1] = (w*L*delta/2)*1j

    #Main content part 1
    for i in range(2,(N//2)):
        M[i,i-2] = -1/(2*delta)
        M[i,i] = 1/(2*delta)
        M[i,N+i-1] = (w*L)*1j
    
    #Condition no.3
    M[N//2,N//2-2] = -1/delta
    M[N//2,N//2-1] = 1/delta
    M[N//2,3*N//2-2] = (w*L*delta/2)*1j
    M[N//2,3*N//2-1] = (w*L*delta/2)*1j
    
    #Condition no.4
    M[N//2+1,N//2] = -1/delta
    M[N//2+1,N//2+1] = 1/delta
    M[N//2+1,3*N//2] = (w*L1*delta/2)*1j
    M[N//2+1,3*N//2+1] = ((w*(L1 + alpha1*(delta)**2))*delta/2)*1j
    
    k = 1
    for i in range((N//2)+2,N):
        M[i,i-2] = -1/(2*delta)
        M[i,i] = 1/(2*delta)
        #M[i,N+i-1] = (w*(L1 + alpha1*((i*delta - l/2)**2)))*1j
        M[i,N+i-1] = (w*(L1 + alpha1*((k*delta)**2)))*1j
        k+=1
        
    #Condition no.5
    M[N,N//2-1] = -1
    M[N,N//2] = 1
    
    #Condition no.6
    M[N+1,3*N//2-1] = -1
    M[N+1,3*N//2] = 1
     
    #Condition no.7
    M[N+2,N-1] = 1/delta
    M[N+2,N-2] = -1/delta
    M[N+2, 2*N-1] = (w*(L1 + alpha1*((l/2)**2))*delta/2)*1j
    M[N+2, 2*N-2] = (w*(L1 + alpha1*((l/2 - delta)**2))*delta/2)*1j

    #Main content part 2
    for i in range(N+3,3*(N//2)+1): 
        M[i,i-2] = -1/(2*delta)
        M[i,i] = 1/(2*delta)
        M[i,i-N-1] = (w*C)*1j
 
    #Condition no.8
    M[3*N//2+1, N-1] = 1
    M[3*N//2+1, 2*N-1] = -ZL
    
    k = 1
    for i in range(3*(N//2) + 2, 2*N):
        M[i,i-2] = -1/(2*delta)
        M[i,i] = 1/(2*delta)
        #M[i,i-N-1] = (w*(C1 + alpha2*((l/2 - (i-N-1)*delta)**2)))*1j
        M[i,i-N-1] = (w*(C1 + alpha2*((k*delta)**2)))*1j
        k+=1
    
    variables_v = np.linalg.solve(M,solutions_v)
    
    Zin = variables_v[N//2]/variables_v[3*N//2]
    
    Gamma_in = (Zin - Z0)/(Zin + Z0)
    
    return abs(Gamma_in)
    

#define constants
l = 0.1
N = 1000
F = N
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
ZL = 150
Gamma_G = ((RG-Z0)/(RG+Z0))
Gamma_L = ((ZL-Z0)/(ZL+Z0))
epsilon = 0.0001
#C_match = (0.2*10**10*75*np.sqrt(2))**-1
#L_match = ((75*np.sqrt(2))**2)*C_match
L_match = 5.3033*(10**-8)
C_match = 4.714*(10**-12)
f_range = np.linspace(5*10**9,15*10**9,F)
w_range = 2*np.pi*f_range

#First perform sanity check, to see that N=400 gives a good approximation of the original problem
Gamma_test = np.zeros(F)

for i in range(F):
    Gamma_test[i] = calculate_gamma(L_match,C_match,0,0,w_range[i])

set_low = 0
set_high = 0

for k in range(F//2):
    if(Gamma_test[F//2-1-k]>0.1 and not set_low):
        lower_bound = f_range[F//2-1-k]
        set_low = 1
    if(Gamma_test[F//2+k]>0.1 and not set_high):
        higher_bound = f_range[F//2+k]
        set_high = 1

sym_range = min(2*abs(higher_bound - 10**10),2*abs(lower_bound - 10**10))

print(f"Test bandwidth is: {np.round(sym_range/10**9, 2)} GHz")

plt.grid()
plt.xlabel('GHz')
plt.ylabel("Gamma_in")
plt.plot(f_range,Gamma_test)
plt.show()

#Let's get to our best result

win_v = np.zeros(F)

for k in range(F):
    win_v[k] = calculate_gamma(0.976*L_match,1.55*C_match,754*L_match,135*C_match,w_range[k])
    #win_v[k] = calculate_gamma(0.91*L_match,1.25*C_match,571*L_match,165*C_match,w_range[k]) # yields 4.67
    #win_v[k] = calculate_gamma(0.91*L_match,1.23*C_match,500*L_match,145*C_match,w_range[k]) - yields 4.47
    #win_v[k] = calculate_gamma(0.845*L_match,0.89*C_match,184*L_match,200*C_match,w_range[k])

set_low = 0
set_high = 0

for k in range(F//2):
    if(win_v[F//2-1-k]>0.1 and not set_low):
        lower_bound = f_range[F//2-1-k]
        set_low = 1
    if(win_v[F//2+k]>0.1 and not set_high):
        higher_bound = f_range[F//2+k]
        set_high = 1

sym_range = min(2*abs(higher_bound - 10**10),2*abs(lower_bound - 10**10))

print(f'\nBest results:\nhigher bound: {higher_bound/10**9} GHz\nlower_bound: {lower_bound/10**9} GHz')
print(f'Symmetric bandwidth is: {(np.round(sym_range/10**9,2))} GHz')

plt.xlim((5*10**9, 15*10**9))
plt.ylim((0,0.2))
plt.grid()
plt.plot(f_range,win_v)
plt.show()

#Let's get to our best result with diferent N and F - notice that it's identical to the code above

l = 0.1
N = 400
F = N//2
delta = l/N
Beta = 2*(np.pi)/lamb
f_range = np.linspace(5*10**9,15*10**9,F)
w_range = 2*np.pi*f_range

win_v_2 = np.zeros(F)

for k in range(F):
    win_v_2[k] = calculate_gamma(0.976*L_match,1.55*C_match,754*L_match,135*C_match,w_range[k])
    #win_v[k] = calculate_gamma(0.91*L_match,1.25*C_match,571*L_match,165*C_match,w_range[k]) # yields 4.67
    #win_v[k] = calculate_gamma(0.91*L_match,1.23*C_match,500*L_match,145*C_match,w_range[k]) - yields 4.47
    #win_v[k] = calculate_gamma(0.845*L_match,0.89*C_match,184*L_match,200*C_match,w_range[k])

set_low = 0
set_high = 0

for k in range(F//2):
    if(win_v_2[F//2-1-k]>0.1 and not set_low):
        lower_bound = f_range[F//2-1-k]
        set_low = 1
    if(win_v_2[F//2+k]>0.1 and not set_high):
        higher_bound = f_range[F//2+k]
        set_high = 1

sym_range = min(2*abs(higher_bound - 10**10),2*abs(lower_bound - 10**10))

print(f'\nBest results:\nhigher bound: {higher_bound/10**9} GHz\nlower_bound: {lower_bound/10**9} GHz')
print(f'Symmetric bandwidth is: {(np.round(sym_range/10**9,2))} GHz')

plt.xlim((5*10**9, 15*10**9))
plt.ylim((0,0.2))
plt.grid()
plt.plot(f_range,win_v_2)
plt.show()