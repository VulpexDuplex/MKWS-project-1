# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 15:53:24 2017

@author: Lukasz Osinski
"""
import csv
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def airm(temp, pressure, water_mol, ch4_mol, plot):
    print("--------------------------------------\n")
    size=10000 #number of iterations in the simulationg advancing loop
    
    #Mixture to ignite
    gas = ct.Solution('gri30.xml')
    
    n2=0.78/0.21
    ar=0.01/0.21
    d="CH4:"+str(ch4_mol)+", O2:1, N2:"+str(n2)+", AR:"+str(ar)+", H2O:"+str(water_mol)  

    print (d+"\n")
         
    gas.TPX = temp, pressure, d #setting gas temperature as 'temp, pressure as 'pressure' and molar fractions as in 'd'
    gas.name = "Methane-air-water mixture"
    if plot==2:
        print (gas())

    #adding reactor
    r = ct.Reactor(contents=gas, name='reactor', volume=0.02)
    
    print ('Water mass fraction = %10.3f %%\n' % (r.thermo['H2O'].Y*100))
    print ('Stoichiometry  index (phi) = %10.3f\n' % ( (r.thermo['CH4'].X/r.thermo['O2'].X)/(0.5) ) )
    
    #preparing simulation
    times = np.zeros(size)
    data = np.zeros((size,7))   
    time = 0.0 #starting time
    counter=size #setting 'counter' to 'size' (number of the for loop iterations)
    #different basic time steps for different initial temperatures
    if temp>1400:
        current_step = 1.e-5 #shorter timestep for high temperatures
    elif temp>1100:
        current_step = 1.e-4 #shorter timestep for high temperatures
    else:
        current_step = 1.e-2 #basic time step
    dP=0
    
    f=open('data.txt', 'w') #clearing file
    #simulation
    sim = ct.ReactorNet([r])

    #advance in time - the for loop
    print('%10s %10s %10s %14s' % ('t [s]','T [K]','P [Pa]','u [J/kg]'))    
    f.write("Temperatura poczatkowa wynosi %10.3f\n\n" % (temp))
    f.write ('%10s %10s %10s %10s %10s %10s %10s %10s' % ('t [s]','step [s]','dT [K]','T [K]', 'dP [Pa]', 'P [Pa]', 'licznik', 'iteration\n'))
    stepchange=0
    #different shorter time steps for different initial temperatures
    if temp>1400:
        short_step=1.e-6
    elif temp>1100:
        short_step=2.5*1.e-6 #10 times shorter for temepratures above 1200
    else:
        short_step=1.e-4 #shorter time step [s]       

    itafex=((10*1.e-2)/short_step) #going through this many iterations at specified 
    #timestep='short_step' should be equal to advancing the simulation by 0.1 seconds
    for n in range(size):
        time+= current_step
        sim.advance(time)
        times[n] = time # time in s
        data[n,0] = r.T #temperature [K]
        data[n,1] = r.thermo.P-pressure #relative pressure [Pa]
        data[n,2:5] = r.thermo['O2','H2','CH4'].X #mole fractions (X) of O2, H and CH4 [%]
        #data[n,5] = ignitest(time) #wydatek wodoru
        if n>0: #because data[n-1,1] for n=0 is data[-1,1] and that is impossible
            dP=data[n,1]-data[n-1,1] #current pressure minus pressure from previous iteration, in [Pa]
            data[n,6] = dP/current_step/1.e+6 #dP/dt  - maximum rate of pressure rise  in [Mpa/s]
            
        if r.T>temp*1.5 and n<counter+itafex and stepchange==4: #for the itafex number of interations 
        #between explosion and going back to 1.e-2 tstep, that is equal to 0.0225 [s]
            f.write('%10.7f %10.3e %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n' % (sim.time,
            current_step, data[n,0]-data[n-1,0], r.T, dP, r.thermo.P, counter, n)) 
        elif r.T>temp*1.5 and n>=counter+itafex: #setting time step back to 1.e-2, 0.00225 seconds after the explosion
            f.write('%10.7f %10.3e %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n' % (sim.time,
            current_step, data[n,0]-data[n-1,0], r.T, dP, r.thermo.P, counter, n)) 
            if stepchange==4:
                f.write('setting step back to 1.e-2\n')
                current_step=1.e-2
                stepchange=2
                break
            
        if r.T<=temp*1.5: #writing to file reactor state at each iteration before the explosion 
            f.write('%10.7f %10.3e %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n' % (sim.time,
            current_step, data[n,0]-data[n-1,0], r.T, dP, r.thermo.P, counter, n))      
        elif r.T>temp*1.5 and counter==size: #during the peak, one iteration where counter is still set as size
            f.write('%10.7f %10.3e %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n' % (sim.time,
            current_step, data[n,0]-data[n-1,0], r.T, dP, r.thermo.P, counter, n)) 
            f.write('Setting current_step to 5.e-4\n')
            current_step=short_step
            counter=n
            stepchange=4
         
        
        if n==0: #printing the state of reactor at the beginning of sumiulaton
            print('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T, r.thermo.P, r.thermo.u)) 
        if plot==2: #if argument 'plot' of function "airm" is set to 2, state of the reactor will be printed once a 1000 iterations
            if n % 1000 == 3:
                print('%10.3e %10.3f %10.3f %14.6e' % (sim.time, r.T, r.thermo.P, r.thermo.u))
        
                
    f.close() #closing the data file
                   
    #extracing maximal values out of data tables
    Tmax=max(data[:,0]) #maximal temperature
    Tend=data[n,0] #temperature at the end of simulation
    Pmax=max(data[:,1]) #maximal pressure (pressure computed - initial pressure)
    dPmax=max(data[:,6]) #maximum rate of pressure rise 
    index_max = np.argmax(data[:,6])
    adt=times[index_max]
    if Pmax<(pressure*1.05-pressure):
        adt=None
    print ('Tmax = %s [K]\nTend = %s [K]' % (Tmax, Tend))
    print ('Explosion Pmax = %s [Pa]\n' % (Pmax))
    print ('Explosion dP/dt max = %s [MPa/s]\n' % (dPmax))
    print ('Autoiginition delay time = %s [s]\n' % (adt))
    output=[Tmax, Pmax, dPmax, adt] #this is what this function ("airm") returns at the end
    
    #calculating CH4 volume fraction in mixture
    o2=1
    n2=0.78/0.21
    ar=0.01/0.21
    ch4m_data=(ch4_mol)/(ch4_mol+o2+n2+ar)*100
    print ('CH4 volume fraction = %s %%' % (ch4m_data))     
          
    #plotting graphs if argument 'plot' of "airm" function is set to 1 or 2
    if plot==1 or plot==2:
        n=n+1 #because "0:n" means "from 0 to n-1"
        plt.clf()
        plt.figure(figsize=(6,6))
        plt.subplot(2, 2, 1)
        plt.plot(times[:], data[:,0])
        plt.xlabel('Time [s]')
        plt.ylabel('Temperature (K)')
        plt.subplot(2, 2, 2)
        plt.plot(times[0:n], data[0:n,1]/1.e+6)
        plt.xlabel('Time [s]')
        plt.ylabel('Pressure (MPa)')
        plt.subplot(2, 2, 3)
        plt.plot(times[0:n], data[0:n,2])
        plt.xlabel('Time [s]')
        plt.ylabel('O2 mole fraction [%]')
        plt.subplot(2, 2, 4)
        plt.plot(times[0:n],data[0:n,4])
        plt.xlabel('Time [s]')
        plt.ylabel('CH4 mle fraction [%]')
        plt.tight_layout()
        plt.show()
        
        plt.clf()
        plt.figure(figsize=(6,6))
        plt.subplot(2, 2, 1)
        plt.plot(times[0:n], data[0:n,3])
        plt.xlabel('Time [s]')
        plt.ylabel('H mole fraction [s]')
        plt.tight_layout()
        plt.show()

    return output #Tmax, Pmax, (dP/dt)max

"""
f=open('data.txt', 'a') #appends data.txt
print ('Tmax, Pmax, dp/dt max, autignition delay time')
f.write(' T=1500, p=1atm, water=1mol, chr=0.5 mol phi=0.5\n')
print (airm(temp=1500, pressure=ct.one_atm, water_mol=1, ch4_mol=0.25, plot=2))
f.close()
"""
#
#THE FIRST LOOP - in f(CH4 mole fraction)
#
lp=41 #number of iterations for 1st loop
#preparing space for data
pdata = np.zeros((lp,4)) #without water
ch4m_data = np.zeros(lp) #without water
ch4m_w_data = np.zeros(lp) #with water
ch4m_2w_data = np.zeros(lp) #with water
pdata_w = np.zeros((lp,4)) #with 1 mol of water
pdata_2w = np.zeros((lp,4)) #with 2 moles of water
o2=1 
n2=0.78/0.21
ar=0.01/0.21
wt=1 #mole fraction of water in mixture
#in this loop CH4 mole fraction is changing from 0 to lp*i/2, with and without addition of water
for i in range(lp):
    pdata[i,:]=airm(temp=1000, pressure=ct.one_atm, water_mol=0, ch4_mol=i/20, plot=0) #without water
    pdata_w[i,:]=airm(temp=1000, pressure=ct.one_atm, water_mol=wt, ch4_mol=i/20, plot=0) #with water
    pdata_2w[i,:]=airm(temp=1000, pressure=ct.one_atm, water_mol=2*wt, ch4_mol=i/20, plot=0) #with water
    ch4m_data[i]=(i/20)/(i/20+o2+n2+ar)*100      #mole=volume fraction of CH4 in mixture, when there is no water
    ch4m_w_data[i]=(i/20)/(i/20+o2+n2+ar+wt)*100 #mole=volume fraction of CH4 in mixture, when there is 1 mole water
    ch4m_2w_data[i]=(i/20)/(i/20+o2+n2+ar+2*wt)*100 #mole=volume fraction of CH4 in mixture, when there are 2 moles of water

#
#THE SECOND LOOP - in f(T)
#

l2p=26 #number of iterations for 2nd loop
#preparing space for data
pdata2 = np.zeros((l2p,4)) 
pdata2_w = np.zeros((l2p,4)) #with 1 moles of water
pdata2_2w = np.zeros((l2p,4)) #with 2 moles of  water
temps = np.zeros(l2p)
#loop through start temperatures
f=open('data.txt', 'a') #clears data.txt
f.write(' S Y M U L A C J A\n')
for j in range(l2p):
    current_temp=850+j*50
    pdata2[j,:]=airm(temp=current_temp, pressure=1*ct.one_atm, water_mol=0, ch4_mol=0.5, plot=0)
    pdata2_w[j,:]=airm(temp=current_temp, pressure=1*ct.one_atm, water_mol=1*wt, ch4_mol=0.5, plot=0)
    pdata2_2w[j,:]=airm(temp=current_temp, pressure=1*ct.one_atm, water_mol=2*wt, ch4_mol=0.5, plot=0)
    temps[j]=current_temp
    print('Start temperature = %10.3f C\n' % (current_temp))
        
f.close()

"""
plt.clf()
plt.figure(figsize=(20,20))
plt.subplot(2, 2, 1)
plt.plot(temps, pdata2[:,1])
plt.xlabel('Temperature [K]')
plt.ylabel('P_max [Pa]')
plt.subplot(2, 2, 2)
plt.plot(temps, pdata2[:,0])
plt.xlabel('Temperature [K]')
plt.ylabel('T_max [K]')
plt.subplot(2, 2, 3)
plt.plot(temps, pdata2[:,2])
plt.xlabel('Temperature [K]')
plt.ylabel('dP/dt max [MPa/s]')
plt.tight_layout()
plt.show()
"""

#
#THE THIRD LOOP - in f(pressure)
#
l3p=26 #number of loop iterations
#preparing space for data
pdata3 = np.zeros((l3p,4))
pdata3_w = np.zeros((l3p,4)) #with 1 mole of water
pdata3_2w = np.zeros((l3p,4)) #with 2 moles of water
press = np.zeros(l3p)
#wt=1 #mole fraction of water in mixture
#loop through start pressures
for k in range(l3p):
    current_press=(0.5*ct.one_atm)+ct.one_atm*k #current start pressure in [Pa]
    pdata3[k,:]=airm(temp=850, pressure=current_press, water_mol=0, ch4_mol=0.5, plot=0)
    pdata3_w[k,:]=airm(temp=850, pressure=current_press, water_mol=wt, ch4_mol=0.5, plot=0)
    pdata3_2w[k,:]=airm(temp=850, pressure=current_press, water_mol=2*wt, ch4_mol=0.5, plot=0)
    press[k]=current_press #in [Pa]
       
"""
plt.clf()
plt.figure(figsize=(20,20))
plt.subplot(2, 2, 1)
plt.plot(press, pdata3[:,1])
plt.xlabel('Pressure [K]')
plt.ylabel('P_max [Pa]')
plt.subplot(2, 2, 2)
plt.plot(press, pdata3[:,0])
plt.xlabel('Pressure [K]')
plt.ylabel('T_max [K]')
plt.tight_layout()
plt.show()
"""

#plotting graphs from data from loops 1-3 in one place 
print ("The first loop\n\n")
plt.clf()
plt.figure(figsize=(20,20)) #setting size of figures
#w/o water
plt.subplot(2, 2, 1)
plt.plot(ch4m_data, pdata[:,1]/1.e+6, label='without water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('P_max [Mpa]')
plt.subplot(2, 2, 2)
plt.plot(ch4m_data, pdata[:,0], label='without water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('T_max [K]')
#with 1 mole of water
plt.subplot(2, 2, 1)
plt.plot(ch4m_w_data, pdata_w[:,1]/1.e+6, label='15.97% of water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('P_max [Mpa]')
plt.subplot(2, 2, 2)
plt.plot(ch4m_data, pdata_w[:,0], label='15.97% of water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('T_max [K]')
#with 2 moles of water
plt.subplot(2, 2, 1)
plt.plot(ch4m_2w_data, pdata_2w[:,1]/1.e+6, label='27.54% of water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('P_max [Mpa]')
plt.legend()
plt.subplot(2, 2, 2)
plt.plot(ch4m_data, pdata_2w[:,0], label='27.54% of water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('T_max [K]')
plt.legend()
plt.tight_layout()
plt.show()

plt.clf()
plt.figure(figsize=(10,10))
plt.subplot(1, 1, 1)
plt.plot(ch4m_data, pdata[:,3]*1.e+3, label='without water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('Autoign delay [ms]')
plt.subplot(1, 1, 1)
plt.plot(ch4m_data, pdata_w[:,3]*1.e+3, label='15.97% of water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('Autoign delay [ms]')
plt.subplot(1, 1, 1)
plt.plot(ch4m_data, pdata_2w[:,3]*1.e+3, label='27.54% of water')
plt.xlabel('CH4 (%vol)')
plt.ylabel('Autoign delay [ms]')
plt.legend()
plt.tight_layout()
plt.show()


csv_file = 'graphdata.csv'
with open(csv_file, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['The first loop'])
    writer.writerow(['Without water'])
    writer.writerow(['CH4 [%]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(lp):
         writer.writerow([ch4m_data[i], pdata[i,0], pdata[i,1]/1.e+6, pdata[i,2], pdata[i,3]])
    writer.writerow(['With 1 mole of water'])
    writer.writerow(['CH4[%]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(lp):
         writer.writerow([ch4m_w_data[i], pdata_w[i,0], pdata_w[i,1]/1.e+6, pdata_w[i,2], pdata_w[i,3]])
    writer.writerow(['With 2 moles of water'])
    writer.writerow(['CH4[%]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(lp):
         writer.writerow([ch4m_2w_data[i], pdata_2w[i,0], pdata_2w[i,1]/1.e+6, pdata_2w[i,2],  pdata_2w[i,3]])


print ("The second loop\n\n")

plt.clf()
plt.figure(figsize=(20,20))
#w/o water
l1=plt.subplot(2, 2, 1)
plt.plot(temps, pdata2[:,1]/1.e+6, label='without water')
plt.xlabel('Temperature [K]')
plt.ylabel('P_max [MPa]')
plt.subplot(2, 2, 2)
plt.plot(temps, pdata2[:,0], label='without water')
plt.xlabel('Temperature [K]')
plt.ylabel('T_max [K]')
#with water
plt.subplot(2, 2, 2)
plt.plot(temps, pdata2_w[:,0], label='15.97% of water')
plt.xlabel('Temperature [K]')
plt.ylabel('T_max [K]')
plt.subplot(2, 2, 1)
plt.plot(temps, pdata2_w[:,1]/1.e+6, label='15.97% of water')
plt.xlabel('Temperature [K]')
plt.ylabel('P_max [MPa]')
#with 2 moles of water
plt.subplot(2, 2, 2)
plt.plot(temps, pdata2_2w[:,0], label='27.54% of water')
plt.xlabel('Temperature [K]')
plt.ylabel('T_max [K]')
plt.legend()
plt.subplot(2, 2, 1)
plt.plot(temps, pdata2_2w[:,1]/1.e+6, label='27.54% of water')
plt.xlabel('Temperature [K]')
plt.ylabel('P_max [MPa]')
plt.legend()
plt.tight_layout()
plt.show()

plt.clf()
plt.figure(figsize=(10,10))
plt.subplot(1, 1, 1)
plt.plot(temps, pdata2[:,3]*1.e+3, label='without water')
plt.xlabel('1000/T [1/K]')
plt.yscale('log')
plt.ylabel('Autoign delay [ms]')
plt.subplot(1, 1, 1)
plt.plot(temps, pdata2_w[:,3]*1.e+3, label='15.97% of water')
plt.xlabel('1000/T [1/K]')
plt.yscale('log')
plt.ylabel('Autoign delay [ms]')
plt.subplot(1, 1, 1)
plt.plot(temps, pdata2_2w[:,3]*1.e+3, label='27.54% of water')
plt.xlabel('1000/T [1/K]')
plt.ylabel('Autoign delay [ms]')
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()


csv_file = 'graphdata.csv'
with open(csv_file, 'a') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['The second loopr'])
    writer.writerow(['Without water'])
    writer.writerow(['T0 [K]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(l2p):
         writer.writerow([temps[i], pdata2[i,0], pdata2[i,1]/1.e+6, pdata2[i,2], pdata2[i,3]])
    writer.writerow(['With 1 mole of water'])
    writer.writerow(['T0 [K]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(l2p):
         writer.writerow([temps[i], pdata2_w[i,0], pdata2_w[i,1]/1.e+6, pdata2_w[i,2], pdata2_w[i,3]])
    writer.writerow(['With 2 moles of water'])
    writer.writerow(['T0 [K]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(l2p):
         writer.writerow([temps[i], pdata2_2w[i,0], pdata2_2w[i,1]/1.e+6, pdata2_2w[i,2], pdata2_2w[i,3]])


print ("The third loop\n\n")
plt.clf() #clear the current figure
#w/o water
plt.figure(figsize=(20,20))
plt.subplot(2, 2, 1)
plt.plot(press/1.e+6, pdata3[:,1]/1.e+6, label='without of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('P_max [MPa]')
plt.subplot(2, 2, 2)
plt.plot(press/1.e+6, pdata3[:,0], label='without of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('T_max [K]')
#with 1 mole of water
plt.subplot(2, 2, 1)
plt.plot(press/1.e+6, pdata3_w[:,1]/1.e+6, label='15.97% of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('P_max [MPa]')
plt.subplot(2, 2, 2)
plt.plot(press/1.e+6, pdata3_w[:,0], label='15.97% of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('T_max [K]')
#with 2 moles of water
plt.subplot(2, 2, 1)
plt.plot(press/1.e+6, pdata3_2w[:,1]/1.e+6, label='27.54% of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('P_max [MPa]')
plt.legend()
plt.subplot(2, 2, 2)
plt.plot(press/1.e+6, pdata3_2w[:,0], label='27.54% of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('T_max [K]')
plt.legend()
plt.tight_layout()
plt.show()       

plt.clf()
plt.figure(figsize=(10,10))
plt.subplot(1, 1, 1)
plt.plot(press/1.e+6, pdata3[:,3]*1.e+3, label='without water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('Autoign delay [ms]')
plt.subplot(1, 1, 1)
plt.plot(press/1.e+6, pdata3_w[:,3]*1.e+3, label='15.97% of water')
plt.xlabel('Pressure [MPa]')
plt.ylabel('Autoign delay [ms]')
plt.subplot(1, 1, 1)
plt.plot(press/1.e+6, pdata3_2w[:,3]*1.e+3, label='27.54% of water')
plt.xlabel('Pressure [MPa]')
plt.yscale('log')
plt.ylabel('Autoign delay [ms]')
plt.legend()
plt.tight_layout()
plt.show()       

csv_file = 'graphdata.csv'
with open(csv_file, 'a') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['The third loop'])
    writer.writerow(['Without water'])
    writer.writerow(['p0 [MPa]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(l3p):
         writer.writerow([press[i]/1.e+6, pdata3[i,0], pdata3[i,1]/1.e+6, pdata3[i,2], pdata3[i,3]])
    writer.writerow(['With 1 mole of water'])
    writer.writerow(['p0 [MPa]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(l3p):
         writer.writerow([press[i]/1.e+6, pdata3_w[i,0], pdata3_w[i,1]/1.e+6, pdata3_w[i,2], pdata3_w[i,3]])
    writer.writerow(['With 2 moles of water'])
    writer.writerow(['p0 [MPa]', 'Tmax [K]', 'Pmax [MPa]', 'dP/dt max [MPa/s]', 'adt [s]'])
    for i in range(l3p):
         writer.writerow([press[i]/1.e+6, pdata3_2w[i,0], pdata3_2w[i,1]/1.e+6, pdata3_2w[i,2], pdata3_2w[i,3]])
