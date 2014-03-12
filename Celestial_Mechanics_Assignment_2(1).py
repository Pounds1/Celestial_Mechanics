# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 10:56:51 2014

@author: user
"""
import matplotlib.pyplot as plt
G = 1.4880826*10**(-34)
N = 9


timestep = 7 #one week in days
counter_1  = timestep-1
counter_2 = 0
#defines attributes of each Particle
class particle_type:
    
    def __init__(self):
        
        self.mass = 0
        self.position = [[0,0,0],[0,0,0],[0,0,0]]
        self.velocity = [[0,0,0],[0,0,0],[0,0,0]]
        self.acceleration = [0,0,0]
        self.position_array_x = []
        self.position_array_y = []
        self.position_array_z = []

#assigns each particle to its attributes
particle_0 = particle_type()
particle_1 = particle_type()
particle_2 = particle_type()
particle_3 = particle_type()
particle_4 = particle_type()
particle_5 = particle_type()
particle_6 = particle_type()
particle_7 = particle_type()
particle_8 = particle_type()
#defines the array of each object Particle
Particles  = [particle_0,particle_1,particle_2,particle_3,particle_4,particle_5,particle_6,particle_7,particle_8]

#initial conditions
Particles[0].position[0][0] = -0.3919024742935512
Particles[0].position[0][1] = -0.05181639857131116
Particles[0].position[0][2] =  0.03152855008825588
Particles[0].velocity[0][0] = -0.001885143088580373
Particles[0].velocity[0][1] = -0.02663480354510647
Particles[0].velocity[0][2] = -0.002003022843923769
Particles[0].mass = 3.302*10**23
Particles[1].position[0][0] =  0.6835846910259694
Particles[1].position[0][1] =  0.2499480571599963
Particles[1].position[0][2] = -0.03605212331793326
Particles[1].velocity[0][0] = -0.006929784113964885
Particles[1].velocity[0][1] =  0.01894334581364429
Particles[1].velocity[0][2] =  0.0006584075451767032
Particles[1].mass = 48.685*10**23
Particles[2].position[0][0] =  0.09861071802361319
Particles[2].position[0][1] =  0.9836946517100500
Particles[2].position[0][2] = -0.00005732770495032725
Particles[2].velocity[0][0] = -0.01740133448679011
Particles[2].velocity[0][1] =  0.001632119558092957
Particles[2].velocity[0][2] =  0.0000001730821486604745
Particles[2].mass = 5.97219*10**24
Particles[3].position[0][0] = -0.1878080938397302
Particles[3].position[0][1] =  1.578841079634226
Particles[3].position[0][2] =  0.03758310419458916
Particles[3].velocity[0][0] = -0.01336168038088148
Particles[3].velocity[0][1] = -0.0004834846141641441
Particles[3].velocity[0][2] =  0.0003183673509975010
Particles[3].mass = 6.4185*10**23
Particles[4].position[0][0] = -5.439497306246083
Particles[4].position[0][1] = -0.1838809862035173
Particles[4].position[0][2] =  0.1225733966675222 
Particles[4].velocity[0][0] =  0.0001619864299930425
Particles[4].velocity[0][1] = -0.007193708963433175
Particles[4].velocity[0][2] =  0.00002613839059360070
Particles[4].mass = 1898.13*10**24
Particles[5].position[0][0] =  7.481228210520391
Particles[5].position[0][1] = -6.430452204421326
Particles[5].position[0][2] = -0.1853447452247906 
Particles[5].velocity[0][0] =  0.003334201313672095
Particles[5].velocity[0][1] =  0.004217870915174572
Particles[5].velocity[0][2] = -0.0002057696886504684
Particles[5].mass = 5.68319*10**26
Particles[6].position[0][0] =  6.031279452943370
Particles[6].position[0][1] = -10.861281772444996
Particles[6].position[0][2] = -0.1472864479578877
Particles[6].velocity[0][0] =  0.003712287568360519
Particles[6].velocity[0][1] =  0.001029937426243943
Particles[6].velocity[0][2] = -0.00004438637402367573
Particles[6].mass = 86.8103*10**24
Particles[7].position[0][0] =  9.636826813796644
Particles[7].position[0][1] = -20.860281459115482
Particles[7].position[0][2] =   0.3669173632042307
Particles[7].velocity[0][0] =   0.002954274889140831
Particles[7].velocity[0][1] =   0.001020050176036712
Particles[7].velocity[0][2] =  -0.00008868810149344669
Particles[7].mass = 102.41*10**24
Particles[8].position[0][0] =   0.002293318515953204
Particles[8].position[0][1] =   0.004296130319630977
Particles[8].position[0][2] =  -0.00007649869320601100
Particles[8].velocity[0][0] =  -0.000001347282947998592 
Particles[8].velocity[0][1] =   0.000005517069296256866
Particles[8].velocity[0][2] =   0.00000003907162561092452
Particles[8].mass = 1.988544*10**30

#Leap Frog Method for the rest of the timesteps.
def leap_frog():
    for i in range(0,N):
               Particles[i].position[2][0] = Particles[i].position[0][0] + 2*timestep*Particles[i].velocity[1][0]
               Particles[i].position[2][1] = Particles[i].position[0][1] + 2*timestep*Particles[i].velocity[1][1]
               Particles[i].position[2][2] = Particles[i].position[0][2] + 2*timestep*Particles[i].velocity[1][2]    
               Particles[i].velocity[2][0] = Particles[i].velocity[0][0] + 2*timestep*Particles[i].acceleration[0]
               Particles[i].velocity[2][1] = Particles[i].velocity[0][1] + 2*timestep*Particles[i].acceleration[1]
               Particles[i].velocity[2][2] = Particles[i].velocity[0][2] + 2*timestep*Particles[i].acceleration[2]
               i = i + 1
#goes through and finds the accelerations due to each Particle on each other.
#based soley on the posisitons of each Partilce
def compute_grav_accel():
           for j in range(0,N):
            for i in range(j+1, N):
                    dx = Particles[i].position[1][0]-Particles[j].position[1][0]
                    dy = Particles[i].position[1][1]-Particles[j].position[1][1]
                    dz = Particles[i].position[1][2]-Particles[j].position[1][2]
                    r_3 = (dx**2+dy**2+dz**2)**(1.5)
                    dx = dx/r_3
                    dy = dy/r_3
                    dz = dz/r_3
                    Particles[j].acceleration[0] = Particles[j].acceleration[0] + G*Particles[i].mass*dx
                    Particles[i].acceleration[0] = Particles[i].acceleration[0] - G*Particles[j].mass*dx
                    Particles[j].acceleration[1] = Particles[j].acceleration[1] + G*Particles[i].mass*dy
                    Particles[i].acceleration[1] = Particles[i].acceleration[1] - G*Particles[j].mass*dy
                    Particles[j].acceleration[2] = Particles[j].acceleration[2] + G*Particles[i].mass*dz
                    Particles[i].acceleration[2] = Particles[i].acceleration[2] - G*Particles[j].mass*dz                         
                    i = i + 1
                    
                        
            j = j + 1
def reset_accelerations():
    for i in range(0,N):
        Particles[i].acceleration[0] = 0
        Particles[i].acceleration[1] = 0 
        Particles[i].acceleration[2] = 0
        i = i + 1
#adds the positions for x,y,z into arrays for plot    
def position_to_array():
            for i in range(0,N):        
                Particles[i].position_array_x.insert(1,Particles[i].position[1][0])
                Particles[i].position_array_y.insert(1,Particles[i].position[1][1])
                Particles[i].position_array_z.insert(1,Particles[i].position[1][2])
                i = i + 1 
def swap_positions():
           for i in range(0,N):
                #old positions/velocitys become the current   
                Particles[i].position[0][0] = Particles[i].position[1][0]
                Particles[i].position[0][1] = Particles[i].position[1][1]
                Particles[i].position[0][2] = Particles[i].position[1][2]
                Particles[i].velocity[0][0] = Particles[i].velocity[1][0]
                Particles[i].velocity[0][1] = Particles[i].velocity[1][1]
                Particles[i].velocity[0][2] = Particles[i].velocity[1][2]           
               #current positions/velocitys become the new
                Particles[i].position[1][0] = Particles[i].position[2][0]
                Particles[i].position[1][1] = Particles[i].position[2][1]
                Particles[i].position[1][2] = Particles[i].position[2][2]
                Particles[i].velocity[1][0] = Particles[i].velocity[2][0]
                Particles[i].velocity[1][1] = Particles[i].velocity[2][1]
                Particles[i].velocity[1][2] = Particles[i].velocity[2][2]
    
                i = i + 1
#Euler Method for the fist timestep.
def euler_method():
    for i in range(0,N):
            Particles[i].position[1][0] = Particles[i].position[0][0] + timestep*Particles[i].velocity[0][0]
            Particles[i].position[1][1] = Particles[i].position[0][1] + timestep*Particles[i].velocity[0][1]
            Particles[i].position[1][2] = Particles[i].position[0][2] + timestep*Particles[i].velocity[0][2]
            Particles[i].velocity[1][0] = Particles[i].velocity[0][0] + timestep*Particles[i].acceleration[0]
            Particles[i].velocity[1][1] = Particles[i].velocity[0][1] + timestep*Particles[i].acceleration[1]
            Particles[i].velocity[1][2] = Particles[i].velocity[0][2] + timestep*Particles[i].acceleration[2]
            i = i + 1


# Sets the first value of the position arrays with the initial positions
for i in range(0,N):        
            Particles[i].position_array_x.insert(0,Particles[i].position[0][0])
            Particles[i].position_array_y.insert(0,Particles[i].position[0][1])
            Particles[i].position_array_z.insert(0,Particles[i].position[0][2])
            i = i + 1
for j in range(0,N):
    for i in range(j+1, N):
        dx = Particles[i].position[0][0]-Particles[j].position[0][0]
        dy = Particles[i].position[0][1]-Particles[j].position[0][1]
        dz = Particles[i].position[0][2]-Particles[j].position[0][2]
        r_3 = (dx**2+dy**2+dz**2)**(1.5)
        dx = dx/r_3
        dy = dy/r_3
        dz = dz/r_3
        Particles[j].acceleration[0] = Particles[j].acceleration[0] + G*Particles[i].mass*dx
        Particles[i].acceleration[0] = Particles[i].acceleration[0] - G*Particles[j].mass*dx
        Particles[j].acceleration[1] = Particles[j].acceleration[1] + G*Particles[i].mass*dy
        Particles[i].acceleration[1] = Particles[i].acceleration[1] - G*Particles[j].mass*dy
        Particles[j].acceleration[2] = Particles[j].acceleration[2] + G*Particles[i].mass*dz
        Particles[i].acceleration[2] = Particles[i].acceleration[2] - G*Particles[j].mass*dz                         
        i = i + 1                    
    j = j + 1
while counter_1 < timestep:
        euler_method()        
        #adds the positions for x,y,z into arrays for plot
        for i in range(0,N):        
                Particles[i].position_array_x.insert(1,Particles[i].position[1][0])
                Particles[i].position_array_y.insert(1,Particles[i].position[1][1])
                Particles[i].position_array_z.insert(1,Particles[i].position[1][2])
                i = i + 1              
        counter_1 = counter_1 + 2      
while counter_1 > timestep:
    for counter_2 in range(1,52):       
       reset_accelerations() 
       compute_grav_accel() 
       leap_frog() 
       position_to_array() 
       swap_positions()
       counter_2 = counter_2 + 1
    counter_1 = counter_1 - 5                      
              

plot_1 = plt.plot(Particles[0].position_array_x, Particles[0].position_array_y, '.')
plot_2 = plt.plot(Particles[1].position_array_x, Particles[1].position_array_y, '.')
plot_3 = plt.plot(Particles[2].position_array_x, Particles[2].position_array_y, '.')
plot_4 = plt.plot(Particles[3].position_array_x, Particles[3].position_array_y, '.')
plot_5 = plt.plot(Particles[4].position_array_x, Particles[4].position_array_y, '.')
plot_6 = plt.plot(Particles[5].position_array_x, Particles[5].position_array_y, '.')
plot_7 = plt.plot(Particles[6].position_array_x, Particles[6].position_array_y, '.')
plot_8 = plt.plot(Particles[7].position_array_x, Particles[7].position_array_y, '.')
plot_9 = plt.plot(Particles[8].position_array_x, Particles[8].position_array_y, '*')
plt.show(plot_1)
plt.show(plot_2)
plt.show(plot_3)
plt.show(plot_4)
plt.show(plot_5)
plt.show(plot_6)
plt.show(plot_7)
plt.show(plot_8)
plt.show(plot_9)


