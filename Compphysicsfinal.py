#Earth and Moon orbiting the Sun
#
#
#Bobby Flores
#April 30,2014


from math import *
from pylab import *
import numpy

"""This model is performed in units of astronomical units(AU) and years(yr),
which is the distance from the Earth to the Sun and the time it takes for
the Earth to orbit the Sun, respectively. This ensures that the numbers
distance, time, and velocity values are much small and more managable.
In yr/AU units 4*pi**2 is equal to the gravitational constant G multiplied
by the mass of the Sun, which is used frequently in the differential equations."""

"""For variables e=Earth,m=Moon,s=Sun,vx=velocity in x direction,vy=velocity in
y direction,M=mass,r=distance"""

#Mass of the Earth, Moon, and Sun
Me = 5.97e24
Mm = 7.35e22
Ms = 1.99e30

#G*Ms in yr/AU units
C=4*pi**2

"""The following functions are differential equations that represent the
acceleration in the x and y directions of the earth and moon due to the gravitional
forces between the Earth, Moon, and Sun. Each of the equations combines the
acceleration terms due to the other two respective celestial bodies. All four
functions take the input parameters of the component position (x or y) of the
Earth and Moon and the distance between the Earth in moon at each given time.
Additionally, functions for Earth's acceleration takes the parameter of Earth's
orbit radius around the Sun and the functions for the Moon take the Moon's orbit
radius around the Sun."""
def Aex(xe,xm,re,rem):
    """Acceleration of the Earth in the x direction
Ax = d**2x/dt**2 = dVx/dt = -(G*Ms*xe)/re**3-(G*Mm*(xe-xm))/rem**3 """
    return -(C*xe)/re**3-(C*(Mm/Ms)*(xe-xm))/rem**3

def Aey(ye,ym,re,rem):
    """Acceleration of the Earth in the y direction
Ay = d**2y/dt**2 = dVy/dt = -(G*Ms*ye)/re**3-(G*Mm*(ye-ym))/rem**3 """
    return -(C*ye)/re**3-(C*(Mm/Ms)*(ye-ym))/rem**3

def Amx(xe,xm,rm,rem):
    """Acceleration of the Moon in the x direction
Ax = d**2x/dt**2 = dVx/dt = -(G*Ms*xm)/rm**3-(G*Me*(xm-xe))/rem**3"""
    return -(C*xm)/rm**3-(C*(Me/Ms)*(xm-xe))/rem**3

def Amy(ye,ym,rm,rem):
    """Acceleration of the Moon in the x direction
Ay = d**2y/dt**2 = dVy/dt = -(G*Ms*ym)/rm**3-(G*Me*(ym-ye))/rem**3"""
    return -(C*ym)/rm**3-(C*(Me/Ms)*(ym-ye))/rem**3

##Intial conditions
#Both the Earth and Moon are placed at their maximum radius and minimum velocity
#initial velocity and position of Earth
vxei = 0
vyei = 6.17871
xei = 1.0167
yei = 0
#intial velocity and position of Moon
vxmi = 0
vymi = vyei+.203355
xmi = xei+.0027106+.000061
ymi = 0

time = 1
steps = 10000

def earth_moon_orbit(vxei,vyei,xei,yei,vxmi,vymi,xmi,ymi, time,steps):
    """This function tracks the position and velocity of the earth as it
orbits the sun and the moon as it orbits the earth and sun. This function
assumes that the Sun is at rest and does not take into account the effects
of the other planets. The orbit is tracked using the Euler-Cromer method for
solving differential equations. In this case, the method is used to solve the
ordinary differential equations for the x and y components of the acceleration
due to gravity. This function takes the parameters for the initial position of
the earth and moon, initial velocities or the earth and moon,time duration,
and number of steps."""
    #position and velocity of earth
    vxe = [0]*(steps+1)
    vye = [0]*(steps+1)
    xe = [0]*(steps+1)
    ye = [0]*(steps+1)
    vxe[0] = vxei
    vye[0] = vyei
    xe[0] = xei
    ye[0] = yei
    re = [0]*(steps+1)
    ve = [0]*(steps+1)

    #position and velocity of moon
    vxm = [0]*(steps+1)
    vym = [0]*(steps+1)
    xm = [0]*(steps+1)
    ym = [0]*(steps+1)
    vxm[0] = vxmi
    vym[0] = vymi
    xm[0] = xmi
    ym[0] = ymi
    rm = [0]*(steps+1)
    vm = [0]*(steps+1)
    
    tstep=(time/steps)
    rme=[0]*(steps+1)

    #Euler-Cromer Method of integration
    for i in range(steps):
        #distance of earth to sun, moon to sun, and earth to moon at each timestep
        re[i] = sqrt(xe[i]**2+ye[i]**2)
        rm[i] = sqrt(xm[i]**2+ym[i]**2)
        rme[i] = sqrt((xe[i]-xm[i])**2+(ye[i]-ym[i])**2)

        #total velocity of Earth and Moon
        ve[i] = sqrt(vxe[i]**2+vye[i]**2)
        vm[i] = sqrt(vye[i]**2+vye[i]**2)

        #components of velocity and postion of earth at each timestep
        vxe[i+1] = vxe[i]+Aex(xe[i],xm[i],re[i],rme[i])*tstep
        xe[i+1] = xe[i]+vxe[i+1]*tstep
        vye[i+1] = vye[i]+Aey(ye[i],ym[i],re[i],rme[i])*tstep
        ye[i+1] = ye[i]+vye[i+1]*tstep
        
        #velocity and position of moon
        vxm[i+1] = vxm[i]+Amx(xe[i],xm[i],rm[i],rme[i])*tstep
        xm[i+1] = xm[i]+vxm[i+1]*tstep
        vym[i+1] = vym[i]+Amy(ye[i],ym[i],rm[i],rme[i])*tstep
        ym[i+1] = ym[i]+vym[i+1]*tstep
        
        if re[i]<(1/214):
            #the the radius is smaller than the radius of the sun
            return "The Earth crashed into the Sun!"

    #Eccentricities of Earth around the Sun and Moon around the Earth
    em=(max(rme)-min(rme[:-1]))/(max(rme)+min(rme[:-1]))
    ee=(max(re)-min(re[:-1]))/(max(re)+min(re[:-1]))
    print("Moon",em,max(rme),min(rme[:-1]))
    print("Earth",ee,max(re),min(re[:-1]))

    #Path of Earth and Moon around the Sun
    plot(xe,ye,"b")
    plot(xm,ym,"r")
    xlabel("x")
    ylabel("y")
    title("Earth and Moon path around the Sun")
    legend(("Earth","Moon"))
    xlim(.5,.7)
    ylim(.7,.9)
    show()

    #position of Earth and Moon around the Sun at each time interval
    Xe = [0]*(steps+1)
    Ye = [0]*(steps+1)
    Xm = [0]*(steps+1)
    Ym = [0]*(steps+1)
    for i in range(steps+1): #1/20 of calculate position values are plotted
        if i%20==0:
            Xe[i] = xe[i]
            Ye[i] = ye[i]
            Xm[i] = xm[i]
            Ym[i] = ym[i]        
    plot(Xe,Ye,"bo")
    plot(Xm,Ym,"ro")
    xlabel("x")
    ylabel("y")
    xlim(.95,1.025)
    ylim(0,.3)
    legend(("Earth","Moon"))
    title("Position of Earth and Moon around Sun")
    show()
    
    
    #distance between the earth and moon as a function of time
    t=linspace(0,1,steps+1)
    plot(t,rme)
    ylim(.002,.003)
    xlabel("Time (yrs)")
    ylabel("Distance (AU)")
    title("Distance between the earth and moon as a function of time")
    show()
    
    #orbit distance of earth and moon around the sun as a function of time
    plot(t,re,"b")
    plot(t,rm,"r")
    ylim(.97,1.03)
    xlabel("Time (yrs)")
    ylabel("Distance (AU)")
    legend(("Earth","Moon"))
    title("Orbit distance of Earth and Moon around the Sun as a function of time")
    show()

    
    
print(earth_moon_orbit(vxei,vyei,xei,yei,vxmi,vymi,xmi,ymi, time,steps))

