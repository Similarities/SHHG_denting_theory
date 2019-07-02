__author__ = 'julia'

import numpy as np
import math
import matplotlib.pyplot as plt
import itertools

# a0
# Ln - length of plasma gradient
# Pi = sqrt(RZmecosTheta/2Amp)
c = 3*10E8 # speed of light
me = 0.511 *10E6/(c**2) #mass electron in MeV*s**2/(m**2)

mp = 1883*me
# Z = charge state
# A = mass number (average to be given)
# R = plasma reflectin coefficent

# Theta = inc angle in rad 0.79rad for 45degree
class denting:
    def __init__(self,a0, Ln, R, Z, A, Theta):
        self.a0 = a0
        self.Ln = Ln
        self.R = R
        self.Z = Z
        self.A = A
        self.Theta = Theta

        self.c = 3*1E8
        self.me = 0.511 *1E6/(self.c**2)
        self.mp = 1883 *self.me

        self.Pi = float
        self.tau = 20
        self.w0 = 6*1E-6
        self.lambdaL = 800*1E-9


    def coefficient(self):
        Pi = (self.R * self.Z * self.me * math.cos(self.Theta) / (2*self.A*self.mp))

        return Pi

    def denting_vs_Ln (self):
        Ln = list(range(1,100))

        xi_list=list(itertools.repeat(0, 99))
        print(xi_list, "xi list")

        print (Ln, "length of plasma gradient in nm")

        self.Pi = self.coefficient()


        for x in range(0,len(Ln)):

            xi_list[x] = 2*Ln[x]*10E-8*math.log(1+self.Pi/(2*math.cos(self.Theta))*0.5*self.a0**2)
            #print(xi_list[x], "denting as function of plasma gradient in 10E-8 m, for a0: ", self.a0)

            Ln[x] = Ln[x] *10
        self.plot_results(Ln,xi_list,"L in 10E-9m", "dening in nm")


    def plot_results(self,x,y,name_x, name_y):


        plt.plot(x[:],y[:], linewidth=2)
        plt.xlabel(name_x)
        plt.ylabel(name_y)
        #plt.legend()
        plt.show()




    def denting_vs_a0(self):
        self.Pi = self.coefficient()

        a0_list = list(range(1,10))

        xi_list=list(itertools.repeat(0, 9))

        for x in range(0,len(a0_list)):

            a0_list[x] = a0_list[x] * 0.5
            print(a0_list, "a0 in")

            xi_list[x] = 2*self.Ln*10E-8*math.log(1+self.Pi/(2*math.cos(self.Theta))*0.5*(a0_list[x]**2))


        print(a0_list, xi_list, ' outcome')
        self.plot_results(a0_list, xi_list, 'a0', 'denting in m')




    def denting_constant(self):

        self.Pi = self.coefficient()

        xi = 2*self.Ln*10E-8*math.log(1+self.Pi/(2*math.cos(self.Theta))*0.5*(self.a0**2))
        return xi




    def denting_jana_modded(self):
        D_x_t = list(range(0,20))
        t_list = list(range(0,20))
        x_coos= list(range(0,20))
        print(x_coos, "xcoos erst")
        norm_xcoos = self.w0/10
        print(norm_xcoos,  "norm to w0")


        print(t_list, 'tlist')
        xi = self.denting_constant()
        print(xi)


        for x in range(0, len(D_x_t)):

            C1= 2*(x_coos[x]*norm_xcoos)**2/(self.w0**2)
            print(C1, "c1 coefficient")
            V2 = 2*((t_list[x]**2)/(self.tau**2))
            print(V2, "v2 coefficent")

            D_x_t[x] = -1*((xi* math.exp(2*(-V2-C1)))-xi)

        self.plot_results(t_list, D_x_t, "t in fs", 'Denting per time')

    def refocus_length_of_N(self):

        R_zf = self.w0**2/(4*self.denting_constant())

        print(R_zf, "curvature in m")
        N_HHG = list(range(10,50))
        zf_N = list(range(10,50))

        for x in range(0,len(N_HHG)):
            numerator = (4*self.denting_constant()*self.w0**2)
            denumerator = ((2*self.denting_constant())**2) + self.lambdaL/((math.pi*N_HHG[x])**2)

            zf_N[x] = numerator/denumerator

        self.plot_results(N_HHG,zf_N, "harmonic order N", 'focus position after target')

    def divergence_HHG_with_denting(self):

        Theta_HHG = list(range(2,50))
        N_HHG = list(range(2,50))

        for x in range(len(N_HHG)):
            C1 = 1/(math.pi *self.w0)
            Theta_HHG[x] = C1*((4*math.pi*self.denting_constant())**2 + (self.lambdaL/N_HHG[x])**2)**0.5

            #rad?

        print("denting for the given parameters", self.denting_constant())
        self.plot_results(N_HHG, Theta_HHG, "harmonic order N", "HHG divergence in [rad?]")





#denting_vs_Ln(3.5, 0.8, 8, 20, 0.79 )

Denting=denting(0.5, 5000, 0.8, 8, 10, 0.79)
#Denting.denting_vs_Ln()
#Denting.denting_vs_a0()
#Denting.denting_jana_modded()
#Denting.refocus_length_of_N()
Denting.divergence_HHG_with_denting()