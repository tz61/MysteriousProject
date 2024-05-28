from matplotlib import pyplot as plt
import numpy as np
def IVP1func(t, y):
    return y**2-t**2
IVP1_t0 = 0
IVP1_y0 = 1
def IVP2func(t,y):
    return y**3-t**2*y
IVP2_t0 = 0
IVP2_y0 = 1
# alias for methods
EULER_METHOD = 0
IMP_EULER_METHOD = 1
RUNGE_4TH_METHOD = 2
ABM_METHOD = 3
MS_METHOD = 4
HAMMING_METHOD = 5
POWER_SERIES_METHOD = 6
# alias for IVP
IVP1 = 0
IVP2 = 1
Iter_Method=['Euler Method','Improved Euler Method','Runge-Kutta 4th Order Method','Adams-Bashforth-Moulton Method','Milne-Simpson Method','Hamming Method','Power Series Method']
IVP=[[IVP1func,IVP1_t0,IVP1_y0],[IVP2func,IVP2_t0,IVP2_y0]]
class NumericalSols:
    def __init__(self, IVPQuestion, step, lower_bound, upper_bound, method_num, overflow_cap=50):
        self.t_0 = IVP[IVPQuestion][1]
        self.y_0 = IVP[IVPQuestion][2]
        self.IVP_func = IVP[IVPQuestion][0]
        assert(lower_bound <= 0 and upper_bound >= 0)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.y_overflow_cap = overflow_cap
        self.step = step
        self.t_all = []
        self.y_all = []
        self.method_num = method_num
        print("start IVP%d with %s" % (IVPQuestion+1, Iter_Method[method_num]))
        if method_num == 0:
            self.eulerMethod()
        elif method_num == 1:
            self.impEuler()
        elif method_num == 2:
            self.runge4thMethod()
        elif method_num == 3:
            self.ABM()
        elif method_num == 4:
            self.MS()
        elif method_num == 5:
            self.Hamming()        
        elif method_num == 6:
            self.powerSeries()
        else:
            raise NotImplementedError
        print("end IVP%d with %s" % (IVPQuestion+1, Iter_Method[method_num]))

    def assertOverflow(self, t_n, y_n):
        if abs(y_n) > self.y_overflow_cap:
            raise OverflowError("Discontinuous at t=%.3f, y=%.3f" % (t_n, y_n))
    # method_num: 0
    def eulerMethod(self):
        # less than 0
        i = 1
        y_n = t_n = 0
        while t_n >= self.lower_bound:
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0
            
            self.y_all.insert(0, y_n)
            self.t_all.insert(0, t_n)
            try:
                y_n = y_n + self.IVP_func(t_n, y_n) * (-self.step)
                t_n = t_n + (-self.step)
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break
        
        # larger than 0   
        i = 1
        y_n = t_n = 0
        while t_n <= self.upper_bound:
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0
            self.y_all.append(y_n)
            self.t_all.append(t_n)
            try:
                y_n = y_n + self.IVP_func(t_n, y_n) * self.step
                t_n = t_n + self.step
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break

    # 1
    def impEuler(self):
        i = 1
        y_n = t_n = 0
        while t_n >= self.lower_bound:
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0
            self.y_all.insert(0, y_n)
            self.t_all.insert(0, t_n)
            # go left in t
            try:
                y_n_euler = y_n + self.IVP_func(t_n, y_n) * (-self.step)    
                y_n = y_n + 0.5 * (-self.step) * (self.IVP_func(t_n, y_n) + self.IVP_func(t_n + (-self.step), y_n_euler))
                t_n = t_n + (-self.step)
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break
        
        i = 1
        y_n = t_n = 0
        while t_n <= self.upper_bound:
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0
            self.y_all.append(y_n)
            self.t_all.append(t_n)
            # go right in t
            try:
                y_n_euler = y_n + self.IVP_func(t_n, y_n) * self.step
                y_n = y_n + 0.5 * self.step * (self.IVP_func(t_n, y_n) + self.IVP_func(t_n + self.step, y_n_euler))
                t_n = t_n + self.step
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break
    
    # 2
    def runge4thMethod(self):
        # go left in t
        i = 1
        y_n = t_n = 0
        while t_n >= self.lower_bound:
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0
            self.y_all.insert(0, y_n)
            self.t_all.insert(0, t_n)
            try:
                k1 = self.step * self.IVP_func(y_n, t_n)
                k2 = self.step * self.IVP_func(y_n + 0.5 * k1, t_n + 0.5 * self.step)
                k3 = self.step * self.IVP_func(y_n + 0.5 * k2, t_n + 0.5 * self.step)
                k4 = self.step * self.IVP_func(y_n + k3, t_n + self.step)
                y_n = y_n - (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
                t_n = t_n + (-self.step)
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break
        # go right in t
        i = 1
        y_n = t_n = 0
        while t_n <= self.upper_bound:
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0
            self.y_all.append(y_n)
            self.t_all.append(t_n)
            try:
                k1 = self.step * self.IVP_func(y_n, t_n)
                k2 = self.step * self.IVP_func(y_n + 0.5 * k1, t_n + 0.5 * self.step)
                k3 = self.step * self.IVP_func(y_n + 0.5 * k2, t_n + 0.5 * self.step)
                k4 = self.step * self.IVP_func(y_n + k3, t_n + self.step)
                y_n = y_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
                t_n = t_n + self.step
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break
    
    ### method 3, 4, 5 are linear multiple step methods, and the code below just implemented the case with t>=0
    # 3
    def ABM(self):
        y_n = self.y_0
        t_n = self.t_0
        i=0
        while i<4:
            self.y_all.insert(0, y_n)
            self.t_all.insert(0, t_n)
            try:
                k1 = self.step * self.IVP_func(y_n, t_n)
                k2 = self.step * self.IVP_func(y_n + 0.5 * k1, t_n + 0.5 * self.step)
                k3 = self.step * self.IVP_func(y_n + 0.5 * k2, t_n + 0.5 * self.step)
                k4 = self.step * self.IVP_func(y_n + k3, t_n + self.step)
                y_n = y_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
                t_n = t_n + self.step
                i=i+1
                self.assertOverflow(t_n, y_n)
            except OverflowError as e:
                print(e)
                break
        k=3
        while t_n <= self.upper_bound:
            try:
                t_0=self.t_all[k-3]
                t_1=self.t_all[k-2]
                t_2=self.t_all[k-1]
                t_3=self.t_all[k]
                y_0=self.y_all[k-3]
                y_1=self.y_all[k-2]
                y_2=self.y_all[k-1]
                y_3=self.y_all[k]
                dy=y_3+(self.step/24)*(-9*self.IVP_func(t_0,y_0)+37*self.IVP_func(t_1,y_1)-59*self.IVP_func(t_2,y_2)+55*self.IVP_func(t_3,y_3))
                t_n=t_3+self.step
                y_n=y_3+(self.step/24)*(self.IVP_func(t_1,y_1)-5*self.IVP_func(t_2,y_2)+19*self.IVP_func(t_3,y_3)+9*self.IVP_func(t_n,dy))
                self.assertOverflow(t_n, y_n)
                self.t_all.append(t_n)
                self.y_all.append(y_n)
                k=k+1
            except OverflowError as e:
                print(e)
                break
            
    # 4
    def MS(self):
        y_n = self.y_0
        t_n = self.t_0
        i=0     
        while i<4:
            self.y_all.insert(0, y_n)
            self.t_all.insert(0, t_n)
            k1 = self.step * self.IVP_func(y_n, t_n)
            k2 = self.step * self.IVP_func(y_n + 0.5 * k1, t_n + 0.5 * self.step)
            k3 = self.step * self.IVP_func(y_n + 0.5 * k2, t_n + 0.5 * self.step)
            k4 = self.step * self.IVP_func(y_n + k3, t_n + self.step)
            y_n = y_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            t_n = t_n + self.step
            i=i+1
        k=3
        while (t_n <= self.upper_bound):
            t_0=self.t_all[k-3]
            t_1=self.t_all[k-2]
            t_2=self.t_all[k-1]
            t_3=self.t_all[k]
            y_0=self.y_all[k-3]
            y_1=self.y_all[k-2]
            y_2=self.y_all[k-1]
            y_3=self.y_all[k]
            dy=y_0+4*(self.step/3)*(2*self.IVP_func(t_1,y_1)-self.IVP_func(t_2,y_2)+2*self.IVP_func(t_3,y_3))
            t_n=t_3+self.step
            y_n=y_2+(self.step/3)*(self.IVP_func(t_2,y_2)+4*self.IVP_func(t_3,y_3)+self.IVP_func(t_n,dy))
            if(abs(y_n)>5): break
            self.t_all.append(t_n)
            self.y_all.append(y_n)
            k=k+1
            
    # 5
    def Hamming(self):
        y_n = self.y_0
        t_n = self.t_0
        i=0     
        while i<4:
            self.y_all.insert(0, y_n)
            self.t_all.insert(0, t_n)
            k1 = self.step * self.IVP_func(y_n, t_n)
            k2 = self.step * self.IVP_func(y_n + 0.5 * k1, t_n + 0.5 * self.step)
            k3 = self.step * self.IVP_func(y_n + 0.5 * k2, t_n + 0.5 * self.step)
            k4 = self.step * self.IVP_func(y_n + k3, t_n + self.step)
            y_n = y_n + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            t_n = t_n + self.step
            i=i+1
        k=3
        while (t_n <= self.upper_bound):
            t_0=self.t_all[k-3]
            t_1=self.t_all[k-2]
            t_2=self.t_all[k-1]
            t_3=self.t_all[k]
            y_0=self.y_all[k-3]
            y_1=self.y_all[k-2]
            y_2=self.y_all[k-1]
            y_3=self.y_all[k]
            dy=y_0+4*(self.step/3)*(2*self.IVP_func(t_1,y_1)-self.IVP_func(t_2,y_2)+2*self.IVP_func(t_3,y_3))
            t_n=t_3+self.step
            y_n=-y_1/8+9*y_3/8+3*(self.step/8)*(-self.IVP_func(t_2,y_2)+2*self.IVP_func(t_3,y_3)+self.IVP_func(t_n,dy))
            if(abs(y_n)>5): break
            self.t_all.append(t_n)
            self.y_all.append(y_n)
            k=k+1  

    # 6
    def powerSeries(self):
        # pass
        a = np.zeros(10000, dtype=float)
        b = np.zeros(10000, dtype=float)
        c = np.zeros(10000, dtype=float)
        
        n_pwr = 10
        # manully calculate
        a[0] = self.y_0
        self.help_bc(a, b, c, 0)
        a[1] = c[0]
        self.help_bc(a, b, c, 1)
        a[2] = (c[1]-a[0])/2
        self.help_bc(a, b, c, 2)
        a[3] = (c[2]-b[0]-a[1])/3 
        self.help_bc(a, b, c, 3)
        a[4] = (c[3]-b[1]-a[2]+1)/4 

        for i in range(4,n_pwr):
            self.help_bc(a, b, c, i)
            a[i+1] = (c[i]-b[i-2]-a[i-1])/(i+1)

        for i in range(n_pwr):    
            print(i, "a=", a[i], "b=", b[i], "c=", c[i])
        self.help_radius(a, n_pwr)

        self.help_PwrSrs(a, b, c, -1, n_pwr)
        self.help_PwrSrs(a, b, c, 1, n_pwr)
    
    def help_bc(self, a, b, c, i):
        for j in range(i+1):
            b[i] += a[j]*a[i-j]
        for j in range(i+1):
            c[i] += a[j]*b[i-j]

    def help_PwrSrs(self, a, b, c, dir, n):
        i = 1
        y_n = t_n = 0
        while (t_n>=self.lower_bound and t_n<=self.upper_bound):
            if i == 1:
                i = 0
                y_n = self.y_0
                t_n = self.t_0

            if(dir==1):
                self.y_all.append(y_n)
                self.t_all.append(t_n)
            elif(dir==-1):
                self.y_all.insert(0, y_n)
                self.t_all.insert(0, t_n)
                
            y_n = 0
            for i in range(n):
                y_n += a[i] * (t_n ** i)
            if(abs(y_n)>50): break

            t_n = t_n + dir* self.step

    def help_radius(self, a, n):
        for i in range(n-1):
            # print(i, "r=", abs(a[i]/a[i+1]))
            pass

    # Draw
    def draw(self,color_="red",label_=""):
        plt.xlabel('t', fontsize=19)
        plt.ylabel('y', fontsize=19)
        plt.plot(self.t_all, self.y_all, c=color_, label=label_)

    def drawlinear(self,n):
        plt.xlabel('t', fontsize=19)
        plt.ylabel('y', fontsize=19)
        if (n==3):
            plt.plot(self.t_all, self.y_all, c="red",label="ABM")
        if (n==4):
            plt.plot(self.t_all, self.y_all, c="blue", label="MS")
        if (n==5):
            plt.plot(self.t_all, self.y_all, c="green", label="Hamming")
def testDiffMethods(IVPQuestion, step, lower_bound, upper_bound):
    """
    Test different methods according to the given parameters by drawing plots.
    """
    # eulerMethod
    euler_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 0)
    plt.plot(euler_meth.t_all, euler_meth.y_all, c="green", label="Euler method")
    
    # impEuler
    im_eu_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 1)
    plt.plot(im_eu_meth.t_all, im_eu_meth.y_all, c="red", label="Improve Euler method")
    
    # runge4thMethod
    runge_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 2)
    plt.plot(runge_meth.t_all, runge_meth.y_all, c="blue", label="Runge 4th order")
    
    # ABM
    mulstep_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 3)
    plt.plot(mulstep_meth.t_all, mulstep_meth.y_all, c="brown", label="ABM")
    
    # MS
    mulstep_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 4)
    plt.plot(mulstep_meth.t_all, mulstep_meth.y_all, c="black", label="MS")

    # Hamming
    mulstep_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 5)
    plt.plot(mulstep_meth.t_all, mulstep_meth.y_all, c="purple", label="Hamming")

    
    # powerSeries method
    power_meth = NumericalSols(IVPQuestion, step, lower_bound, upper_bound, 6)
    plt.plot(power_meth.t_all, power_meth.y_all, c="orange", label="Power series")
def testDiffSteps(IVPQuestion, lower_bound, upper_bound, method_num):
    """
    Test one method with different steps by drawing plots.
    """
    # step == 0.2
    meth1 = NumericalSols(IVPQuestion, 0.2, lower_bound, upper_bound, method_num)
    plt.plot(meth1.t_all, meth1.y_all, c="green", label="step=0.2")
    
    # step == 0.1
    meth2 = NumericalSols(IVPQuestion, 0.1, lower_bound, upper_bound, method_num)
    plt.plot(meth2.t_all, meth2.y_all, c="red", label="step=0.1")
    
    # step == 0.05
    meth3 = NumericalSols(IVPQuestion, 0.05, lower_bound, upper_bound, method_num)
    plt.plot(meth3.t_all, meth3.y_all, c="blue", label="step=0.05")
def plotDiffMethods(IVPQuestion, step, lower_bound, upper_bound):
    # fig=plt.figure(num=1,figsize=(8,6))
    plt.xlabel('t', fontsize=19)
    plt.ylabel('y', fontsize=19)

    # ax1 = fig.add_subplot(311)
    # ax1.set_title("311")
    testDiffMethods(IVPQuestion, step, lower_bound, upper_bound)

    plt.legend()
    plt.show()
def plotDiffSteps(IVPQuestion, lower_bound, upper_bound, method_num):
    plt.xlabel('t', fontsize=19)
    plt.ylabel('y', fontsize=19)

    testDiffSteps(IVPQuestion, lower_bound, upper_bound, method_num)

    plt.legend()
    plt.show()  