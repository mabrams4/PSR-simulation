#PSRSimulationCodeSingle

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np

# Code written during 2021 summer as a research intern at the Technion Institute of Technology
# Author: Matan Abrams

# This program simulates the conditions in a Perfectly Stirred Reactor for a ramjet combustion engine. 
# It takes as inputs an air temperature (Kelvin), inner and outer inlet diameters (mm), air flow rate (kg/s), fuel flow rate (g/s),
# outlet diameter (mm) and flame holder diameter (mm), and returns the 
# residence time, flow time, Damkohler number, and equivalence ratio of the modeled reaction. 

# Perfectly Stirred Reactor Class
class PSR:
    
    # initializes variables to None as they will be set later
    def __init__(self, solution_name):
        self.solution = ct.Solution(solution_name)
        self.temp_air = None
        self.pressure_air = None
        self.inner_d = None
        self.outer_d = None
        self.flame_holder_d = None
        self.mdot_air = None
        self.mdot_fuel = None
        self.exit_d = None
    
    # sets the outlet diameter for the reaction
    def set_exit_diameter(self):
        try:
            exit_d = float(input('Enter exit diameter (mm): ')) / 1000
            self.exit_d = exit_d
        except ValueError:
            print('Not a valid number. Please try again')
            
        return exit_d
    
    # sets the inner and outer inlet diameters for the reaction
    def set_inner_outer_inlet_diameters(self):
        try:
            inner = float(input("Enter inner diameter (mm): ")) / 1000
            outer = float(input("Enter outer diameter (mm): ")) / 1000
            self.inner_d = inner
            self.outer_d = outer
        except ValueError:
            print("Not a valid number. Please try again")

        return inner, outer
    
    # sets the flame holder diameter for the reaction
    def set_flame_holder_diameter(self):
        try:
            flame_holder_d = float(input("Enter flame holder diameter (mm): ")) / 1000
            self.flame_holder_d = flame_holder_d 
        except ValueError:
            print("Not a valid number. Please try again")

        return flame_holder_d
    
    # sets the air flow rate for the reaction
    def set_air_flow_rate(self):
        try:
            mdot_air = float(input("Enter air flow rate [kg/s]: "))
            self.mdot_air = mdot_air
        except ValueError:
            print("Not a valid number. Please try again")

        return mdot_air
           
    # sets the air temperature for the reaction       
    def set_air_temp(self):
        try:
            temp_air = float(input("Enter air temperature [K]: "))
            self.temp_air = temp_air
        except ValueError:
            print("Not a valid number. Please try again")

        return temp_air
        
    # sets the fuel flow rate for the reaction
    def set_fuel_flow_rate(self):
        try:
            mdot_fuel = float(input("Enter fuel flow rate [g/s]: ")) / 1000
            self.mdot_fuel = mdot_fuel
        except ValueError:
            print("Not a valid number. Please try again")

        return mdot_fuel
    
    # calculates the mass flow rate of the reaction
    def calculate_mdot(self, M):
        gas = self.solution
        P = self.pressure_air
        T = self.temp_air
        g = gas.cp_mass / gas.cv_mass
        R = gas.cp_mass - gas.cv_mass
        A = np.pi * ((self.outer_d / 2)**2 - (self.inner_d / 2)**2)
        
        term1 = A * P / np.sqrt(T)
        term2 = np.sqrt(g / R) * M
        term3 = (1 + ((g - 1) / 2) * M**2)
        term4 = -(g + 1) / (2 * (g - 1))
        mdot = term1 * term2 * term3**term4
        return mdot
     
    # iterates to narrow in on the correct value using the commonly known bisection method    
    def bisection_method(self):
        tol = 1e-4
        n = 0; a = 0; b = 1;
        max_iter = 100
        diff = 2 * tol
        m_real = self.mdot_air + self.mdot_fuel
        while n <= max_iter and abs(diff) > tol:
            M_guess = (a + b) / 2
            m_guess = self.calculate_mdot(M_guess)
            diff = m_real - m_guess
            if diff > 0:
                a = M_guess
            else:
                b = M_guess
            n += 1
            if n > max_iter:
                raise ValueError('reached max iterations')
                
        return M_guess 

    # calculates the air pressure of the reaction
    def calculate_air_pressure(self):
        gas = self.solution
        mdot = self.mdot_air + self.mdot_fuel
        P_calc = 8e5
        tol = 1
        max_iter = 100
        n = 0
        diff = 2 * tol
        while n <= max_iter:
            P_guess=P_calc
            gas.TPX = self.temp_air, P_guess,'O2:1.0, N2:3.76'
            air_mass_frac = gas.Y
            gas.TPX = 293, P_guess, 'POSF10264:1'
            fuel_mass_frac = gas.Y
    
            total_mass_frac = [0] * len(fuel_mass_frac)
            for i in range(len(total_mass_frac)):
                total_mass_frac[i] = (air_mass_frac[i] * self.mdot_air) + (fuel_mass_frac[i] * self.mdot_fuel)
                
            gas.TPY = self.temp_air, P_guess, total_mass_frac
            phi = gas.equivalence_ratio()
            gas.set_equivalence_ratio(phi, 'POSF10264:1', 'O2:1.0, N2:3.76')
            gas.equilibrate('HP') 
            
            A = np.pi * (self.exit_d / 2)**2
            g = gas.cp_mass / gas.cv_mass
            R = gas.cp_mass - gas.cv_mass 
            term1 = mdot * np.sqrt(gas.T) / A
            term2 = np.sqrt(R / g)
            term3 = (g + 1) / 2
            term4 = -(g + 1) / (2 * (g - 1))
            P_calc = term1 * term2 / term3**term4
            diff = abs(P_calc - P_guess)
            n += 1
            if abs(diff) < tol:
                self.pressure_air = P_calc
                return P_calc
        if n > max_iter:
            raise ValueError('reached max iterations')

    # simulates the reaction 
    def solve(self):
        inner_d, outer_d = self.set_inner_outer_inlet_diameters()
        flame_holder_d = self.set_flame_holder_diameter()
        self.set_exit_diameter()
        temp_air = self.set_air_temp()
        mdot_air = self.set_air_flow_rate()
        mdot_fuel = self.set_fuel_flow_rate()
        pressure_air = self.calculate_air_pressure()
        gas = self.solution
        
        gas.TPX = temp_air, pressure_air,'O2:1.0, N2:3.76'
        air_mass_frac = gas.Y
        gas.TPX = 293, pressure_air, 'POSF10264:1'
        fuel_mass_frac = gas.Y

        total_mass_frac = [0] * len(fuel_mass_frac)
        for i in range(len(total_mass_frac)):
            total_mass_frac[i] = (air_mass_frac[i] * mdot_air) + (fuel_mass_frac[i] * mdot_fuel)
            
        gas.TPY = temp_air, pressure_air, total_mass_frac
        phi = gas.equivalence_ratio()
        gas.set_equivalence_ratio(phi, 'POSF10264:1', 'O2:1.0, N2:3.76')
        inlet = ct.Reservoir(gas)
        
        g = gas.cp_mass / gas.cv_mass
        R = gas.cp_mass - gas.cv_mass        
        M = self.bisection_method()
        static_temp = self.temp_air / (1 + ((g - 1) / 2) * M**2)
        a = np.sqrt(g * R * static_temp)
        V = M * a
        t_flow = flame_holder_d / V
        
        gas.equilibrate('HP')  
        exhaust = ct.Reservoir(gas)
        reactor = ct.IdealGasReactor(gas)
        reactor.volume = 1.0
        temp_add = gas.T
        theta = (reactor.T - temp_air) / (temp_add - temp_air)
        
        def mdot(t):
            return reactor.mass / residence_time
        
        residence_time = 0.1  # starting residence time
        inlet_mfc = ct.MassFlowController(inlet, reactor, mdot=mdot)       
        outlet_mfc = ct.PressureController(reactor, exhaust, master=inlet_mfc, K=0.01)
        sim = ct.ReactorNet([reactor])
        states = ct.SolutionArray(gas, extra=['tres']) 

        theta_arr = []
        while theta > 0.5:  
            sim.set_initial_time(0.0)  # reset the integrator
            sim.advance_to_steady_state()
            theta = (reactor.T - temp_air) / (temp_add - temp_air) * 100
            theta_arr.append(theta)
            states.append(reactor.thermo.state, tres=residence_time)
            residence_time *= 0.9  # decrease the residence time for the next iteration
            
        # Plot results
        f, ax1 = plt.subplots(1, 1)
        ax1.plot(states.tres, states.heat_release_rate, '.-', color='C0')
        ax2 = ax1.twinx()
        ax1.set_xscale('log')
        ax2.plot(states.tres[:-1], theta_arr[:-1], '.-', color='C1')
        ax1.set_xlabel('residence time [s]')
        ax1.set_ylabel('heat release rate [W/m$^3$]', color='C0')
        ax2.set_ylabel('theta [%]', color='C1')
        f.tight_layout()
        plt.show()
        
        p_bar = pressure_air / 1e5
        mdot_fuel_gs = mdot_fuel * 1000
        print(f'P = {round(p_bar , 2)} [bar] T = {round(temp_air)} [K] ' 
              f'mdot_air = {round(mdot_air)} [kg/s] mdot_fuel = {round(mdot_fuel_gs)} [g/s] '
              f'T_adiabatic = {round(temp_add)} [K]')
          
        return residence_time, t_flow, phi
    
# runs the program based on sim which encodes the various properties of the reacting gas and prints out the results    
if __name__ == "__main__":
    sim = PSR("HyChemHighTdetailed.cti")
    t_res, t_flow, phi = sim.solve()
    Da = t_flow / t_res
    t_res *= 1e6
    t_flow *= 1000
    print(f'residence_time = {round(t_res, 2)} [microseconds]')
    print(f't_flow = {round(t_flow, 2)} [ms]')
    print(f'Da = {round(Da, 2)}; equiv_ratio = {round(phi, 2)}')
    



