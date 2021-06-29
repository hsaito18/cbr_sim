"""
Four stroke engine cycle simulator created for PHYS 203 final project. Uses simple cooling model to
track temperature increease over time in a race.
Uses "track_data2.csv" that contains RPM, coolant temperature, and oil temperature data.

Author: Henrique Saito
"""

import matplotlib.pyplot as plt
import math
import pandas as pd

# Constants
OPTIMAL_AIR_FUEL_RATIO = 14.7
OCTANE_MOLAR_MASS = 114.23
AIR_MOLAR_MASS = 28.97
AIR_C = 29.19
ENTHALPY_COMBUSTION_OCTANE = 5430000
CV_CO2 = 28.46
CV_N2 = 20.8
a = 0.137
ao2 = 0.1382
aco2 = 0.3640
b = 0.0000387
bco2 = 0.00004267
bo2 = 0.00003186
R = 8.3145
ENGINE_MASS = 5
ENGINE_HEAT_CAPACITY = 910
ENGINE_K = 230
ENGINE_THICKNESS = 0.02

"""
Kinematics of piston and crankshaft. Just to plot the P-V graph in the correct shape.
Largely inspired by https://skill-lync.com/projects/otto-cycle-simulator-using-python-8\
"""
def engine_kinematics(bore, stroke, con_rod, cr, start_crank, end_crank):
    # geometric parameters
    r = con_rod / (stroke / 2)

    # volume
    v_s = (math.pi) / 4 * pow(bore, 2) * stroke
    v_c = v_s / (cr - 1)
    sc = math.radians(start_crank)
    ec = math.radians(end_crank)
    n_theta = 100
    dtheta = (ec - sc) / (n_theta - 1)
    vol = []
    for k in range(0, n_theta):
        theta = sc + (k * dtheta)
        term1 = 0.5 * (cr - 1)
        term2 = r + 1 - math.cos(theta)
        term3 = pow(pow(r, 2) - pow(math.sin(theta), 2), 0.5)
        vol.append((1 + term1 * (term2 - term3)) * v_c)
    return vol

"""
Simulate one four stroke engine cycle, and return the final engine temperature.
Some ideal-gas elements were inspired by https://skill-lync.com/student-projects/otto-cycle-simulator-using-python-8
"""
def four_stroke(p1, t1, rps, speed, cool_temp):
    if rps == 0:
        return [t1, p1]
    time_diff = 1.0 / rps
    gamma = 7.0 / 5.0
    # geometric parameters
    bore = 0.067
    stroke = 0.0425
    con_rod = 0.0918
    cr = 12.2

    # volume
    v_s = math.pi / 4 * pow(bore, 2) * stroke
    v_c = v_s / (cr - 1)
    v1 = v_s + v_c
    v2 = v_c
    sa = 2 * math.pi / 4 * bore ** 2 * math.pi * bore * stroke / 2

    # number of moles of air
    n = p1 * v1 / (R * t1)

    # separating composition of air
    n_o2 = n * 0.20946
    n_n2 = n * 0.78084

    # 1 -> 2
    t2n2 = ((v1 - n_n2 * b) / (v2 - n_n2 * b)) ** (2.0 / 5.0) * t1
    t2o2 = ((v1 - n_o2 * bo2) / (v2 - n_o2 * bo2)) ** (2.0 / 5.0) * t1

    # temperature lost due to cooling from Newton's law of cooling
    t2 = (t2n2 * n_n2 + t2o2 * n_o2) / (n_n2 + n_o2)

    p2n2 = n_n2 * R * t2 / (v2 - n_n2 * b) - a * n_n2 ** 2 / (v2 ** 2)
    p2o2 = n_o2 * R * t2 / (v2 - n_o2 * bo2) - ao2 * n_o2 ** 2 / (v2 ** 2)
    p2 = p2n2 + p2o2

    # 2 -> 3
    # knowing that 14.7:1 is the ideal ratio w.r.t. mass.
    air_fuel_ratio_moles = OPTIMAL_AIR_FUEL_RATIO / AIR_MOLAR_MASS * OCTANE_MOLAR_MASS
    num_moles_octane = n / air_fuel_ratio_moles
    num_moles_co2 = num_moles_octane * 9
    num_moles_n2 = n * 0.78
    total_moles_left = num_moles_n2 + num_moles_co2

    # using heat of combustion we find excess energy released in combustion
    excess_energy = num_moles_octane * ENTHALPY_COMBUSTION_OCTANE
    deltaE_n2 = excess_energy * num_moles_n2 / total_moles_left
    deltaE_co2 = excess_energy * num_moles_co2 / total_moles_left

    t3n2 = deltaE_n2 / (num_moles_n2 * CV_N2) + t2
    t3co2 = deltaE_co2 / (num_moles_co2 * CV_CO2) + t2

    t3 = (t3n2 * num_moles_n2 + t3co2 * num_moles_co2) / (num_moles_n2 + num_moles_co2)

    # 3 -> 4
    v_compression = engine_kinematics(bore, stroke, con_rod, cr, 180, 0)
    constant = p1 * pow(v1, gamma)

    p_compression = []
    for v in v_compression:
        p_compression.append(constant / pow(v, gamma))
    # state3
    v3 = v2

    # p3v3/t3=p2v2/t2  |   p2v2/t2=rhs | t2=p2v2/rhs
    rhs = p2 * v2 / t2
    p3 = rhs * t3 / v3

    v_expansion = engine_kinematics(bore, stroke, con_rod, cr, 0, 180)

    constant = p3 * pow(v3, gamma)

    p_expansion = []
    for v in v_expansion:
        p_expansion.append(constant / pow(v, gamma))

    # 4 -> 1
    v4 = v1
    t4n = t3 * ((v3 - num_moles_n2 * b) / (v4 - num_moles_n2 * b)) ** (2.0 / 5.0)
    t4c = t3 * ((v3 - num_moles_co2 * bco2) / (v4 - num_moles_co2 * bco2)) ** (2.0 / 5.0)
    t4 = (t4n * num_moles_n2 + t4c * num_moles_co2) / (num_moles_n2 + num_moles_co2)

    mean_temp = ((t1 + t2)/2 * 0.5 + (t3+t3)/2 * 0.5 + (t4+t1)/2 * 0.5 + t1*0.5)/2
    q = ENGINE_K * sa * (mean_temp-cool_temp) / ENGINE_THICKNESS
    heat_transfer = q * time_diff * 2
    engine_temp_f = t1 + heat_transfer / (ENGINE_HEAT_CAPACITY * ENGINE_MASS)

    p4n2 = num_moles_n2 * R * t4 / (v4 - num_moles_n2 * b) - a * num_moles_n2 * 2 / (v4 ** 2)
    p4co2 = num_moles_co2 * R * t4 / (v4 - num_moles_co2 * bco2) - aco2 * num_moles_co2 ** 2 / (v4 ** 2)
    p4 = p4n2 + p4co2

    """
    # plot
    plt.plot([v2, v3], [p2, p3])
    plt.plot(v_compression, p_compression)
    plt.plot(v_expansion, p_expansion)
    plt.plot([v4, v1], [p4, p1])
    plt.xlabel('Volume (cubic meters)')
    plt.ylabel('Pressure (Pa)')
    plt.title('P-V Diagram for van der Waal gas Otto cycle')
    plt.show()
    """
    return [engine_temp_f, p1]


temperature = [323]
pressure = [101325]
data = pd.read_csv("track_data2.csv")
time = data['Time (s)']
rpm = data['RPM']
coolant_temps = data['Coolant']
spds = data['Speed']
oil_temps = data['Oil']

for i in range(time.size - 1):
    dt = time[i+1] - time[i]
    rps = rpm[i]/60.0
    v = spds[i]
    ct = coolant_temps[i] + 273
    num_four_strokes = rps * dt / 2
    for j in range(num_four_strokes.astype(int)):
        results = four_stroke(pressure[-1], temperature[-1], rps, v, ct)
        # print(results[0])
        pressure[-1] = results[1]
        temperature[-1] = results[0]
    results = four_stroke(pressure[-1], temperature[-1], rps, v, ct)
    temperature.append(results[0])
    pressure.append(results[0])

for i in range(len(temperature)):
    temperature[i] -= 273
plt.plot(time, temperature, label='Simulated Temperature')
plt.plot(time, coolant_temps, label='Coolant Temperature')
plt.plot(time, oil_temps, label='Oil Temperature')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (C)")
plt.title("Engine Temperature Over Time")
plt.legend()
plt.show()







