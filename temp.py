import numpy as np
import math
import time

# Define the values of constants
DENSITY_MUSCLE = 1140                       # kg/m^3, density of muscle
SPEC_HEAT_CAPA_MUSCLE = 3570                # J/kg.K, specific heat capacity of muscle
CONS_A = DENSITY_MUSCLE * SPEC_HEAT_CAPA_MUSCLE
THERMAL_CONDUCT = 0.498                     # W/m.K, thermal conductivity of muscle
CONS_B = THERMAL_CONDUCT
PERF_RATE_MUSCLE = 2 * 3.97 * pow(10, -4)   # /s, perfusion rate of blood
SPEC_HEAT_CAPA_BLOOD = 3850                 # J/kg.K, specific heat capacity of blood
CONS_C = PERF_RATE_MUSCLE * SPEC_HEAT_CAPA_BLOOD

WAV_LEN_X = 10 * pow(10, -3)                # m, wave length in direction of x
WAV_LEN_Y = 10 * pow(10, -3)                # m, ...y
WAV_LEN_Z = 30 * pow(10, -3)                # m, ...z
STEP_LEN_X = 0.5 * pow(10, -3)              # m, step length in direction of x
STEP_LEN_Y = 0.5 * pow(10, -3)              # m, ...y
STEP_LEN_Z = 0.5 * pow(10, -3)              # m, ...z
X_SCALE = int(WAV_LEN_X / STEP_LEN_X + 1)   # The total index of list x
Y_SCALE = int(WAV_LEN_Y / STEP_LEN_Y + 1)   # ...y
Z_SCALE = int(WAV_LEN_Z / STEP_LEN_Z + 1)   # ...z

HEAT_TIME = 20                              # Duration of heating period
TOTAL_TIME = 60                             # Duration of total period
DELTA_T = 0.5                                 # s, the interval of time
TIME_INDEX = int(TOTAL_TIME / DELTA_T + 1)  # Total index of time list

PI = math.pi                                # 3.14 value of Pi

INIT_TEMP = 37                              # Initial tissue temperature


# Data initialization
t_star_0 = basicTemp_fft
q_star_0 = Q_fft
#q_time = inv_fft_func(Q_fft)

t_time = [np.zeros((Z_SCALE, Y_SCALE, X_SCALE))] * TIME_INDEX
t_time[0] = np.ones((Z_SCALE, Y_SCALE, X_SCALE)) * INIT_TEMP

t_star_freq = [np.zeros((Z_SCALE, Y_SCALE, X_SCALE))] * TIME_INDEX
t_star_freq[0] = t_star_0


# Function for FFT
def fft_func(in_time):

    out_freq = abs(np.fft.fft2(in_time))

    return out_freq


# Function for inverse FFT
def inv_fft_func(in_freq):

    out_time = abs(np.fft.ifft2(in_freq))

    return out_time


# Function for creation of the array containing nu coef
def func_nu_square(k, j, i):

    delta_x = 1 / WAV_LEN_X / (WAV_LEN_X / STEP_LEN_X)
    delta_y = 1 / WAV_LEN_Y / (WAV_LEN_Y / STEP_LEN_Y)
    delta_z = 1 / WAV_LEN_Z / (WAV_LEN_Z / STEP_LEN_Z)

    m = i - WAV_LEN_X / STEP_LEN_X / 2
    n = j - WAV_LEN_Y / STEP_LEN_Y / 2
    p = k - WAV_LEN_Z / STEP_LEN_Z / 2

    nv_square_value = pow(m, 2)*pow(delta_x, 2) + pow(n, 2)*pow(delta_y, 2) + pow(p, 2)*pow(delta_z, 2)

    return nv_square_value


# Solution of equation while heating
def solution_one_loop_soni(t_star_in, q_star, nu_square):

    coef_exp = -4 * pow(PI, 2) * CONS_B * DELTA_T / CONS_A                  # Coefficient, to simplify calculation
    coef_deno = 4 * pow(PI, 2) * CONS_B                                     # Coefficient, ...
    array1 = coef_exp * nu_square
    array_exp = np.exp(array1)                                              # Value of exp member

    t_star_out_coef1 = t_star_in * array_exp                                # First coefficient of t_star_out

    t_star_out_coef21 = q_star / (coef_deno * nu_square)
    t_star_out_coef22 = array_exp - 1
    t_star_out_coef2 = t_star_out_coef21 * t_star_out_coef22                # Second coefficient of t_star_out

    t_star_out = t_star_out_coef1 - t_star_out_coef2

    return t_star_out


# Solution of equation while cooling
def solution_one_loop_cool(t_star_in, nu_square):

    coef_exp = -4 * pow(PI, 2) * CONS_B * DELTA_T / CONS_A                  # Coefficient, to simplify calculation
    coef_deno = 4 * pow(PI, 2) * CONS_B                                     # Coefficient, ...
    array1 = coef_exp * nu_square
    array_exp = np.exp(array1)                                              # Value of exp member

    t_star_out_coef1 = t_star_in * array_exp                                # First coefficient of t_star_out

    t_star_out_coef21 = 0 / (coef_deno * nu_square)                         # Q* is 0 during cooling period
    t_star_out_coef22 = array_exp - 1
    t_star_out_coef2 = t_star_out_coef21 * t_star_out_coef22                # Second coefficient of t_star_out

    t_star_out = t_star_out_coef1 - t_star_out_coef2

    return t_star_out


# Calculate the array of T by loop
def solution_iter():

    for i in range(int(HEAT_TIME / DELTA_T)):
        t_star_out = solution_one_loop_soni(t_star_freq[i], q_star_0, array_nu_square)
        t_star_out[int((Z_SCALE - 1) / 2), int((Y_SCALE - 1) / 2), int((X_SCALE - 1) / 2)] = 0
        t_star_freq[i+1] = t_star_out
        t_temp = inv_fft_func(t_star_out)
        t_temp = t_temp.reshape((Z_SCALE * Y_SCALE * X_SCALE, 1))
        for x in range(t_temp.size):
            if t_temp[x] < INIT_TEMP:
                t_temp[x] = INIT_TEMP
        t_temp = t_temp.reshape((Z_SCALE, Y_SCALE, X_SCALE))
        t_time[i+1] = t_temp

    for j in range(int(HEAT_TIME / DELTA_T), int(TOTAL_TIME / DELTA_T)):
        t_star_out = solution_one_loop_cool(t_star_freq[j], array_nu_square)
        t_star_out[(int(Z_SCALE - 1) / 2), int((Y_SCALE - 1) / 2), int((X_SCALE - 1) / 2)] = 0
        t_star_freq[j+1] = t_star_out
        t_temp = inv_fft_func(t_star_out)
        t_temp = t_temp.reshape((Z_SCALE * Y_SCALE * X_SCALE, 1))
        for x in range(t_temp.size):
            if t_temp[x] < INIT_TEMP:
                t_temp[x] = INIT_TEMP
        t_temp = t_temp.reshape((Z_SCALE, Y_SCALE, X_SCALE))
        t_time[j+1] = t_temp


    return 0


start1 = time.time()
array_nu_square = np.fromfunction(func_nu_square, (Z_SCALE, Y_SCALE, X_SCALE))
stop1 = time.time()
print("The generation of coefficient of nu cost: " + str(stop1 - start1))

# start2 = time.time()
# t_star_1 = solution_one_loop(t_star_0, q_star_0, array_nu_square)
# t_time[1] = inv_fft_func(t_star_1)
# t_star_freq[1] = t_star_1
# stop2 = time.time()
# print("The generation of coefficient of one loop: " + str(stop2 - start2))

start2 = time.time()
solution_iter()
stop2 = time.time()
print("The calculation time of T costs: " + str(stop2 - start2))

"""
a = t_star_freq[0]
b = t_star_freq[1][31,2,1]
at = np.fft.ifft2(a)
ct = np.fft.ifft(b)
"""