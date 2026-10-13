"""
n2o_tank_warmup.py - 1D
    air + sun  -->  bare Al 6061-T6 tank wall  -->  N2O
"""
import numpy as np
import CoolProp.CoolProp as CP
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt

SIGMA = 5.670e-8                # W/m2-K4, Stefan-Boltzmann constant
N2O = "NitrousOxide"  

#  GIVEN 
D_o    = 6.00 * 0.0254          # m, tank OD (6 in)
t_w    = 0.25 * 0.0254          # m, wall (0.25 in)
rho_Al, c_Al = 2700.0, 896.0    # Al 6061-T6: kg/m3, J/kg-K (MatWeb)
T0, T_goal = 283.15, 288.0      # K, load and launch temperatures

#  DTEG 2026 App. G / 6.3.2 
T_AIR = {"307.15 K air": 307.15, "293.65 K air": 293.65}   # avg high (33-35 C), avg low (19-22 C): midpoints
wind  = 5.2                     # m/s, average wind
P_atm = 90.9e3                  # Pa, air pressure at 890 m site elevation

#  ASSUMED
fill  = 0.85                    # liquid volume fraction at T0
q_sun = 320.0                   # W/m2, clear-sky sun averaged around a vertical cylinder
FINISH = {"bare Al": (0.35, 0.09)}   # no paint, mill finish: solar absorptivity alpha, IR emissivity eps

D_i = D_o - 2 * t_w
perim = np.pi * D_o
m_Al = rho_Al * np.pi / 4 * (D_o**2 - D_i**2)

rho_l0 = CP.PropsSI("D", "T", T0, "Q", 0, N2O)            # saturated liquid density at T0
rho_v0 = CP.PropsSI("D", "T", T0, "Q", 1, N2O)            # saturated vapor density at T0
rho_bar = fill * rho_l0 + (1 - fill) * rho_v0
m_N2O = rho_bar * np.pi / 4 * D_i**2

# CoolProp (T, rho)
def u_N2O(T):
    return CP.PropsSI("U", "T", T, "D", rho_bar, N2O)    # J/kg

def C_per_m(T, dT=0.05):
    return m_N2O * (u_N2O(T + dT) - u_N2O(T - dT)) / (2 * dT) + m_Al * c_Al

# (283.15 K -> 288 K)
C_closed_form = m_N2O * (u_N2O(T_goal) - u_N2O(T0)) / (T_goal - T0) + m_Al * c_Al

T_full = brentq(lambda T: CP.PropsSI("D", "T", T, "Q", 0, N2O) - rho_bar, T0, 309.0)

def h_conv(T_s, T_air):
    Tf = 0.5 * (T_s + T_air)
    k  = CP.PropsSI("L", "T", Tf, "P", P_atm, "Air")                                       # W/m-K
    nu = CP.PropsSI("V", "T", Tf, "P", P_atm, "Air") / CP.PropsSI("D", "T", Tf, "P", P_atm, "Air")  # m2/s
    Pr = CP.PropsSI("Prandtl", "T", Tf, "P", P_atm, "Air")
    Re = wind * D_o / nu
  
    # Nusselt number (Churchill-Bernstein, crossflow over a cylinder):
    Nu = 0.3 + 0.62 * Re**0.5 * Pr**(1/3) / (1 + (0.4 / Pr)**(2/3))**0.25 \
         * (1 + (Re / 282000)**(5/8))**(4/5)
    return Nu * k / D_o

def air_side(T, T_air, alpha, eps):
    h_c = h_conv(T, T_air)
    h_r = eps * SIGMA * (T**2 + T_air**2) * (T + T_air)
    T_sa = T_air + alpha * q_sun / (h_c + h_r)
    UA = perim * (h_c + h_r)
    return UA, T_sa, h_c, h_r

# Time to 288 K, T at 15 min, load temperature
if __name__ == "__main__":
    print(f"per meter: m_Al = {m_Al:.3f} kg/m, m_N2O = {m_N2O:.3f} kg/m (rho_bar {rho_bar:.1f} kg/m3), "
          f"c_eff = {(C_closed_form - m_Al*c_Al)/m_N2O:.0f} J/kg-K, C' = {C_closed_form:.0f} J/m-K, "
          f"liquid-full at {T_full:.1f} K\n")
    cols = ["Tank, Air Temp", "h_c", "h_r", "T_sa", "UA'", "tau", "t288 est", "t288 ODE", "T15", "T30", "T0 for 288@15"]
    w = [max(len(c), 7) for c in cols]; w[0] = 22; w[9] = 14
    hdr = "  ".join(f"{c:<{w[0]}s}" if i == 0 else f"{c:>{w[i]}s}" for i, c in enumerate(cols))
    print(hdr + "\n" + "-" * len(hdr))
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for air_name, T_air in T_AIR.items():
        for fin, (alpha, eps) in FINISH.items():
            UA, Tsa, hc, hr = air_side(0.5 * (T0 + T_goal), T_air, alpha, eps)
            tau = C_closed_form / UA
            t_est = tau * np.log((Tsa - T0) / (Tsa - T_goal)) / 60
            T0_req = Tsa - (Tsa - T_goal) * np.exp(15 * 60 / tau)

            def rhs(t, y):
                UA_, Tsa_, _, _ = air_side(y[0], T_air, alpha, eps)
                return [UA_ * (Tsa_ - y[0]) / C_per_m(y[0])]
            hit = lambda t, y: y[0] - T_goal                   # event: T = 288 K
            full = lambda t, y: y[0] - (T_full - 1e-3)         # event: liquid-full, stop
            full.terminal = True
            s = solve_ivp(rhs, [0, 1800], [T0], events=[hit, full], dense_output=True, max_step=10, rtol=1e-8)
            t_ode = s.t_events[0][0] / 60
            t_stop = s.t[-1]
            T15 = s.sol(900)[0]
            T30 = f"{s.sol(1800)[0]:.1f}" if t_stop >= 1800 else f"full @ {t_stop/60:.0f} min"
            name = f"{fin}, {air_name}"
            vals = [f"{hc:.1f}", f"{hr:.1f}", f"{Tsa:.1f}", f"{UA:.2f}", f"{tau/60:.1f}", f"{t_est:.1f}",
                    f"{t_ode:.1f}", f"{T15:.1f}", T30, f"{T0_req:.1f}"]
            print(f"{name:<{w[0]}s}  " + "  ".join(f"{v:>{w[i+1]}s}" for i, v in enumerate(vals)))
            tt = np.linspace(0, t_stop, 181)
            ax.plot(tt / 60, s.sol(tt)[0], label=name)

    print("\nunits: h W/m2-K | T K | UA' W/m-K | tau, t288 min")
    ax.axhline(T_goal, c="k", lw=0.8); ax.axvline(15, c="gray", ls=":")
    ax.axhline(T_full, c="r", lw=0.8, ls=":"); ax.text(0.3, T_full + 0.3, "liquid-full", fontsize=8)
    ax.set_xlabel("time after fill [min]"); ax.set_ylabel("N2O temperature [K]")
    ax.set_title("Exposed 6 in x 0.25 Al tank, IREC sun, 5.2 m/s wind", fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig("n2o_tank_warmup.png", dpi=150)
    print("plot saved: n2o_tank_warmup.png")
