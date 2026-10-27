from math import pi, sqrt

# INPUTS: edit only the values in this section.
AREA_MODE = 'CAD'                 # CAD = measured STEP restrictions; FGSE = worksheet design orifices
MODEL_FRACTION = 0.5              # Fraction of complete injector represented in Fluent
CD_ETHANOL = 0.75                 # FGSE assumed discharge coefficient, not experimentally calibrated
CD_N2O = 0.75                     # FGSE assumed discharge coefficient, not experimentally calibrated
P_ETHANOL_PA = 5_171_070.0        # FGSE absolute pressure upstream of injector [Pa]
P_N2O_PA = 5_171_070.0            # FGSE absolute pressure upstream of injector [Pa]
P_CHAMBER_PA = 1_723_690.0        # FGSE absolute chamber pressure [Pa]
P_NOZZLE_EXIT_PA = 101_325.0      # Absolute exit pressure [Pa]; display only, not used for injector flow
P_OPERATING_PA = 101_325.0        # Fluent operating pressure [Pa]
RHO_ETHANOL = 785.2               # FGSE fuel density [kg/m^3]
RHO_N2O = 744.3344972             # FGSE N2O density [kg/m^3], liquid-model assumption
ETHANOL_HOLE_COUNT = 8            # Full-model equivalent number of radial holes in STEP
ETHANOL_HOLE_DIAMETER_MM = 0.92   # STEP radial hole diameter [mm]
N2O_ANNULUS_INNER_R_MM = 8.000    # STEP narrow annulus inner radius [mm]
N2O_ANNULUS_OUTER_R_MM = 8.305    # STEP narrow annulus outer radius [mm]
ETHANOL_FACE_DIAMETER_MM = 8.0    # Full-model equivalent upstream inlet diameter [mm]
N2O_FACE_DIAMETER_MM = 5.0        # Full-model equivalent upstream inlet diameter [mm]
NOZZLE_EXIT_RADIUS_MM = 29.08793 # Full-model nozzle exit radius [mm]
FGSE_ETHANOL_AREA_MM2 = 9.786473082  # FGSE full-engine design orifice area [mm^2]
FGSE_N2O_AREA_MM2 = 30.15459725     # FGSE full-engine design orifice area [mm^2]
FGSE_ETHANOL_TARGET_KGS = 0.5400540054  # FGSE target full-engine fuel flow [kg/s]
FGSE_N2O_TARGET_KGS = 1.620162016      # FGSE target full-engine oxidizer flow [kg/s]


def circle_area(diameter_mm):
    return pi * diameter_mm ** 2 / 4


def orifice_mdot(cd, area_mm2, rho, p_in):
    dp = p_in - P_CHAMBER_PA
    if not (0 < cd <= 1 and area_mm2 > 0 and rho > 0 and dp > 0):
        raise ValueError('Require Cd in (0,1], positive area/density, and injector pressure > chamber pressure.')
    return cd * area_mm2 * 1e-6 * sqrt(2 * rho * dp)


def solve():
    if AREA_MODE not in ('CAD', 'FGSE') or not 0 < MODEL_FRACTION <= 1:
        raise ValueError('AREA_MODE must be CAD or FGSE, and MODEL_FRACTION must be in (0,1].')
    if N2O_ANNULUS_OUTER_R_MM <= N2O_ANNULUS_INNER_R_MM:
        raise ValueError('Annulus outer radius must exceed its inner radius.')
    cad_ethanol_area = ETHANOL_HOLE_COUNT * circle_area(ETHANOL_HOLE_DIAMETER_MM)
    cad_n2o_area = pi * (N2O_ANNULUS_OUTER_R_MM**2 - N2O_ANNULUS_INNER_R_MM**2)
    if AREA_MODE == 'CAD':
        ethanol_area, n2o_area = cad_ethanol_area, cad_n2o_area
    else:
        ethanol_area, n2o_area = FGSE_ETHANOL_AREA_MM2, FGSE_N2O_AREA_MM2
    ethanol_full = orifice_mdot(CD_ETHANOL, ethanol_area, RHO_ETHANOL, P_ETHANOL_PA)
    n2o_full = orifice_mdot(CD_N2O, n2o_area, RHO_N2O, P_N2O_PA)
    ethanol_half, n2o_half = ethanol_full * MODEL_FRACTION, n2o_full * MODEL_FRACTION
    total = ethanol_half + n2o_half
    area_ethanol_face = MODEL_FRACTION * circle_area(ETHANOL_FACE_DIAMETER_MM) * 1e-6
    area_n2o_face = MODEL_FRACTION * circle_area(N2O_FACE_DIAMETER_MM) * 1e-6
    area_exit = MODEL_FRACTION * pi * (NOZZLE_EXIT_RADIUS_MM * 1e-3)**2
    v_ethanol = ethanol_half / (RHO_ETHANOL * area_ethanol_face)
    v_n2o = n2o_half / (RHO_N2O * area_n2o_face)
    v_exit = (ethanol_half / RHO_ETHANOL + n2o_half / RHO_N2O) / area_exit
    return {
        'ethanol_half_kg_s': ethanol_half,
        'n2o_half_kg_s': n2o_half,
        'outlet_half_kg_s': total,
        'ethanol_inlet_velocity_m_s': v_ethanol,
        'n2o_inlet_velocity_m_s': v_n2o,
        'outlet_bulk_velocity_m_s': v_exit,
        'ethanol_fraction': ethanol_half / total,
        'oxidizer_fuel_ratio': n2o_half / ethanol_half,
        'ethanol_pct_vs_fgse': 100 * (ethanol_full / FGSE_ETHANOL_TARGET_KGS - 1),
        'n2o_pct_vs_fgse': 100 * (n2o_full / FGSE_N2O_TARGET_KGS - 1),
    }


if __name__ == '__main__':
    r = solve()
    print(f'Area mode: {AREA_MODE} | injector ΔP: {(P_ETHANOL_PA - P_CHAMBER_PA)/1e6:.6f} MPa (ethanol); {(P_N2O_PA - P_CHAMBER_PA)/1e6:.6f} MPa (N2O)')
    print(f'Ethanol:  {r["ethanol_half_kg_s"]:.6f} kg/s half | inlet {r["ethanol_inlet_velocity_m_s"]:.3f} m/s | vs FGSE {r["ethanol_pct_vs_fgse"]:+.1f}%')
    print(f'N2O:      {r["n2o_half_kg_s"]:.6f} kg/s half | inlet {r["n2o_inlet_velocity_m_s"]:.3f} m/s | vs FGSE {r["n2o_pct_vs_fgse"]:+.1f}%')
    print(f'Outlet:   {r["outlet_half_kg_s"]:.6f} kg/s net | mean liquid-only speed {r["outlet_bulk_velocity_m_s"]:.3f} m/s')
    print(f'Mixture:  ethanol Y={r["ethanol_fraction"]:.4f} | N2O/ethanol={r["oxidizer_fuel_ratio"]:.3f}')
    print(f'Fluent:   inlet gauge reference={P_ETHANOL_PA-P_OPERATING_PA:.0f} Pa ethanol; {P_N2O_PA-P_OPERATING_PA:.0f} Pa N2O | outlet={P_NOZZLE_EXIT_PA-P_OPERATING_PA:.0f} Pa gauge')
    print('Limit: liquid-orifice screening only; not validated for N2O flashing, injector losses, or nozzle gas dynamics.')
